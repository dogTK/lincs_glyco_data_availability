"""Which glycan-related programs do mechanism-of-action classes move in HepG2?

The compound-level values are SigCom LINCS Level 5 signatures (characteristic direction coefficients of the
CycleGAN-predicted RNA-seq-like profiles); "signature value" below means these coefficients.

1. Annotate the 4,704 HepG2 compounds with Drug Repurposing Hub MOA/target (BRD -> InChIKey14 -> name)
2. MOA class x module: module eigengenes of each class (>= 5 compounds) vs all other compounds
3. Three programs defined by GlycoEnzOnto pathway membership (lysosomal glycosaminoglycan degradation,
   nucleotide sugar synthesis, terminal sialylation): program score shift of every MOA class
4. HDAC inhibitors: module profile, perturbation-strength-matched control of the lysosomal score,
   per-condition percentile, and gene-level response
5. Figure 2 (ATC class x module GSEA) and Figure 3A, 3B, 3D
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'results/20_moa_glyco_programs'
OUT.mkdir(parents=True, exist_ok=True)
ME_CSV = ROOT / 'results/02_wgcna/module_eigengenes.csv'
ATC_GSEA = ROOT / 'results/03_atc_pathway_gsea/atc_L1_module_gsea.csv'
HUB_DRUGS = ROOT / 'data/external/repurposing_drugs_20200324.txt'
HUB_SAMPLES = ROOT / 'data/external/repurposing_samples_20200324.txt'
WIDE_PQ = ROOT / 'data/processed/hepg2_glyco_genes_wide.parquet'
FIGDIR = ROOT / 'results/figures'
META = ['VALUE', 'canonical_smiles', 'cell', 'cmapid', 'compound_alias', 'dose', 'inchi_key', 'pertid', 'pertname',
        'timepoint', 'sample_id']

N_PERM = 10_000
MIN_SET = 5
RNG = np.random.default_rng(42)

HYDROXAMATE = ['PCI-24781', 'PCI-34051', 'vorinostat', 'trichostatin-a', 'scriptaid', 'dacinostat',
               'JNJ-26481585', 'pyroxamide', 'belinostat', 'tubastatin-a']
HX_IDS = HYDROXAMATE + ['BRD-K68202742', 'BRD-A19037878']  # trichostatin A and vorinostat under Broad IDs
BENZAMIDE = ['mocetinostat', 'entinostat', 'TC-H-106', 'Merck60']
# GlycoEnzOnto pathway members detected in HepG2 (Supplementary Table S1)
SIDE_GENE_SETS = {
    'lysosomal_degradation': ['ARSB', 'ARSK', 'GALNS', 'GLB1', 'GNS', 'GUSB', 'HEXA', 'HEXB', 'HYAL1', 'IDS', 'IDUA',
                              'NEU1', 'NEU3', 'NEU4'],
    'nucleotide_sugar': ['ALG5', 'AMDHD2', 'CMAS', 'DPM1', 'DPM2', 'DPM3', 'FPGT', 'GALE', 'GALK1', 'GALK2', 'GALT',
                         'GCK', 'GFPT1', 'GFPT2', 'GMDS', 'GMPPA', 'GMPPB', 'GNE', 'GNPDA1', 'GNPDA2', 'GNPNAT1', 'GPI',
                         'HK1', 'HK2', 'HK3', 'MPI', 'NAGK', 'NANS', 'PAPSS1', 'PAPSS2', 'PGM1', 'PGM2', 'PGM3', 'PMM1',
                         'PMM2', 'RENBP', 'UAP1', 'UGDH', 'UGP2', 'UXS1'],
    'sialylation': ['ST3GAL1', 'ST3GAL2', 'ST3GAL3', 'ST3GAL4', 'ST3GAL5', 'ST3GAL6', 'ST6GAL1', 'ST6GAL2',
                    'ST6GALNAC1', 'ST6GALNAC2', 'ST6GALNAC3', 'ST6GALNAC4', 'ST6GALNAC5', 'ST6GALNAC6', 'ST8SIA1',
                    'ST8SIA2', 'ST8SIA3', 'ST8SIA4', 'ST8SIA5', 'ST8SIA6'],
}
# lysosomal glycan-degrading enzymes shown at the gene level: the GAG degradation program plus the
# N-linked glycan lysosomal degradation module (Module 3, Table 1)
LYSOSOMAL_PANEL = SIDE_GENE_SETS['lysosomal_degradation'] + ['AGA', 'FUCA1', 'MAN2B1', 'MANBA']
MODULE_LABELS = {'ME1': 'M1 ER quality control', 'ME2': 'M2 Glucuronidation', 'ME3': 'M3 N-glycan lysosomal degr.',
                 'ME5': 'M5 Dolichol precursor', 'ME6': 'M6 HA degradation', 'ME7': 'M7 UDP-GlcNAc',
                 'ME9': 'M9 KS/GAG lysosomal degr.', 'ME10': 'M10 UDP-glucose', 'ME14': 'M14 Ganglioside core',
                 'ME18': 'M18 GPI anchor'}


def savefig(fig, name):
    """成果物ごとのフォルダに PNG/PDF/SVG を保存する"""
    outdir = FIGDIR / name
    outdir.mkdir(parents=True, exist_ok=True)
    for ext in ('png', 'pdf', 'svg'):
        fig.savefig(outdir / f'{name}.{ext}', dpi=600, bbox_inches='tight')


def rbc(a, b):
    u, p = mannwhitneyu(a, b)
    return 2 * u / (len(a) * len(b)) - 1, p


# ---------------------------------------------------------------- 1. annotation
def annotate():
    me = pd.read_csv(ME_CSV).rename(columns={'pertname': 'compound_id'})
    ik = pd.read_parquet(WIDE_PQ, columns=['pertname', 'inchi_key']).dropna().drop_duplicates('pertname')
    ph = me[['compound_id']].merge(ik.rename(columns={'pertname': 'compound_id'}), how='left')
    dr = pd.read_csv(HUB_DRUGS, sep='\t', comment='!')
    sa = pd.read_csv(HUB_SAMPLES, sep='\t', comment='!')
    sa['brd'] = sa['broad_id'].str[:13]
    sa['ik14'] = sa['InChIKey'].str[:14]
    cols = ['pert_iname', 'moa', 'target', 'clinical_phase']
    ann = sa[['pert_iname', 'brd', 'ik14']].drop_duplicates().merge(dr[cols], on='pert_iname')

    by_brd = ann.drop_duplicates('brd').set_index('brd')[cols]
    by_ik = ann.dropna(subset=['ik14']).drop_duplicates('ik14').set_index('ik14')[cols]
    by_name = dr.assign(k=dr['pert_iname'].str.lower()).drop_duplicates('k').set_index('k')[cols]

    hit_brd = by_brd.reindex(ph['compound_id']).reset_index(drop=True)
    hit_ik = by_ik.reindex(ph['inchi_key'].str[:14]).reset_index(drop=True)
    hit_name = by_name.reindex(ph['compound_id'].str.lower()).reset_index(drop=True)
    hub = hit_brd.combine_first(hit_ik).combine_first(hit_name)
    hub['match'] = np.select([hit_brd['pert_iname'].notna(), hit_ik['pert_iname'].notna(),
                              hit_name['pert_iname'].notna()], ['brd', 'inchikey14', 'name'], None)
    ph = pd.concat([ph.reset_index(drop=True), hub], axis=1).merge(me, on='compound_id')
    ph.to_csv(OUT / 'compounds_with_hub.csv', index=False)
    print(f'Hub MOA coverage: {ph["moa"].notna().mean():.1%} ({ph["moa"].notna().sum()} compounds); match:',
          ph['match'].value_counts().to_dict())
    return ph


def compound_values():
    """compound-level median signature value over all conditions, and the sample-level table"""
    wide = pd.read_parquet(WIDE_PQ)
    genes = [c for c in wide.columns if c not in META]
    z = wide[genes].apply(pd.to_numeric, errors='coerce')
    return z.groupby(wide['pertname']).median(), wide[['pertname', 'dose', 'timepoint']].join(z)


def explode_moa(ph):
    x = ph.assign(category=ph['moa'].str.split('|')).explode('category').dropna(subset=['category'])
    return {c: set(g['compound_id']) for c, g in x.groupby('category') if g['compound_id'].nunique() >= MIN_SET}


# ---------------------------------------------------------------- 2. MOA class x module
def moa_modules(ph, classes):
    me_cols = [c for c in ph.columns if c.startswith('ME')]
    rows = []
    for cat, ids in classes.items():
        g, rest = ph[ph['compound_id'].isin(ids)], ph[~ph['compound_id'].isin(ids)]
        for m in me_cols:
            r, p = rbc(g[m], rest[m])
            rows.append(dict(moa=cat, n=len(g), module=m, r_rb=r, p_value=p))
    res = pd.DataFrame(rows)
    res['q_value'] = res.groupby('module')['p_value'].transform(lambda p: multipletests(p, method='fdr_bh')[1])
    res.to_csv(OUT / 'moa_module_shift.csv', index=False)
    print(f'\nMOA class x module ({len(classes)} classes), q < 0.05:')
    print(res[res['q_value'] < 0.05].sort_values(['module', 'r_rb']).round(4).to_string(index=False))
    return res


# ---------------------------------------------------------------- 3. MOA class x program
def moa_side_shift(ph, cz, classes):
    scores = pd.DataFrame({side: cz[[g for g in genes if g in cz.columns]].mean(1)
                           for side, genes in SIDE_GENE_SETS.items()})
    d = ph[['compound_id']].merge(scores, left_on='compound_id', right_index=True)
    rows = []
    for cat, ids in classes.items():
        g, rest = d[d['compound_id'].isin(ids)], d[~d['compound_id'].isin(ids)]
        for side in SIDE_GENE_SETS:
            r, p = rbc(g[side], rest[side])
            rows.append(dict(moa=cat, n=len(g), side=side, r_rb=r, p_value=p))
    res = pd.DataFrame(rows)
    res['q_value'] = res.groupby('side')['p_value'].transform(lambda p: multipletests(p, method='fdr_bh')[1])
    res.to_csv(OUT / 'moa_side_shift.csv', index=False)
    print(f'\nMOA classes: shift of program scores vs other compounds, q < 0.05')
    for side in SIDE_GENE_SETS:
        s = res[(res['side'] == side)]
        print(f'  {side}: min q = {s["q_value"].min():.3f}, nominal p < 0.05 in {(s["p_value"] < 0.05).sum()} classes')
        print(s[s['q_value'] < 0.05][['moa', 'n', 'r_rb', 'q_value']].round(4).to_string(index=False))
    return res, scores


# ---------------------------------------------------------------- 4. HDAC inhibitors
def hdaci(ph, cz, samples, scores):
    me_cols = [c for c in ph.columns if c.startswith('ME')]
    is_h = ph['moa'].fillna('').str.contains('HDAC inhibitor').values

    # module profile
    rows = []
    for c in me_cols:
        r, p = rbc(ph.loc[is_h, c], ph.loc[~is_h, c])
        rows.append(dict(module=c, hdaci_median=ph.loc[is_h, c].median(), others_median=ph.loc[~is_h, c].median(),
                         r_rb=r, p_value=p))
    mp = pd.DataFrame(rows)
    mp['q_value'] = multipletests(mp['p_value'], method='fdr_bh')[1]
    mp.to_csv(OUT / 'hdaci_module_profile.csv', index=False)
    print(f'\nHDACi (n={is_h.sum()}) module eigengene profile vs other compounds')
    print(mp.round(4).to_string(index=False))

    # perturbation strength (sum |ME|) does not explain the lysosomal increase
    lys = scores['lysosomal_degradation'].reindex(ph['compound_id']).values
    sbin = pd.qcut(ph[me_cols].abs().sum(axis=1), 10, labels=False).values
    draws = np.column_stack([RNG.choice(lys[(sbin == b) & ~is_h], N_PERM) for b in sbin[is_h]])
    null = np.median(draws, axis=1)
    obs = np.median(lys[is_h])
    p = ((null >= obs).sum() + 1) / (N_PERM + 1)
    print(f'strength-matched permutation (lysosomal score): observed median = {obs:.5f}, '
          f'null median = {np.median(null):.5f}, p = {p:.4f}')

    # per condition: percentile of each HDACi sample's lysosomal score among all samples at the same time point
    lys_genes = SIDE_GENE_SETS['lysosomal_degradation']
    s = samples[['pertname', 'dose', 'timepoint']].assign(lysosomal_score=samples[lys_genes].mean(axis=1))
    s['pct_within_time'] = s.groupby('timepoint')['lysosomal_score'].rank(pct=True)
    hd = s[s['pertname'].isin(HX_IDS + BENZAMIDE)].copy()
    hd['chemotype'] = np.where(hd['pertname'].isin(BENZAMIDE), 'benzamide', 'hydroxamate')
    hd.sort_values(['chemotype', 'pertname', 'timepoint', 'dose']).to_csv(OUT / 'hdaci_per_condition_lysosomal.csv',
                                                                          index=False)
    print('\nHDACi per condition: lysosomal score percentile among all samples at the same time point')
    print(hd.groupby(['chemotype', 'timepoint']).agg(n=('pct_within_time', 'size'),
                                                     median_pct=('pct_within_time', 'median')).round(3).to_string())

    # gene level: hydroxamates vs compounds other than HDAC inhibitors
    hx = cz.loc[cz.index.isin(HX_IDS)]
    rest = cz.loc[~cz.index.isin(set(ph.loc[is_h, 'compound_id']) | set(HX_IDS))]
    rows = []
    for g in LYSOSOMAL_PANEL:
        r, p = rbc(hx[g], rest[g])
        rows.append(dict(gene=g, hydroxamate_median=hx[g].median(), others_median=rest[g].median(), r_rb=r, p_value=p))
    gl = pd.DataFrame(rows)
    gl['q_value'] = multipletests(gl['p_value'], method='fdr_bh')[1]
    gl.to_csv(OUT / 'hdaci_lysosomal_gene_level.csv', index=False)
    print('\nLysosomal enzymes: hydroxamate HDACi vs other compounds')
    print(gl.round(3).to_string(index=False))
    return gl


# ---------------------------------------------------------------- 5. figures
def style():
    plt.rcParams.update({'axes.labelsize': 16, 'axes.titlesize': 16, 'xtick.labelsize': 14, 'ytick.labelsize': 14,
                         'legend.fontsize': 14, 'font.size': 14, 'svg.fonttype': 'none'})


def heatmap(mat, q, rows, cols, name, cbar_label, vmax, figsize):
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(mat, cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='auto')
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            stars = '***' if q[i, j] < 0.001 else '**' if q[i, j] < 0.01 else '*' if q[i, j] < 0.05 else ''
            ax.text(j, i, stars, ha='center', va='center', fontsize=13,
                    color='white' if abs(mat[i, j]) > 0.6 * vmax else 'black')
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=0)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows)
    fig.colorbar(im, ax=ax, label=cbar_label, shrink=0.6)
    fig.tight_layout()
    savefig(fig, name)
    plt.close(fig)


def plot_atc():
    """Figure 2: GSEA of the ten modules for the 12 ATC Level 1 classes"""
    g = pd.read_csv(ATC_GSEA)
    nes = g.pivot(index='pathway', columns='atc_class', values='nes')
    q = g.pivot(index='pathway', columns='atc_class', values='q_value')
    order = [m for m in MODULE_LABELS if m in nes.index]
    cols = nes.loc[order].abs().sum().sort_values(ascending=False).index
    heatmap(nes.loc[order, cols].to_numpy(), q.loc[order, cols].to_numpy(), [MODULE_LABELS[m] for m in order], list(cols),
            'fig2_atc_module_gsea', 'NES', 2.5, (10, 6.5))


def plot_moa_modules(mm):
    """Figure 3A: rank-biserial r of MOA classes with q < 0.05 in at least one module"""
    keep = mm.loc[mm['q_value'] < 0.05, 'moa'].unique()
    r = mm[mm['moa'].isin(keep)].pivot(index='moa', columns='module', values='r_rb')
    q = mm[mm['moa'].isin(keep)].pivot(index='moa', columns='module', values='q_value')
    order = [m for m in MODULE_LABELS if m in r.columns]
    r = r[order].sort_values(['ME9', 'ME3'], ascending=False)
    n = mm.drop_duplicates('moa').set_index('moa')['n']
    heatmap(r.to_numpy(), q.loc[r.index, order].to_numpy(), [f'{c} (n={n[c]})' for c in r.index],
            [m.replace('ME', 'M') for m in order], 'fig3a_moa_module_shift', 'Rank-biserial r', 1, (11, 5))


def plot_programs(side):
    """Figure 3B: shift of each program score per MOA class (one dot per class); HDAC inhibitors in red"""
    programs = [('lysosomal_degradation', 'Lysosomal\ndegradation'), ('nucleotide_sugar', 'Nucleotide\nsugar synthesis'),
                ('sialylation', 'Terminal\nsialylation')]
    rng = np.random.default_rng(0)
    fig, ax = plt.subplots(figsize=(7, 6))
    for k, (prog, _) in enumerate(programs):
        d = side[side['side'] == prog]
        x = k + rng.uniform(-0.18, 0.18, len(d))
        sig = (d['q_value'] < 0.05).values
        hd = (d['moa'] == 'HDAC inhibitor').values
        ax.scatter(x[~sig & ~hd], d['r_rb'][~sig & ~hd], s=22, c='lightgray')
        ax.scatter(x[sig & ~hd], d['r_rb'][sig & ~hd], s=40, c='#2166ac')
        ax.scatter(np.full(hd.sum(), k), d['r_rb'][hd], s=150, c='#b2182b', edgecolors='black', zorder=3)
    ax.axhline(0, color='gray', lw=1)
    ax.set_xticks(range(len(programs)))
    ax.set_xticklabels([lab for _, lab in programs])
    ax.set_ylabel('Shift vs other compounds (rank-biserial r)')
    ax.set_xlim(-0.5, len(programs) - 0.5)
    fig.tight_layout()
    savefig(fig, 'fig3b_moa_program_shift')
    plt.close(fig)


def plot_gene_level(gl):
    """Figure 3D: lysosomal glycan-degrading enzymes, hydroxamate HDAC inhibitors vs other compounds"""
    g = gl.sort_values('r_rb')
    fig, ax = plt.subplots(figsize=(6, 8))
    y = np.arange(len(g))
    ax.barh(y, g['r_rb'], color=['#b2182b' if q < 0.05 else 'lightgray' for q in g['q_value']])
    ax.set_yticks(y)
    ax.set_yticklabels(g['gene'])
    ax.axvline(0, color='gray', lw=1)
    ax.set_xlabel('Shift vs other compounds (r)')
    fig.tight_layout()
    savefig(fig, 'fig3d_hdaci_gene_level')
    plt.close(fig)


def main():
    ph = annotate()
    classes = explode_moa(ph)
    cz, samples = compound_values()
    mm = moa_modules(ph, classes)
    side, scores = moa_side_shift(ph, cz, classes)
    gl = hdaci(ph, cz, samples, scores)
    style()
    plot_atc()
    plot_moa_modules(mm)
    plot_programs(side)
    plot_gene_level(gl)


if __name__ == '__main__':
    main()
