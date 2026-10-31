"""PHH similarity axis x Drug Repurposing Hub annotation.

1. Annotate the 4,704 HepG2 compounds with Repurposing Hub MOA/target (BRD -> InChIKey14 -> name)
2. Ranked enrichment along the PHH similarity ranking (same weighted-KS method as notebook 12,
   but N_PERM=10,000 so that the p-value floor does not cap the FDR)
3. HDAC inhibitors: per-condition (dose x time) PHH similarity
4. HDAC inhibitors: gene-level response of PHH-similar / PHH-dissimilar marker genes
5. HDAC inhibitors: module eigengene profile and perturbation-strength-matched control
6. All MOA classes: shift of lysosomal-degradation / nucleotide-sugar / sialylation gene scores
7. Figure 4A, 4B, 4D
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'results/20_hdac_phh_repurposing_hub'
OUT.mkdir(parents=True, exist_ok=True)
SIM_CSV = ROOT / 'results/12_phh_drug_class_enrichment/phh_ranked_compounds.csv'
ME_CSV = ROOT / 'results/02_wgcna/module_eigengenes.csv'
HUB_DRUGS = ROOT / 'data/external/repurposing_drugs_20200324.txt'
HUB_SAMPLES = ROOT / 'data/external/repurposing_samples_20200324.txt'
WIDE_PQ = ROOT / 'data/processed/hepg2_glyco_genes_wide.parquet'
PHH_SIG = ROOT / 'results/04_phh_similarity/phh_rank_signature.parquet'

N_PERM = 10_000
MIN_SET = 5
RNG = np.random.default_rng(42)

HYDROXAMATE = ['PCI-24781', 'PCI-34051', 'vorinostat', 'trichostatin-a', 'scriptaid', 'dacinostat',
               'JNJ-26481585', 'pyroxamide', 'belinostat', 'tubastatin-a']
BENZAMIDE = ['mocetinostat', 'entinostat', 'TC-H-106', 'Merck60']
PHH_SIMILAR_GENES = ['GNS', 'HEXA', 'HEXB', 'NEU1', 'NEU3', 'NEU4', 'GUSB', 'ARSK', 'GLB1',
                     'AGA', 'FUCA1', 'MAN2B1', 'MANBA']
PHH_DISSIMILAR_GENES = ['GFPT1', 'GNPNAT1', 'UAP1', 'UGP2', 'GPI', 'HK2',
                        'ST3GAL1', 'ST6GAL1', 'ST6GALNAC1', 'ST6GALNAC2']
# Genes driving the PHH-similar / PHH-dissimilar pathway enrichments (Supplementary Table S9)
SIDE_GENE_SETS = {
    'lysosomal_degradation': ['ARSK', 'GNS', 'GUSB', 'HEXA', 'HEXB', 'NEU1', 'NEU3', 'NEU4'],
    'nucleotide_sugar': ['CMAS', 'DPM1', 'GALE', 'GALK1', 'GALT', 'GFPT1', 'GFPT2', 'GNE', 'GNPNAT1', 'GPI',
                         'HK2', 'PAPSS2', 'PGM1', 'PGM3', 'PMM2', 'RENBP', 'UAP1', 'UGDH', 'UGP2'],
    'sialylation': ['ST3GAL1', 'ST6GAL1', 'ST6GALNAC1', 'ST6GALNAC2', 'ST6GALNAC3', 'ST6GALNAC5',
                    'ST8SIA3', 'ST8SIA4', 'ST8SIA5'],
}
FIGDIR = ROOT / 'results/figures'


def savefig(fig, name):
    """成果物ごとのフォルダに PNG/PDF/SVG を保存する"""
    outdir = FIGDIR / name
    outdir.mkdir(parents=True, exist_ok=True)
    for ext in ('png', 'pdf', 'svg'):
        fig.savefig(outdir / f'{name}.{ext}', dpi=600, bbox_inches='tight')


# ---------------------------------------------------------------- 1. annotation
def annotate():
    ph = pd.read_csv(SIM_CSV, low_memory=False)[['compound_id', 'inchi_key', 'similarity_spearman', 'rank']]
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
    return pd.concat([ph.reset_index(drop=True), hub], axis=1)


# ---------------------------------------------------------------- 2. ranked enrichment
def ranked_enrichment(ranked, long_frames):
    ranked = ranked.sort_values('similarity_spearman', ascending=False).reset_index(drop=True)
    weights = np.abs(ranked['similarity_spearman'].values)
    ids = ranked['compound_id'].values
    N = len(ids)

    def es_batch(masks):  # masks: (B, N) int8
        hit = masks * weights
        p_hit = hit / hit.sum(1, keepdims=True)
        nh = masks.sum(1, keepdims=True)
        p_miss = (1 - masks) / (N - nh)
        run = np.cumsum(p_hit - p_miss, axis=1)
        idx = np.abs(run).argmax(1)
        return run[np.arange(len(run)), idx], idx

    null_cache = {}

    def null(size):
        if size not in null_cache:
            out = []
            for _ in range(N_PERM // 1000):
                m = np.zeros((1000, N), dtype=np.int8)
                for i in range(1000):
                    m[i, RNG.choice(N, size, replace=False)] = 1
                out.append(es_batch(m)[0])
            null_cache[size] = np.concatenate(out)
        return null_cache[size]

    rows = []
    for ctype, frame in long_frames.items():
        for cat, members in frame.groupby('category')['compound_id']:
            mask = np.isin(ids, members.unique()).astype(np.int8)
            nh = int(mask.sum())
            if nh < MIN_SET:
                continue
            (es,), (peak,) = es_batch(mask[None, :])
            nl = null(nh)
            p = ((nl >= es).sum() + 1) / (N_PERM + 1) if es >= 0 else ((nl <= es).sum() + 1) / (N_PERM + 1)
            le = ids[(np.arange(N) <= peak) & (mask == 1)] if es >= 0 else ids[(np.arange(N) >= peak) & (mask == 1)]
            rows.append(dict(category_type=ctype, category=cat, set_size=nh, ES=es,
                             NES=es / np.abs(nl).mean(), p_value=p,
                             direction='up' if es >= 0 else 'down',
                             leading_edge=';'.join(le)))
    res = pd.DataFrame(rows)
    res['q_value'] = res.groupby('category_type')['p_value'].transform(lambda p: multipletests(p, method='fdr_bh')[1])
    return res.sort_values(['category_type', 'p_value'])


# ---------------------------------------------------------------- 3/4. per-sample data
def per_sample(ph_annot):
    wide = pd.read_parquet(WIDE_PQ)
    phh = pd.read_parquet(PHH_SIG).set_index('gene')['rank_phh']
    genes = [g for g in phh.index if g in wide.columns]
    z = wide[genes].apply(pd.to_numeric, errors='coerce')
    sample_sim = z.apply(lambda r: spearmanr(r.values, phh[genes].values, nan_policy='omit')[0], axis=1)
    samples = wide[['pertname', 'inchi_key', 'dose', 'timepoint']].assign(phh_sim=sample_sim.values)

    # sanity check: compound-level median-rank similarity reproduces the published scores
    med_rank = z.rank(axis=1).groupby(wide['pertname']).median()
    recomputed = med_rank.apply(lambda r: spearmanr(r.values, phh[genes].values)[0], axis=1)
    chk = ph_annot.set_index('compound_id')['similarity_spearman'].reindex(recomputed.index)
    print(f'reproduction check: spearman(recomputed, published) = {spearmanr(recomputed, chk, nan_policy="omit")[0]:.4f}')

    # per-condition percentile of each HDACi sample within all samples at the same time point
    samples['pct_within_time'] = samples.groupby('timepoint')['phh_sim'].rank(pct=True)
    hd = samples[samples['pertname'].isin(HYDROXAMATE + BENZAMIDE + ['BRD-K68202742', 'BRD-A19037878'])].copy()
    hd['chemotype'] = np.where(hd['pertname'].isin(BENZAMIDE), 'benzamide', 'hydroxamate')
    hd.sort_values(['chemotype', 'pertname', 'timepoint', 'dose']).to_csv(OUT / 'hdaci_per_condition_phh_similarity.csv', index=False)
    summ = hd.groupby(['chemotype', 'timepoint']).agg(n=('phh_sim', 'size'), median_sim=('phh_sim', 'median'),
                                                      median_pct=('pct_within_time', 'median'))
    print('\nHDACi per-condition PHH similarity (pct = percentile among all samples at same time point)')
    print(summ.round(3).to_string())

    # gene level: compound median z-score, HDACi (hydroxamate) vs all other compounds
    cz = z.groupby(wide['pertname']).median()
    hx_ids = set(HYDROXAMATE) | {'BRD-K68202742', 'BRD-A19037878'}
    hx = cz.loc[cz.index.isin(hx_ids)]
    rest = cz.loc[~cz.index.isin(hx_ids | set(BENZAMIDE))]
    rows = []
    for side, gl in [('PHH-similar (lysosomal)', PHH_SIMILAR_GENES), ('PHH-dissimilar (NS/sialyl)', PHH_DISSIMILAR_GENES)]:
        for g in gl:
            if g not in cz.columns:
                continue
            u, p = mannwhitneyu(hx[g], rest[g])
            rows.append(dict(side=side, gene=g, hydroxamate_median_z=hx[g].median(), others_median_z=rest[g].median(),
                             benzamide_median_z=cz.loc[cz.index.isin(BENZAMIDE), g].median(),
                             r_rb=2 * u / (len(hx) * len(rest)) - 1, p_value=p))
    gl = pd.DataFrame(rows)
    gl['q_value'] = multipletests(gl['p_value'], method='fdr_bh')[1]
    gl.to_csv(OUT / 'hdaci_marker_gene_response.csv', index=False)
    print('\nMarker genes: hydroxamate HDACi vs other compounds (compound-median z)')
    print(gl.round(3).to_string(index=False))
    return cz


# ---------------------------------------------------------------- 5. HDACi module profile + strength control
def hdaci_modules(ph):
    me_cols = [c for c in ph.columns if c.startswith('ME')]
    is_h = ph['moa'].fillna('').str.contains('HDAC inhibitor')
    h, rest = ph[is_h], ph[~is_h]
    rows = []
    for c in me_cols:
        u, p = mannwhitneyu(h[c], rest[c])
        rows.append(dict(module=c, hdaci_median=h[c].median(), others_median=rest[c].median(),
                         r_rb=2 * u / (len(h) * len(rest)) - 1, p_value=p))
    mp = pd.DataFrame(rows)
    mp['q_value'] = multipletests(mp['p_value'], method='fdr_bh')[1]
    mp.to_csv(OUT / 'hdaci_module_profile.csv', index=False)
    print(f'\nHDACi (n={len(h)}) module eigengene profile vs other compounds')
    print(mp.round(4).to_string(index=False))

    # perturbation strength (sum |ME|) is not what places HDACi on the PHH-similar side
    strength = ph[me_cols].abs().sum(1)
    print(f'spearman(PHH similarity, strength) = {spearmanr(ph["similarity_spearman"], strength)[0]:.3f}')
    sbin = pd.qcut(strength, 10, labels=False)
    sim = ph['similarity_spearman'].values
    draws = np.column_stack([RNG.choice(sim[(sbin == b).values & ~is_h.values], N_PERM) for b in sbin[is_h]])
    null = np.median(draws, axis=1)
    obs = np.median(sim[is_h.values])
    p = ((null >= obs).sum() + 1) / (N_PERM + 1)
    print(f'strength-matched permutation: observed median rho = {obs:.3f}, null median = {np.median(null):.3f}, p = {p:.4f}')


# ---------------------------------------------------------------- 6. which MOA classes move which side
def moa_side_shift(ph, cz):
    scores = pd.DataFrame({side: cz[[g for g in genes if g in cz.columns]].mean(1)
                           for side, genes in SIDE_GENE_SETS.items()})
    d = ph[['compound_id', 'moa']].merge(scores, left_on='compound_id', right_index=True)
    x = d.assign(category=d['moa'].str.split('|')).explode('category')
    rows = []
    for cat, g in x.groupby('category'):
        if len(g) < MIN_SET:
            continue
        rest = d[~d['compound_id'].isin(g['compound_id'])]
        for side in SIDE_GENE_SETS:
            u, p = mannwhitneyu(g[side], rest[side])
            rows.append(dict(moa=cat, n=len(g), side=side, r_rb=2 * u / (len(g) * len(rest)) - 1, p_value=p))
    res = pd.DataFrame(rows)
    res['q_value'] = res.groupby('side')['p_value'].transform(lambda p: multipletests(p, method='fdr_bh')[1])
    res.to_csv(OUT / 'moa_side_shift.csv', index=False)
    print(f'\nMOA classes (n >= {MIN_SET}, {res["moa"].nunique()} classes): shift of side scores vs other compounds, q < 0.05')
    for side in SIDE_GENE_SETS:
        s = res[(res['side'] == side) & (res['q_value'] < 0.05)].sort_values('r_rb')
        print(f'  {side}: up {(s["r_rb"] > 0).sum()}, down {(s["r_rb"] < 0).sum()}')
        if len(s):
            print(s[['moa', 'n', 'r_rb', 'q_value']].round(4).to_string(index=False))
    return res


# ---------------------------------------------------------------- 7. figures
def plot_figures(enr, side):
    plt.rcParams.update({'axes.labelsize': 16, 'axes.titlesize': 16, 'xtick.labelsize': 14,
                         'ytick.labelsize': 14, 'legend.fontsize': 14, 'font.size': 14, 'svg.fonttype': 'none'})

    # Fig 4A: ranked enrichment of Hub MOA classes along the PHH similarity ranking
    e = enr[(enr['category_type'] == 'HUB_MOA') & (enr['p_value'] < 0.05)].sort_values('NES')
    fig, ax = plt.subplots(figsize=(8, 9))
    sc = ax.scatter(e['NES'], np.arange(len(e)), s=e['set_size'] * 10, c=-np.log10(e['q_value']),
                    cmap='Reds', vmin=0, edgecolors='black', linewidths=0.5)
    ax.set_yticks(np.arange(len(e)))
    ax.set_yticklabels([f'{c} (n={n})' for c, n in zip(e['category'], e['set_size'])])
    ax.axvline(0, color='gray', lw=1)
    ax.set_xlabel('NES (positive: PHH-similar)')
    fig.colorbar(sc, ax=ax, label='−log10(q)', shrink=0.5)
    fig.tight_layout()
    savefig(fig, 'fig4a_hub_moa_phh_ranked_enrichment')
    plt.close(fig)

    # Fig 4B: shift of each program score per MOA class (one dot per class); HDAC inhibitors in red
    programs = [('lysosomal_degradation', 'Lysosomal\ndegradation'), ('nucleotide_sugar', 'Nucleotide\nsugar synthesis'),
                ('sialylation', 'Terminal\nsialylation')]
    rng = np.random.default_rng(0)
    fig, ax = plt.subplots(figsize=(6, 9))
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
    savefig(fig, 'fig4b_moa_program_shift')
    plt.close(fig)


def plot_gene_level(genes):
    """Figure 4D: per-gene shift (hydroxamate HDAC inhibitors vs other compounds), PHH-similar vs PHH-dissimilar genes"""
    plt.rcParams.update({'axes.labelsize': 16, 'xtick.labelsize': 14, 'ytick.labelsize': 13, 'font.size': 14,
                         'svg.fonttype': 'none'})
    g = genes.copy()
    g['group'] = np.where(g['side'].str.startswith('PHH-similar'), 'Lysosomal degradation',
                          np.where(g['gene'].str.startswith('ST'), 'Terminal sialylation', 'Nucleotide sugar synthesis'))
    order = ['Lysosomal degradation', 'Nucleotide sugar synthesis', 'Terminal sialylation']
    g['group'] = pd.Categorical(g['group'], order)
    g = g.sort_values(['group', 'r_rb'], ascending=[False, True])
    colors = {'Lysosomal degradation': '#b2182b', 'Nucleotide sugar synthesis': '#2166ac', 'Terminal sialylation': '#5e9fd1'}
    fig, ax = plt.subplots(figsize=(6, 10))
    y = np.arange(len(g))
    ax.barh(y, g['r_rb'], color=[colors[k] if q < 0.05 else 'lightgray' for k, q in zip(g['group'], g['q_value'])])
    ax.set_yticks(y)
    ax.set_yticklabels(g['gene'])
    ax.axvline(0, color='gray', lw=1)
    ax.set_xlabel('Shift vs other compounds (r)')
    for k in order:
        ax.barh([np.nan], [np.nan], color=colors[k], label=k)
    ax.legend(loc='upper center', bbox_to_anchor=(0.4, -0.09), fontsize=11, frameon=False)
    fig.tight_layout()
    savefig(fig, 'fig4d_hdaci_gene_level')
    plt.close(fig)

def main():
    ph = annotate()
    me = pd.read_csv(ME_CSV).rename(columns={'pertname': 'compound_id'})
    ph.merge(me, on='compound_id', how='left').to_csv(OUT / 'phh_with_hub.csv', index=False)
    print(f'Hub MOA coverage: {ph["moa"].notna().mean():.1%} ({ph["moa"].notna().sum()} compounds); match:',
          ph['match'].value_counts().to_dict())

    long_frames = {}
    for col, ctype in [('moa', 'HUB_MOA'), ('target', 'HUB_TARGET')]:
        f = ph[['compound_id', col]].dropna()
        f = f.assign(category=f[col].str.split('|')).explode('category')[['compound_id', 'category']].drop_duplicates()
        long_frames[ctype] = f
    enr = ranked_enrichment(ph, long_frames)
    enr.to_csv(OUT / 'phh_ranked_enrichment_hub.csv', index=False)
    print('\nRanked enrichment (Hub), q < 0.25:')
    print(enr[enr['q_value'] < 0.25].drop(columns='leading_edge').round(4).to_string(index=False))

    cz = per_sample(ph)
    ph_me = ph.merge(me, on='compound_id')
    hdaci_modules(ph_me)
    side = moa_side_shift(ph, cz)
    plot_figures(enr, side)
    plot_gene_level(pd.read_csv(OUT / 'hdaci_marker_gene_response.csv'))


if __name__ == '__main__':
    main()
