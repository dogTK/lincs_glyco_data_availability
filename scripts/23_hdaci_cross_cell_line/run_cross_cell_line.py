"""Do HDAC inhibitors raise the lysosomal glycan-degradation program outside HepG2?

For LINCS cell lines with many compounds (and the second liver cancer line HUH7), compare hydroxamate HDAC
inhibitors with all non-HDAC compounds on the three program scores used in script 20 (lysosomal degradation,
nucleotide sugar synthesis, terminal sialylation). Compound-level values are medians of z-scores over all
conditions, computed in Snowflake. Outputs: results/23_hdaci_cross_cell_line and Figure 3C.
"""
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cryptography.hazmat.primitives import serialization
from scipy.stats import mannwhitneyu
from snowflake.connector import connect
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts/20_moa_glyco_programs'))
from run_moa_programs import SIDE_GENE_SETS, LYSOSOMAL_PANEL, HUB_DRUGS, HUB_SAMPLES  # noqa: E402

OUT = ROOT / 'results/23_hdaci_cross_cell_line'
OUT.mkdir(parents=True, exist_ok=True)
CACHE = ROOT / 'data/processed/cross_cell_compound_median_z.parquet'
CELLS = ['HEPG2', 'HUH7', 'VCAP', 'A549', 'A375', 'PC3', 'MCF7', 'HT29', 'HCC515', 'HA1E']
MIN_HDACI = 5


GENES = sorted({g for gs in SIDE_GENE_SETS.values() for g in gs} | set(LYSOSOMAL_PANEL))


def credentials():
    """Snowflake key-pair credentials from the environment, else from the default connection in
    ~/.snowflake/config.toml"""
    if 'SNOWFLAKE_USER' in os.environ:
        return dict(user=os.environ['SNOWFLAKE_USER'], account=os.environ['SNOWFLAKE_ACCOUNT'],
                    private_key_path=os.environ.get('SNOWFLAKE_PRIVATE_KEY_PATH', '~/.ssh/snowflake_rsa_key.pem'),
                    warehouse=os.environ.get('SNOWFLAKE_WAREHOUSE'), role=os.environ.get('SNOWFLAKE_ROLE'))
    import tomllib
    return tomllib.load(open(os.path.expanduser('~/.snowflake/config.toml'), 'rb'))['connections']['default']


def fetch():
    if CACHE.exists():
        df = pd.read_parquet(CACHE)
        if set(GENES) <= set(df.columns):
            return df
    c = credentials()
    key = serialization.load_pem_private_key(open(os.path.expanduser(c['private_key_path']), 'rb').read(), password=None)
    der = key.private_bytes(serialization.Encoding.DER, serialization.PrivateFormat.PKCS8,
                            serialization.NoEncryption())
    con = connect(user=c['user'], account=c['account'], private_key=der, warehouse=c.get('warehouse'),
                  role=c.get('role'), database='RAW', schema='LINCS')
    genes = GENES
    cols = ', '.join(f'MEDIAN(TRY_TO_DOUBLE("{g}")) AS "{g}"' for g in genes)
    cells = ', '.join(f"'{c}'" for c in CELLS)
    q = (f'SELECT "cell", "pertname", MAX("inchi_key") AS "inchi_key", COUNT(*) AS "n_samples", {cols} '
         f'FROM RAW.LINCS.GLYCO_GENES_WIDE WHERE "cell" IN ({cells}) GROUP BY "cell", "pertname"')
    df = pd.read_sql(q, con)
    df.to_parquet(CACHE)
    return df


def annotate(df):
    """Flag HDAC inhibitors (any chemotype) and hydroxamate HDAC inhibitors via Repurposing Hub InChIKey14 / name."""
    dr = pd.read_csv(HUB_DRUGS, sep='\t', comment='!')
    sa = pd.read_csv(HUB_SAMPLES, sep='\t', comment='!')
    hdac_names = set(dr.loc[dr['moa'].fillna('').str.contains('HDAC inhibitor'), 'pert_iname'])
    sa = sa[sa['pert_iname'].isin(hdac_names)]
    hdac_ik = set(sa['InChIKey'].dropna().str[:14])
    hydroxamate = {'PCI-24781', 'PCI-34051', 'vorinostat', 'trichostatin-a', 'scriptaid', 'dacinostat',
                   'JNJ-26481585', 'pyroxamide', 'belinostat', 'tubastatin-a'}
    hx_ik = set(sa.loc[sa['pert_iname'].isin(hydroxamate), 'InChIKey'].dropna().str[:14])
    ik = df['inchi_key'].fillna('').str[:14]
    name = df['pertname'].str.lower()
    df['hdaci'] = ik.isin(hdac_ik) | name.isin({n.lower() for n in hdac_names})
    df['hydroxamate'] = ik.isin(hx_ik) | name.isin({n.lower() for n in hydroxamate})
    return df


def main():
    df = annotate(fetch())
    for side, genes in SIDE_GENE_SETS.items():
        df[side] = df[[g for g in genes if g in df.columns]].mean(axis=1)
    rows = []
    for cell, d in df.groupby('cell'):
        hx = d[d['hydroxamate']]
        rest = d[~d['hdaci']]
        row = dict(cell=cell, n_hydroxamate=len(hx), n_other=len(rest),
                   hydroxamates=';'.join(sorted(hx['pertname'])))
        for side in SIDE_GENE_SETS:
            if len(hx) >= MIN_HDACI:
                u, p = mannwhitneyu(hx[side].dropna(), rest[side].dropna())
                row[f'{side}_r'] = 2 * u / (hx[side].notna().sum() * rest[side].notna().sum()) - 1
                row[f'{side}_p'] = p
            else:
                row[f'{side}_r'] = row[f'{side}_p'] = np.nan
        rows.append(row)
    res = pd.DataFrame(rows).sort_values('n_hydroxamate', ascending=False)
    for side in SIDE_GENE_SETS:
        ok = res[f'{side}_p'].notna()
        res.loc[ok, f'{side}_q'] = multipletests(res.loc[ok, f'{side}_p'], method='fdr_bh')[1]
    res.to_csv(OUT / 'hdaci_cross_cell_line_program_shift.csv', index=False)
    show = ['cell', 'n_hydroxamate', 'n_other'] + [f'{s}_{k}' for s in SIDE_GENE_SETS for k in ('r', 'q')]
    print(res[show].round(3).to_string(index=False))
    plot_cross_cell(res)



def plot_cross_cell(res):
    """Figure 3C: rank-biserial r of hydroxamate HDAC inhibitors per cell line and program (* q < 0.05)"""
    plt.rcParams.update({'axes.labelsize': 16, 'xtick.labelsize': 14, 'ytick.labelsize': 14, 'font.size': 14,
                         'svg.fonttype': 'none'})
    progs = [('lysosomal_degradation', 'Lysosomal\ndegradation'), ('nucleotide_sugar', 'Nucleotide\nsugar synthesis'),
             ('sialylation', 'Terminal\nsialylation')]
    d = res.dropna(subset=['lysosomal_degradation_r']).copy()
    d['order'] = (d['cell'] != 'HEPG2').astype(int)
    d = d.sort_values(['order', 'lysosomal_degradation_r'], ascending=[True, False])
    r = d[[f'{p}_r' for p, _ in progs]].to_numpy()
    q = d[[f'{p}_q' for p, _ in progs]].to_numpy()
    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(r, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
    for i in range(r.shape[0]):
        for j in range(r.shape[1]):
            ax.text(j, i, f'{r[i, j]:.2f}' + ('*' if q[i, j] < 0.05 else ''), ha='center', va='center', fontsize=12,
                    color='white' if abs(r[i, j]) > 0.6 else 'black')
    ax.set_xticks(range(len(progs)))
    ax.set_xticklabels([lab for _, lab in progs], fontsize=12)
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels([c if c != 'HEPG2' else 'HepG2' for c in d['cell']])
    fig.colorbar(im, ax=ax, label='Rank-biserial r', shrink=0.6)
    fig.tight_layout()
    outdir = ROOT / 'results/figures/fig3c_hdaci_cross_cell_line'
    outdir.mkdir(parents=True, exist_ok=True)
    for ext in ('png', 'pdf', 'svg'):
        fig.savefig(outdir / f'fig3c_hdaci_cross_cell_line.{ext}', dpi=600, bbox_inches='tight')
    plt.close(fig)

if __name__ == '__main__':
    main()
