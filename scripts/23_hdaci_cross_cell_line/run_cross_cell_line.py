"""Do HDAC inhibitors raise the lysosomal glycan-degradation program outside HepG2?

For LINCS cell lines with many compounds (and the second liver cancer line HUH7), compare hydroxamate HDAC
inhibitors with all non-HDAC compounds on the three program scores used in script 20 (lysosomal degradation,
nucleotide sugar synthesis, terminal sialylation). Compound-level values are medians of z-scores over all
conditions, computed in Snowflake. Output: Supplementary Table S15.
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from cryptography.hazmat.primitives import serialization
from scipy.stats import mannwhitneyu
from snowflake.connector import connect
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts/20_hdac_phh_repurposing_hub'))
from run_hub_phh_analysis import SIDE_GENE_SETS, HUB_DRUGS, HUB_SAMPLES  # noqa: E402

OUT = ROOT / 'results/23_hdaci_cross_cell_line'
OUT.mkdir(parents=True, exist_ok=True)
CACHE = ROOT / 'data/processed/cross_cell_compound_median_z.parquet'
CELLS = ['HEPG2', 'HUH7', 'VCAP', 'A549', 'A375', 'PC3', 'MCF7', 'HT29', 'HCC515', 'HA1E']
MIN_HDACI = 5


def fetch():
    if CACHE.exists():
        return pd.read_parquet(CACHE)
    # Snowflake key-pair authentication; credentials come from the environment
    key_path = os.path.expanduser(os.environ.get('SNOWFLAKE_PRIVATE_KEY_PATH', '~/.ssh/snowflake_rsa_key.pem'))
    key = serialization.load_pem_private_key(open(key_path, 'rb').read(), password=None)
    der = key.private_bytes(serialization.Encoding.DER, serialization.PrivateFormat.PKCS8,
                            serialization.NoEncryption())
    con = connect(user=os.environ['SNOWFLAKE_USER'], account=os.environ['SNOWFLAKE_ACCOUNT'], private_key=der,
                  warehouse=os.environ.get('SNOWFLAKE_WAREHOUSE'), role=os.environ.get('SNOWFLAKE_ROLE'),
                  database='RAW', schema='LINCS')
    genes = sorted({g for gs in SIDE_GENE_SETS.values() for g in gs})
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


if __name__ == '__main__':
    main()
