"""Are the nucleotide sugar synthesis and sialylation programs unmoved because opposing branches cancel out?

Splits the two PHH-dissimilar programs of script 20 into branches (GlycoEnzOnto nucleotide-sugar synthetic
pathways with >= 3 genes; sialyltransferase families) and single genes, and repeats the per-MOA-class comparison
(Mann-Whitney vs all other compounds, rank-biserial r, BH within each branch / gene) in HepG2.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts/20_hdac_phh_repurposing_hub'))
from run_hub_phh_analysis import SIDE_GENE_SETS, MIN_SET  # noqa: E402

OUT = ROOT / 'results/24_program_branch_analysis'
OUT.mkdir(parents=True, exist_ok=True)
BRANCHES = {
    'UDP-GlcNAc': ['GALE', 'GFPT1', 'GFPT2', 'GNPNAT1', 'PGM3', 'UAP1'],
    'CMP-Neu5Ac': ['CMAS', 'GNE', 'RENBP'],
    'UDP-Gal': ['GALE', 'GALK1', 'GALT'],
    'UDP-Glc': ['GALE', 'GPI', 'HK2', 'PGM1', 'UGP2'],
    'ST6GALNAC family': ['ST6GALNAC1', 'ST6GALNAC2', 'ST6GALNAC3', 'ST6GALNAC5'],
    'ST8SIA family': ['ST8SIA3', 'ST8SIA4', 'ST8SIA5'],
    'ST3GAL1': ['ST3GAL1'],
    'ST6GAL1': ['ST6GAL1'],
}


def class_tests(scores, d, label_col):
    rows = []
    x = d.assign(category=d['moa'].str.split('|')).explode('category')
    for cat, g in x.groupby('category'):
        if len(g) < MIN_SET:
            continue
        rest = d[~d['compound_id'].isin(g['compound_id'])]
        for name in scores.columns:
            a, b = g[name].dropna(), rest[name].dropna()
            u, p = mannwhitneyu(a, b)
            rows.append({'moa': cat, 'n': len(g), label_col: name, 'r_rb': 2 * u / (len(a) * len(b)) - 1, 'p_value': p})
    res = pd.DataFrame(rows)
    res['q_value'] = res.groupby(label_col)['p_value'].transform(lambda p: multipletests(p, method='fdr_bh')[1])
    return res


def main():
    wide = pd.read_parquet(ROOT / 'data/processed/hepg2_glyco_genes_wide.parquet')
    genes = sorted(set(SIDE_GENE_SETS['nucleotide_sugar']) | set(SIDE_GENE_SETS['sialylation']))
    cz = wide[genes].apply(pd.to_numeric, errors='coerce').groupby(wide['pertname']).median()
    ph = pd.read_csv(ROOT / 'results/20_hdac_phh_repurposing_hub/phh_with_hub.csv', low_memory=False)
    ann = ph[['compound_id', 'moa']].dropna()

    branch_scores = pd.DataFrame({b: cz[g].mean(axis=1) for b, g in BRANCHES.items()})
    d = ann.merge(branch_scores, left_on='compound_id', right_index=True)
    br = class_tests(branch_scores, d, 'branch')
    br.to_csv(OUT / 'moa_branch_shift.csv', index=False)

    d = ann.merge(cz, left_on='compound_id', right_index=True)
    gn = class_tests(cz, d, 'gene')
    gn.to_csv(OUT / 'moa_gene_shift.csv', index=False)

    print(f'{br["moa"].nunique()} MOA classes; branch-level (BH within branch)')
    summ = br.groupby('branch').agg(min_q=('q_value', 'min'), n_q05=('q_value', lambda q: (q < 0.05).sum()),
                                    n_nominal=('p_value', lambda p: (p < 0.05).sum()),
                                    max_abs_r=('r_rb', lambda r: r.abs().max()))
    print(summ.loc[list(BRANCHES)].round(3).to_string())
    print('\nBranch hits q < 0.05:')
    print(br[br['q_value'] < 0.05].sort_values('q_value').round(4).to_string(index=False))
    print('\nHDAC inhibitor by branch:')
    print(br[br['moa'] == 'HDAC inhibitor'][['branch', 'r_rb', 'q_value']].round(3).to_string(index=False))
    print('\nGene-level hits q < 0.05 (BH within gene):')
    print(gn[gn['q_value'] < 0.05].sort_values(['gene', 'q_value']).round(4).to_string(index=False))


if __name__ == '__main__':
    main()
