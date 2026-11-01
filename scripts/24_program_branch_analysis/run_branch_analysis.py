"""Are the nucleotide sugar synthesis and sialylation programs unmoved because opposing branches cancel out?

Splits the nucleotide sugar synthesis and terminal sialylation programs of script 20 into branches (GlycoEnzOnto
nucleotide sugar synthetic pathways with >= 3 program genes; sialyltransferase families) and single genes, and repeats the per-MOA-class comparison
(Mann-Whitney vs all other compounds, rank-biserial r, BH within each branch / gene) in HepG2.
Also compares hydroxamate HDAC inhibitors with other compounds for every sialyltransferase detected in HepG2.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts/20_moa_glyco_programs'))
from run_moa_programs import SIDE_GENE_SETS, MIN_SET  # noqa: E402

OUT = ROOT / 'results/24_program_branch_analysis'
OUT.mkdir(parents=True, exist_ok=True)
BRANCHES = {
    'UDP-GlcNAc': ['AMDHD2', 'GALE', 'GFPT1', 'GFPT2', 'GNPDA1', 'GNPDA2', 'GNPNAT1', 'NAGK', 'PGM3', 'UAP1'],
    'UDP-GalNAc': ['GALE', 'GALK2', 'UAP1'],
    'UDP-Glc': ['GALE', 'GCK', 'GPI', 'HK1', 'HK2', 'HK3', 'PGM1', 'PGM2', 'UGP2'],
    'UDP-Gal': ['GALE', 'GALK1', 'GALT'],
    'GDP-Man': ['GCK', 'GMPPA', 'GMPPB', 'HK1', 'HK2', 'HK3', 'MPI', 'PMM1', 'PMM2'],
    'CMP-Neu5Ac': ['CMAS', 'GNE', 'RENBP'],
    'ST3GAL family': [g for g in SIDE_GENE_SETS['sialylation'] if g.startswith('ST3GAL')],
    'ST6GAL family': [g for g in SIDE_GENE_SETS['sialylation'] if g.startswith('ST6GAL') and 'NAC' not in g],
    'ST6GALNAC family': [g for g in SIDE_GENE_SETS['sialylation'] if g.startswith('ST6GALNAC')],
    'ST8SIA family': [g for g in SIDE_GENE_SETS['sialylation'] if g.startswith('ST8SIA')],
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
    ph = pd.read_csv(ROOT / 'results/20_moa_glyco_programs/compounds_with_hub.csv', low_memory=False)
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
    all_sialyltransferases(wide, ph)



def all_sialyltransferases(wide, ph):
    """All sialyltransferases detected in HepG2: hydroxamate HDAC inhibitors vs compounds other than HDAC inhibitors"""
    from run_moa_programs import HYDROXAMATE
    st = sorted(c for c in wide.columns if c.startswith(('ST3GAL', 'ST6GAL', 'ST8SIA')))
    cz = wide[st].apply(pd.to_numeric, errors='coerce').groupby(wide['pertname']).median()
    hx = set(HYDROXAMATE) | {'BRD-K68202742', 'BRD-A19037878'}  # trichostatin A has two Broad IDs
    hdac = set(ph.loc[ph['moa'].fillna('').str.contains('HDAC inhibitor'), 'compound_id'])
    landmark = set((ROOT / 'inputs/landmark_genes_978.txt').read_text().split())
    rows = []
    for g in st:
        a, b = cz.loc[cz.index.isin(hx), g].dropna(), cz.loc[~cz.index.isin(hdac), g].dropna()
        u, p = mannwhitneyu(a, b)
        rows.append(dict(gene=g, landmark=g in landmark, n_hydroxamate=len(a), n_other=len(b),
                         r_rb=2 * u / (len(a) * len(b)) - 1, p_value=p))
    res = pd.DataFrame(rows)
    res['q_value'] = multipletests(res['p_value'], method='fdr_bh')[1]
    res.sort_values('p_value').to_csv(OUT / 'hdaci_all_sialyltransferases.csv', index=False)
    print('\nAll sialyltransferases, hydroxamate HDAC inhibitors vs other compounds (BH across genes):')
    print(res.sort_values('p_value').round(3).to_string(index=False))

if __name__ == '__main__':
    main()
