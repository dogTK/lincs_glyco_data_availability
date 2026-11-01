"""Do microtubule-targeting agents lower the lysosomal glycan-degradation program?

In script 20 the Drug Repurposing Hub class "tubulin polymerization inhibitor" was the only class other than HDAC
inhibitors whose lysosomal degradation score shifted at q < 0.05 (in the opposite direction). This script checks the
same points as for HDAC inhibitors:
1. HepG2: score of each compound, gene-level comparison with all other compounds (HDAC inhibitors excluded),
   and a perturbation-strength-matched permutation (sum |ME| deciles, as in script 20)
2. The same comparison in the nine other cell lines of script 23 (Hub annotation by InChIKey14 / name)
Outputs (results/26_microtubule_agents/): compounds, genes, and cell lines. Not used in the paper (kept for the thesis).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts/20_moa_glyco_programs'))
sys.path.insert(0, str(ROOT / 'scripts/23_hdaci_cross_cell_line'))
from run_moa_programs import (SIDE_GENE_SETS, LYSOSOMAL_PANEL, HUB_DRUGS, HUB_SAMPLES, WIDE_PQ,  # noqa: E402
                                  N_PERM, RNG)
from run_cross_cell_line import fetch, annotate as annotate_hdaci  # noqa: E402

OUT = ROOT / 'results/26_microtubule_agents'
OUT.mkdir(parents=True, exist_ok=True)
MOA = 'tubulin polymerization inhibitor'
LYS = SIDE_GENE_SETS['lysosomal_degradation']
META = ['VALUE', 'canonical_smiles', 'cell', 'cmapid', 'compound_alias', 'dose', 'inchi_key', 'pertid', 'pertname',
        'timepoint', 'sample_id']


def rbc(a, b):
    u, p = mannwhitneyu(a, b)
    return 2 * u / (len(a) * len(b)) - 1, p


def hepg2():
    ph = pd.read_csv(ROOT / 'results/20_moa_glyco_programs/compounds_with_hub.csv')
    moa = ph['moa'].fillna('').str.split('|')
    is_mt = moa.apply(lambda m: MOA in m).values
    is_h = moa.apply(lambda m: 'HDAC inhibitor' in m).values

    wide = pd.read_parquet(WIDE_PQ)
    genes = [c for c in wide.columns if c not in META]
    cz = wide[genes].apply(pd.to_numeric, errors='coerce').groupby(wide['pertname']).median().reindex(ph['compound_id'])
    n_samples = wide.groupby('pertname').size().reindex(ph['compound_id']).values
    lys = cz[LYS].mean(axis=1).values

    comp = pd.DataFrame({'compound_id': ph['compound_id'], 'pert_iname': ph['pert_iname'], 'n_samples': n_samples,
                         'lysosomal_degradation_score': lys,
                         'percentile_among_all': pd.Series(lys).rank(pct=True).values})[is_mt]
    comp.sort_values('lysosomal_degradation_score').to_csv(OUT / 'microtubule_compounds.csv', index=False)
    print(f'HepG2 {MOA} (n={is_mt.sum()}): {(comp["percentile_among_all"] < 0.5).sum()} below the median')

    rows = []
    rest = ~is_mt & ~is_h
    for g in LYSOSOMAL_PANEL:
        r, p = rbc(cz.loc[is_mt, g].values, cz.loc[rest, g].values)
        rows.append(dict(gene=g, microtubule_median_z=cz.loc[is_mt, g].median(), others_median_z=cz.loc[rest, g].median(),
                         r_rb=r, p_value=p))
    gl = pd.DataFrame(rows)
    gl['q_value'] = multipletests(gl['p_value'], method='fdr_bh')[1]
    gl.to_csv(OUT / 'microtubule_gene_level.csv', index=False)
    print(gl.round(4).to_string(index=False))

    me_cols = [c for c in ph.columns if c.startswith('ME')]
    sbin = pd.qcut(ph[me_cols].abs().sum(axis=1), 10, labels=False).values
    pool = ~is_mt & ~is_h
    draws = np.column_stack([RNG.choice(lys[(sbin == b) & pool], N_PERM) for b in sbin[is_mt]])
    null = np.median(draws, axis=1)
    obs = np.median(lys[is_mt])
    p = ((null <= obs).sum() + 1) / (N_PERM + 1)
    print(f'strength-matched permutation: observed median = {obs:.5f}, null median = {np.median(null):.5f}, p = {p:.4f}')


def cross_cell():
    df = annotate_hdaci(fetch())
    dr = pd.read_csv(HUB_DRUGS, sep='\t', comment='!')
    sa = pd.read_csv(HUB_SAMPLES, sep='\t', comment='!')
    names = set(dr.loc[dr['moa'].fillna('').str.split('|').apply(lambda m: MOA in m), 'pert_iname'])
    iks = set(sa.loc[sa['pert_iname'].isin(names), 'InChIKey'].dropna().str[:14])
    df['microtubule'] = df['inchi_key'].fillna('').str[:14].isin(iks) | df['pertname'].str.lower().isin(
        {n.lower() for n in names})
    df['lys'] = df[LYS].mean(axis=1)
    rows = []
    for cell, d in df.groupby('cell'):
        mt, rest = d[d['microtubule']], d[~d['microtubule'] & ~d['hdaci']]
        r, p = rbc(mt['lys'].dropna(), rest['lys'].dropna()) if len(mt) >= 5 else (np.nan, np.nan)
        rows.append(dict(cell=cell, n_microtubule=len(mt), n_other=len(rest), lysosomal_degradation_r=r,
                         lysosomal_degradation_p=p, compounds=';'.join(sorted(mt['pertname']))))
    res = pd.DataFrame(rows).sort_values('lysosomal_degradation_r')
    ok = res['lysosomal_degradation_p'].notna()
    res.loc[ok, 'lysosomal_degradation_q'] = multipletests(res.loc[ok, 'lysosomal_degradation_p'], method='fdr_bh')[1]
    res.to_csv(OUT / 'microtubule_cross_cell_line.csv', index=False)
    print(res.drop(columns='compounds').round(4).to_string(index=False))


if __name__ == '__main__':
    hepg2()
    cross_cell()
