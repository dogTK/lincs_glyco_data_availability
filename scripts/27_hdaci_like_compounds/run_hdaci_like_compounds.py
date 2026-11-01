"""Which compounds outside the HDAC inhibitor class give an HDAC inhibitor-like glycogene response?

1. HDAC inhibitor reference profile: median over the 11 hydroxamate HDAC inhibitor compound IDs of the
   compound-level median z-scores of all 385 glycogenes (HepG2)
2. Spearman correlation of every other compound with the reference profile; compounds annotated as HDAC inhibitors
   in the Drug Repurposing Hub are excluded, and the 95th and 99th percentiles of the remaining compounds are reported
3. For the top compounds: lysosomal degradation score percentile in HepG2 and in the nine other cell lines of
   script 23, and the same values for all SSRIs
Outputs: results/27_hdaci_like_compounds (top compounds, cell lines), Figure 3E and 3F.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts/20_moa_glyco_programs'))
sys.path.insert(0, str(ROOT / 'scripts/23_hdaci_cross_cell_line'))
from run_moa_programs import SIDE_GENE_SETS, HX_IDS, WIDE_PQ, META, savefig  # noqa: E402
from run_cross_cell_line import fetch  # noqa: E402

OUT = ROOT / 'results/27_hdaci_like_compounds'
OUT.mkdir(parents=True, exist_ok=True)
LYS = SIDE_GENE_SETS['lysosomal_degradation']
TOP = 30
SSRI = 'selective serotonin reuptake inhibitor (SSRI)'
# known HDAC inhibitors without the HDAC inhibitor annotation in the Repurposing Hub (among the top 7)
KNOWN_HDACI = {'apicidin': 'apicidin', 'merck-ketone': 'Merck ketone', 'ISOX': 'ISOX', 'HC-toxin': 'HC-toxin'}


def main():
    ph = pd.read_csv(ROOT / 'results/20_moa_glyco_programs/compounds_with_hub.csv').set_index('compound_id')
    wide = pd.read_parquet(WIDE_PQ)
    genes = [c for c in wide.columns if c not in META]
    cz = wide[genes].apply(pd.to_numeric, errors='coerce').groupby(wide['pertname']).median()
    n_samples = wide.groupby('pertname').size()

    ref = cz.loc[cz.index.isin(HX_IDS)].median()
    sim = cz.apply(lambda r: spearmanr(r.values, ref.values)[0], axis=1)
    lys_pct = cz[LYS].mean(axis=1).rank(pct=True)
    hub_hdaci = ph['moa'].fillna('').str.contains('HDAC inhibitor').reindex(cz.index, fill_value=False)
    d = pd.DataFrame({'hdaci_similarity': sim, 'lysosomal_percentile': lys_pct, 'n_samples': n_samples}).join(
        ph[['pert_iname', 'moa', 'clinical_phase']])
    d = d[~(hub_hdaci | d.index.isin(HX_IDS))].sort_values('hdaci_similarity', ascending=False)
    q95, q99 = d['hdaci_similarity'].quantile([0.95, 0.99])
    d['rank'] = range(1, len(d) + 1)
    print(f'{len(d)} compounds outside the HDAC inhibitor class: 95th percentile = {q95:.3f}, 99th = {q99:.3f}')
    top = d.head(TOP).rename_axis('compound_id').reset_index()
    top.to_csv(OUT / 'hdaci_like_top_compounds.csv', index=False)
    print(top.round(3).to_string(index=False))
    launched = d[d['clinical_phase'] == 'Launched']
    print('\nhighest launched drugs:')
    print(launched.head(5).round(3).to_string())

    best = launched.index[0]
    ssris = d[d['moa'].fillna('').str.contains(SSRI, regex=False)].index
    cc = fetch()
    cc['lysosomal_percentile'] = cc[LYS].mean(axis=1).groupby(cc['cell']).rank(pct=True)
    show = cc[cc['pertname'].isin({best, *ssris})][['cell', 'pertname', 'n_samples', 'lysosomal_percentile']]
    show = show.merge(d[['hdaci_similarity']], left_on='pertname', right_index=True, how='left')
    show.sort_values(['pertname', 'cell']).to_csv(OUT / 'paroxetine_ssri_cross_cell_line.csv', index=False)
    print(show.pivot(index='cell', columns='pertname', values='lysosomal_percentile').round(2).to_string())
    print(d.loc[ssris, ['pert_iname', 'hdaci_similarity', 'lysosomal_percentile', 'n_samples']].round(3).to_string())
    plot(d, q99, best, ssris, show)


def plot(d, q99, best, ssris, show):
    plt.rcParams.update({'axes.labelsize': 16, 'xtick.labelsize': 14, 'ytick.labelsize': 14, 'legend.fontsize': 12,
                         'font.size': 14, 'svg.fonttype': 'none'})
    # Fig 4E: similarity to the HDAC inhibitor reference profile, all compounds outside the HDAC inhibitor class
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.hist(d['hdaci_similarity'], bins=60, color='lightgray', edgecolor='gray', linewidth=0.3)
    ax.axvline(q99, color='black', ls='--', lw=1)
    ax.text(q99, ax.get_ylim()[1] * 0.97, ' top 1%', va='top', fontsize=12)
    marks = [(k, lab, '#b2182b') for k, lab in KNOWN_HDACI.items()] + [(best, best, '#e08214')]
    marks += [(k, k, '#2166ac') for k in ssris if k != best]  # LINCS names (escitalopram matches citalopram in the Hub)
    ymax = ax.get_ylim()[1]
    for i, (k, lab, col) in enumerate(sorted(marks, key=lambda m: -d.loc[m[0], 'hdaci_similarity'])):
        x = d.loc[k, 'hdaci_similarity']
        y = ymax * (0.85 - 0.075 * i)
        ax.plot([x, x], [0, y], color=col, lw=1)
        ax.text(x, y, lab, color=col, fontsize=11, ha='right' if x > 0.2 else 'left', va='bottom')
    ax.set_xlabel('Similarity to HDAC inhibitor profile (Spearman ρ)')
    ax.set_ylabel('Compounds')
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    savefig(fig, 'fig3e_hdaci_like_compounds')
    plt.close(fig)

    # Fig 4F: lysosomal degradation score percentile of paroxetine and the other SSRIs in ten cell lines
    order = [best] + sorted([k for k in ssris if k != best], key=lambda k: -d.loc[k, 'hdaci_similarity'])
    fig, ax = plt.subplots(figsize=(7, 6))
    rng = np.random.default_rng(0)
    for i, k in enumerate(order):
        v = show.loc[show['pertname'] == k, 'lysosomal_percentile'].values * 100
        col = '#e08214' if k == best else '#2166ac'
        ax.scatter(i + rng.uniform(-0.15, 0.15, len(v)), v, color=col, s=45, edgecolors='black', linewidths=0.4)
        ax.hlines(np.median(v), i - 0.3, i + 0.3, color='black', lw=2)
    ax.axhline(50, color='gray', ls='--', lw=1)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=40, ha='right')
    ax.set_ylim(0, 100)
    ax.set_ylabel('Lysosomal degradation score\n(percentile among compounds)')
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    savefig(fig, 'fig3f_paroxetine_cross_cell_line')
    plt.close(fig)


if __name__ == '__main__':
    main()
