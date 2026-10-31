"""PHH reference axis: genes driving the similarity of the top-30 PHH-similar compounds.

Uses the current outputs of notebook 04 (phh_similarity_scores.csv, compound_rank_signatures.parquet,
phh_rank_signature.parquet) and the same procedure: for each of the top-30 compounds take the 20 genes
with the smallest / largest |rank_compound - rank_PHH|, keep genes appearing in >= 3 compounds, and run
GlycoEnzOnto over-representation analysis (Fisher one-sided, background = all GlycoEnzOnto genes).
Outputs Supplementary Tables S7-S9 and Figure 2A/2B.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'notebooks/lib'))  # this repository
sys.path.insert(0, str(ROOT / 'src'))  # lincs_glyco_data_availability
import task2_phh_rank_similarity as phh  # noqa: E402

R04 = ROOT / 'results/04_phh_similarity'
OUT = ROOT / 'results/21_phh_axis_contributing_genes'
OUT.mkdir(parents=True, exist_ok=True)
FIGDIR = ROOT / 'results/figures'
N_TOP, N_GENES, MIN_FREQ = 30, 20, 3


def savefig(fig, name):
    """成果物ごとのフォルダに PNG/PDF/SVG を保存する"""
    outdir = FIGDIR / name
    outdir.mkdir(parents=True, exist_ok=True)
    for ext in ('png', 'pdf', 'svg'):
        fig.savefig(outdir / f'{name}.{ext}', dpi=300, bbox_inches='tight')


def read_gmt(path):
    pw = {}
    for line in open(path):
        parts = line.rstrip('\n').split('\t')
        if len(parts) >= 3:
            pw[parts[0].strip('"')] = set(parts[2:])
    return pw


def ora(genes, pathways, background):
    genes = set(genes) & background
    rows = []
    for name, pg in pathways.items():
        pg = pg & background
        a = len(genes & pg)
        if a == 0:
            continue
        b, c = len(genes - pg), len(pg - genes)
        d = len(background - genes - pg)
        p = fisher_exact([[a, b], [c, d]], alternative='greater')[1]
        rows.append(dict(pathway=name, overlap=a, gene_set_size=len(genes), pathway_size=len(pg),
                         fold_enrichment=(a / len(genes)) / (len(pg) / len(background)), pvalue=p,
                         genes=','.join(sorted(genes & pg))))
    df = pd.DataFrame(rows)
    df['qvalue'] = multipletests(df['pvalue'], method='fdr_bh')[1]
    return df.sort_values('pvalue')


def main():
    sim = pd.read_csv(R04 / 'phh_similarity_scores.csv').sort_values('similarity_spearman', ascending=False)
    sim.to_csv(OUT / 'phh_similarity_scores.csv', index=False)
    sig = pd.read_parquet(R04 / 'compound_rank_signatures.parquet')
    phh_sig = pd.read_parquet(R04 / 'phh_rank_signature.parquet')
    top = sim.head(N_TOP)
    closest, farthest = phh.identify_contributing_genes(sig, phh_sig, top['compound_id'].tolist(), n_genes=N_GENES)
    closest.to_csv(OUT / 'top30_genes_closest_to_phh.csv', index=False)
    farthest.to_csv(OUT / 'top30_genes_farthest_from_phh.csv', index=False)

    pathways = read_gmt(ROOT / 'inputs/GlycoEnzOnto/GlycoEnzOnto.gmt')
    background = set().union(*pathways.values())
    res = []
    for label, df in [('PHH_similar', closest), ('PHH_dissimilar', farthest)]:
        freq = df['gene'].value_counts()
        genes = freq[freq >= MIN_FREQ].index
        pd.Series(sorted(genes), name='gene').to_csv(OUT / f'{label}_genes.csv', index=False)
        e = ora(genes, pathways, background).assign(direction=label)
        res.append(e)
        print(f'{label}: {len(genes)} genes (background {len(background)})')
        print(e.head(6)[['pathway', 'overlap', 'pathway_size', 'fold_enrichment', 'pvalue', 'qvalue', 'genes']]
              .round(4).to_string(index=False))
    enr = pd.concat(res)
    enr.to_csv(OUT / 'phh_contributing_genes_pathway.csv', index=False)

    print('\nTop 5 compounds:')
    print(top.head(5)[['compound_id', 'similarity_spearman', 'ci_low', 'ci_high', 'n_conditions_used']].round(3).to_string(index=False))
    print('similarity range %.3f to %.3f, mean %.4f, median %.4f, |rho|<0.1: %.1f%%' % (
        sim.similarity_spearman.min(), sim.similarity_spearman.max(), sim.similarity_spearman.mean(),
        sim.similarity_spearman.median(), 100 * (sim.similarity_spearman.abs() < 0.1).mean()))
    print('n_conditions_used distribution:', sim['n_conditions_used'].value_counts().sort_index().head(6).to_dict())

    plt.rcParams.update({'axes.labelsize': 18, 'axes.titlesize': 16, 'xtick.labelsize': 16,
                         'ytick.labelsize': 14, 'legend.fontsize': 16, 'font.size': 16, 'svg.fonttype': 'none'})
    # Figure 2A: top-30 compounds with bootstrap 95% CI
    t = top.iloc[::-1]
    fig, ax = plt.subplots(figsize=(10, 12))
    ax.errorbar(t['similarity_spearman'], np.arange(len(t)),
                xerr=[t['similarity_spearman'] - t['ci_low'], t['ci_high'] - t['similarity_spearman']],
                fmt='o', color='#b2182b', ecolor='gray', capsize=3)
    ax.set_yticks(np.arange(len(t)))
    ax.set_yticklabels(t['compound_id'])
    ax.set_xlabel('Similarity to PHH (Spearman ρ)')
    fig.tight_layout()
    savefig(fig, 'fig2a_phh_top30_similarity')
    plt.close(fig)

    # Figure 2B: pathway enrichment of PHH-similar / PHH-dissimilar genes (p < 0.05)
    sub = {lab: enr[(enr['direction'] == lab) & (enr['pvalue'] < 0.05)].sort_values('pvalue', ascending=False)
           for lab in ('PHH_similar', 'PHH_dissimilar')}
    fig, axes = plt.subplots(2, 1, figsize=(12, 14), sharex=True,
                             gridspec_kw={'height_ratios': [len(sub['PHH_similar']), len(sub['PHH_dissimilar'])]})
    for ax, label, color in [(axes[0], 'PHH_similar', '#b2182b'), (axes[1], 'PHH_dissimilar', '#2166ac')]:
        e = sub[label]
        ax.barh(e['pathway'], -np.log10(e['pvalue']), color=color)
        ax.axvline(-np.log10(0.05), color='gray', ls='--')
        ax.set_title(label.replace('_', '-'))
    axes[1].set_xlabel('−log10(p)')
    fig.tight_layout()
    savefig(fig, 'fig2b_phh_contributing_genes_pathway')
    plt.close(fig)

    # Supplementary Figure S12: frequency of contributing genes among the top-30 compounds
    fig, axes = plt.subplots(1, 2, figsize=(14, 16))
    for ax, (label, df), color in [(axes[0], ('PHH_similar', closest), '#b2182b'),
                                   (axes[1], ('PHH_dissimilar', farthest), '#2166ac')]:
        freq = df['gene'].value_counts()
        freq = freq[freq >= MIN_FREQ].head(40).iloc[::-1]
        ax.barh(freq.index, freq.values, color=color)
        ax.set_xlabel(f'Compounds (of top {N_TOP})')
        ax.set_title(label.replace('_', '-'))
        ax.tick_params(axis='y', labelsize=11)
    fig.tight_layout()
    savefig(fig, 'figS12_phh_contributing_genes_frequency')
    plt.close(fig)

if __name__ == '__main__':
    main()
