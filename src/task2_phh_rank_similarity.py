"""
PHH Rank Similarity Analysis Module
====================================

Using PHH (Primary Human Hepatocytes) as the normal reference, create a glycogene "rank signature" and
score how similar each compound's glycogene rank signature is to PHH (Spearman correlation).

Author: Claude Code
Date: 2026-01-16
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Tuple, Optional, List, Dict, Set
from scipy.stats import spearmanr, pearsonr
from scipy.cluster.hierarchy import linkage, fcluster
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.manifold import TSNE
from umap import UMAP
import logging

# Logging configuration
logger = logging.getLogger(__name__)

# Random seed
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)


def create_rank_signature(
    df: pd.DataFrame,
    gene_columns: List[str],
    sample_id_col: str = 'sample_id',
    ascending: bool = True
) -> pd.DataFrame:
    """
    Convert z-scores to ranks to create rank signature.

    Parameters
    ----------
    df : pd.DataFrame
        z-score data (long or wide format)
    gene_columns : List[str]
        List of glycogene column names
    sample_id_col : str
        Sample ID column name
    ascending : bool
        True: lower z-score = rank 1 (default)
        False: higher z-score = rank 1

    Returns
    -------
    pd.DataFrame
        Long-format DataFrame containing sample_id, gene, rank

    Notes
    -----
    Rank definition: When ascending=True, the gene with the lowest z-score gets rank=1,
    and the gene with the highest z-score gets rank=N (ascending rank).
    Ties are handled using average rank.
    """
    available_genes = [g for g in gene_columns if g in df.columns]

    if len(available_genes) == 0:
        raise ValueError("No glycogene columns found in dataframe")

    results = []

    for _, row in df.iterrows():
        sample_id = row[sample_id_col]
        gene_values = row[available_genes].astype(float)

        # Rank conversion (method='average' for ties)
        ranks = gene_values.rank(ascending=ascending, method='average')

        for gene, rank in ranks.items():
            results.append({
                'sample_id': sample_id,
                'gene': gene,
                'rank': rank,
                'zscore': gene_values[gene]
            })

    return pd.DataFrame(results)


def create_phh_reference_signature(
    df_phh: pd.DataFrame,
    gene_columns: List[str],
    sample_id_col: str = 'sample_id',
    ascending: bool = True
) -> pd.DataFrame:
    """
    Create representative PHH rank signature (median of ranks method).

    Procedure:
    1. For each sample, rank glycogene z-scores across genes
    2. Calculate median of rank values per gene to create representative signature

    Parameters
    ----------
    df_phh : pd.DataFrame
        PHH z-score data (wide format)
    gene_columns : List[str]
        List of glycogene column names
    sample_id_col : str
        Sample ID column name
    ascending : bool
        Rank definition (True: lower z-score = rank 1)

    Returns
    -------
    pd.DataFrame
        DataFrame containing gene, rank_phh, median_zscore, n_samples_used
        n_samples_used is the number of samples with valid rank values
    """
    available_genes = [g for g in gene_columns if g in df_phh.columns]

    if len(available_genes) == 0:
        raise ValueError("No glycogene columns found in dataframe")

    n_samples = len(df_phh)

    # Rank each sample and collect rank values per gene
    gene_ranks = {gene: [] for gene in available_genes}
    gene_zscores = {gene: [] for gene in available_genes}

    for _, row in df_phh.iterrows():
        gene_values = row[available_genes].astype(float)
        # Rank within this sample
        ranks = gene_values.rank(ascending=ascending, method='average')

        for gene in available_genes:
            if pd.notna(ranks[gene]):
                gene_ranks[gene].append(ranks[gene])
            if pd.notna(gene_values[gene]):
                gene_zscores[gene].append(gene_values[gene])

    # Calculate median rank per gene
    results = []
    for gene in available_genes:
        ranks_list = gene_ranks[gene]
        zscores_list = gene_zscores[gene]

        median_rank = np.median(ranks_list) if ranks_list else np.nan
        median_zscore = np.median(zscores_list) if zscores_list else np.nan

        results.append({
            'gene': gene,
            'rank_phh': median_rank,
            'median_zscore': median_zscore,
            'n_samples_used': len(ranks_list)
        })

    df_phh_sig = pd.DataFrame(results)

    return df_phh_sig[['gene', 'rank_phh', 'median_zscore', 'n_samples_used']]


def create_compound_signatures(
    df_lincs: pd.DataFrame,
    gene_columns: List[str],
    compound_id_col: str = 'pertname',
    inchi_key_col: str = 'inchi_key',
    sample_id_col: str = 'sample_id',
    ascending: bool = True
) -> pd.DataFrame:
    """
    Create representative rank signature for each compound (median of ranks method).

    Procedure:
    1. For each condition (sample), rank glycogene z-scores across genes
    2. Calculate median of rank values per compound x gene to create representative signature

    Parameters
    ----------
    df_lincs : pd.DataFrame
        LINCS z-score data (wide format)
    gene_columns : List[str]
        List of glycogene column names
    compound_id_col : str
        Compound ID column name
    inchi_key_col : str
        InChIKey column name
    sample_id_col : str
        Sample ID column name
    ascending : bool
        Rank definition

    Returns
    -------
    pd.DataFrame
        DataFrame containing compound_id, inchi_key, gene, rank_compound, median_zscore, n_conditions_used
        n_conditions_used is the number of conditions with valid rank values
    """
    available_genes = [g for g in gene_columns if g in df_lincs.columns]

    if len(available_genes) == 0:
        raise ValueError("No glycogene columns found in dataframe")

    results = []

    for compound_id, df_comp in df_lincs.groupby(compound_id_col):
        inchi_key = df_comp[inchi_key_col].iloc[0]
        n_conditions = len(df_comp)

        # Rank each condition and collect rank values per gene
        gene_ranks = {gene: [] for gene in available_genes}
        gene_zscores = {gene: [] for gene in available_genes}

        for _, row in df_comp.iterrows():
            gene_values = row[available_genes].astype(float)
            # Rank within this condition
            ranks = gene_values.rank(ascending=ascending, method='average')

            for gene in available_genes:
                if pd.notna(ranks[gene]):
                    gene_ranks[gene].append(ranks[gene])
                if pd.notna(gene_values[gene]):
                    gene_zscores[gene].append(gene_values[gene])

        # Calculate median rank per gene
        for gene in available_genes:
            ranks_list = gene_ranks[gene]
            zscores_list = gene_zscores[gene]

            median_rank = np.median(ranks_list) if ranks_list else np.nan
            median_zscore = np.median(zscores_list) if zscores_list else np.nan

            results.append({
                'compound_id': compound_id,
                'inchi_key': inchi_key,
                'gene': gene,
                'rank_compound': median_rank,
                'median_zscore': median_zscore,
                'n_conditions_used': len(ranks_list)
            })

    return pd.DataFrame(results)


def calculate_phh_similarity(
    df_compound_sig: pd.DataFrame,
    df_phh_sig: pd.DataFrame,
    min_genes: int = 300,
    bootstrap_n: int = 100,
    bootstrap_fraction: float = 0.8,
    random_seed: int = RANDOM_SEED
) -> pd.DataFrame:
    """
    Calculate similarity (Spearman correlation) of each compound to PHH.

    Parameters
    ----------
    df_compound_sig : pd.DataFrame
        Compound rank signature (compound_id, gene, rank_compound)
    df_phh_sig : pd.DataFrame
        PHH rank signature (gene, rank_phh)
    min_genes : int
        Minimum number of genes (compounds with fewer are excluded)
    bootstrap_n : int
        Number of bootstrap iterations
    bootstrap_fraction : float
        Fraction of genes to sample in bootstrap
    random_seed : int
        Random seed

    Returns
    -------
    pd.DataFrame
        DataFrame containing compound_id, inchi_key, similarity_spearman, pearson_corr,
        ci_low, ci_high, n_genes_used, n_conditions_used
    """
    np.random.seed(random_seed)

    # PHH gene set
    phh_genes = set(df_phh_sig['gene'].values)
    phh_rank_dict = df_phh_sig.set_index('gene')['rank_phh'].to_dict()

    results = []
    compounds = df_compound_sig['compound_id'].unique()

    for compound_id in compounds:
        df_comp = df_compound_sig[df_compound_sig['compound_id'] == compound_id]

        # Common genes
        comp_genes = set(df_comp['gene'].values)
        common_genes = list(phh_genes & comp_genes)

        if len(common_genes) < min_genes:
            continue

        # Get InChIKey
        inchi_key = df_comp['inchi_key'].iloc[0]
        n_conditions = df_comp['n_conditions_used'].max()

        # Create rank vectors
        comp_rank_dict = df_comp.set_index('gene')['rank_compound'].to_dict()

        phh_ranks = np.array([phh_rank_dict[g] for g in common_genes])
        comp_ranks = np.array([comp_rank_dict[g] for g in common_genes])

        # Exclude NaN (use only genes valid in both)
        valid_mask = np.isfinite(phh_ranks) & np.isfinite(comp_ranks)
        phh_ranks_valid = phh_ranks[valid_mask]
        comp_ranks_valid = comp_ranks[valid_mask]

        n_valid_genes = len(phh_ranks_valid)
        if n_valid_genes < min_genes:
            continue

        # Spearman correlation
        rho, _ = spearmanr(phh_ranks_valid, comp_ranks_valid)

        # Pearson correlation
        r, _ = pearsonr(phh_ranks_valid, comp_ranks_valid)

        # Bootstrap CI (sampling with replacement)
        bootstrap_rhos = []
        n_sample = max(1, int(n_valid_genes * bootstrap_fraction))  # Ensure at least 1

        for _ in range(bootstrap_n):
            indices = np.random.choice(n_valid_genes, n_sample, replace=True)
            boot_phh = phh_ranks_valid[indices]
            boot_comp = comp_ranks_valid[indices]
            boot_rho, _ = spearmanr(boot_phh, boot_comp)
            if np.isfinite(boot_rho):
                bootstrap_rhos.append(boot_rho)

        # Empty array check
        if len(bootstrap_rhos) > 0:
            ci_low = np.percentile(bootstrap_rhos, 2.5)
            ci_high = np.percentile(bootstrap_rhos, 97.5)
        else:
            ci_low = np.nan
            ci_high = np.nan

        results.append({
            'compound_id': compound_id,
            'inchi_key': inchi_key,
            'similarity_spearman': rho,
            'pearson_corr': r,
            'ci_low': ci_low,
            'ci_high': ci_high,
            'n_genes_used': n_valid_genes,
            'n_conditions_used': n_conditions
        })

    return pd.DataFrame(results)


def identify_contributing_genes(
    df_compound_sig: pd.DataFrame,
    df_phh_sig: pd.DataFrame,
    top_compound_ids: List[str],
    n_genes: int = 30
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Extract genes with smallest/largest rank differences from PHH for top compounds.

    Parameters
    ----------
    df_compound_sig : pd.DataFrame
        Compound rank signature
    df_phh_sig : pd.DataFrame
        PHH rank signature
    top_compound_ids : List[str]
        List of compound IDs to analyze
    n_genes : int
        Number of genes to extract

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        (genes most similar to PHH, genes least similar to PHH)
    """
    phh_rank_dict = df_phh_sig.set_index('gene')['rank_phh'].to_dict()

    closest_results = []
    farthest_results = []

    for compound_id in top_compound_ids:
        df_comp = df_compound_sig[df_compound_sig['compound_id'] == compound_id]

        if len(df_comp) == 0:
            continue

        # Calculate rank difference
        df_comp = df_comp.copy()
        df_comp['phh_rank'] = df_comp['gene'].map(phh_rank_dict)
        df_comp = df_comp.dropna(subset=['phh_rank'])
        df_comp['delta_rank'] = np.abs(df_comp['rank_compound'] - df_comp['phh_rank'])

        # Genes most similar to PHH (small delta_rank)
        closest = df_comp.nsmallest(n_genes, 'delta_rank')[['gene', 'rank_compound', 'phh_rank', 'delta_rank']]
        closest['compound_id'] = compound_id
        closest_results.append(closest)

        # Genes least similar to PHH (large delta_rank)
        farthest = df_comp.nlargest(n_genes, 'delta_rank')[['gene', 'rank_compound', 'phh_rank', 'delta_rank']]
        farthest['compound_id'] = compound_id
        farthest_results.append(farthest)

    df_closest = pd.concat(closest_results, ignore_index=True) if closest_results else pd.DataFrame()
    df_farthest = pd.concat(farthest_results, ignore_index=True) if farthest_results else pd.DataFrame()

    return df_closest, df_farthest


# =============================================================================
# Visualization functions
# =============================================================================

def plot_top_compounds_barplot(
    df_similarity: pd.DataFrame,
    atc_mapping: Optional[pd.DataFrame] = None,
    n_top: int = 30,
    figsize: Tuple[int, int] = (12, 10),
    output_path: Optional[Path] = None
) -> plt.Figure:
    """
    Create similarity bar plot for Top N compounds.

    Parameters
    ----------
    df_similarity : pd.DataFrame
        Similarity data
    atc_mapping : Optional[pd.DataFrame]
        ATC mapping (inchi_key, atc_class)
    n_top : int
        Number of compounds to display
    figsize : Tuple[int, int]
        Figure size
    output_path : Optional[Path]
        Output path

    Returns
    -------
    plt.Figure
        Created figure
    """
    df_top = df_similarity.nlargest(n_top, 'similarity_spearman').copy()

    # Join ATC mapping if available
    if atc_mapping is not None and 'inchi_key' in df_top.columns:
        df_top = df_top.merge(atc_mapping[['inchi_key', 'atc_class']], on='inchi_key', how='left')
        df_top['atc_class'] = df_top['atc_class'].fillna('Unknown')
        has_atc = True
    else:
        has_atc = False

    fig, ax = plt.subplots(figsize=figsize)

    # Bar colors
    if has_atc:
        unique_atc = df_top['atc_class'].unique()
        color_palette = sns.color_palette('tab20', n_colors=len(unique_atc))
        atc_to_color = {atc: color_palette[i] for i, atc in enumerate(unique_atc)}
        colors = [atc_to_color[atc] for atc in df_top['atc_class']]
    else:
        colors = plt.cm.Reds(np.linspace(0.4, 0.8, n_top))[::-1]

    # Bar plot
    y_pos = range(n_top)
    bars = ax.barh(y_pos, df_top['similarity_spearman'].values[::-1], color=colors[::-1], edgecolor='black', linewidth=0.5)

    # Error bars
    ax.errorbar(
        df_top['similarity_spearman'].values[::-1],
        y_pos,
        xerr=[
            (df_top['similarity_spearman'] - df_top['ci_low']).values[::-1],
            (df_top['ci_high'] - df_top['similarity_spearman']).values[::-1]
        ],
        fmt='none',
        color='black',
        capsize=3,
        capthick=1,
        linewidth=1
    )

    # Labels
    labels = df_top['compound_id'].values[::-1]
    if has_atc:
        labels = [f"{c} [{atc}]" if pd.notna(atc) and atc != 'Unknown' else c
                  for c, atc in zip(df_top['compound_id'].values[::-1], df_top['atc_class'].values[::-1])]

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.set_xlabel('Spearman Correlation (Similarity to PHH)', fontweight='bold')
    ax.set_title(f'Top {n_top} Compounds Most Similar to PHH\n(Glycogene Rank Signature)', fontweight='bold')
    ax.axvline(0, color='black', linestyle='-', linewidth=0.5)
    ax.grid(axis='x', alpha=0.3, linestyle='--')

    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
        fig.savefig(output_path.with_suffix('.pdf'), bbox_inches='tight', facecolor='white')
        fig.savefig(output_path.with_suffix('.svg'), format='svg', bbox_inches='tight', facecolor='white')
        logger.info(f"Saved: {output_path}")

    return fig


def plot_similarity_distribution(
    df_similarity: pd.DataFrame,
    figsize: Tuple[int, int] = (10, 6),
    output_path: Optional[Path] = None
) -> plt.Figure:
    """
    Create histogram/density plot of similarity distribution.

    Parameters
    ----------
    df_similarity : pd.DataFrame
        Similarity data
    figsize : Tuple[int, int]
        Figure size
    output_path : Optional[Path]
        Output path

    Returns
    -------
    plt.Figure
        Created figure
    """
    from scipy.stats import gaussian_kde

    fig, ax = plt.subplots(figsize=figsize)

    rho_values = df_similarity['similarity_spearman'].dropna().values

    # Color by positive/negative
    positive = rho_values[rho_values >= 0]
    negative = rho_values[rho_values < 0]

    bins = np.linspace(rho_values.min(), rho_values.max(), 51)
    ax.hist(positive, bins=bins, color='salmon', alpha=0.7, edgecolor='darkred', linewidth=0.5, label=f'ρ ≥ 0 (n={len(positive)})')
    ax.hist(negative, bins=bins, color='steelblue', alpha=0.7, edgecolor='darkblue', linewidth=0.5, label=f'ρ < 0 (n={len(negative)})')

    # Density curve
    if len(rho_values) > 1:
        kde = gaussian_kde(rho_values)
        x_range = np.linspace(rho_values.min(), rho_values.max(), 200)
        ax2 = ax.twinx()
        ax2.plot(x_range, kde(x_range), color='black', linewidth=2, label='Density')
        ax2.set_ylabel('Density', fontweight='bold')
        ax2.set_ylim(0, None)

    # Zero line
    ax.axvline(0, color='red', linestyle='--', linewidth=2, label='ρ = 0')

    # Statistics
    stats_text = f'Mean: {rho_values.mean():.3f}\nMedian: {np.median(rho_values):.3f}\nSD: {rho_values.std():.3f}'
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    ax.set_xlabel(r'Spearman $\rho$ (Similarity to PHH)', fontweight='bold')
    ax.set_ylabel('Number of Compounds', fontweight='bold')
    ax.set_title('Distribution of PHH Similarity Scores', fontweight='bold', pad=15)
    ax.legend(loc='upper right')
    ax.grid(alpha=0.3, linestyle='--')

    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
        fig.savefig(output_path.with_suffix('.pdf'), bbox_inches='tight', facecolor='white')
        fig.savefig(output_path.with_suffix('.svg'), format='svg', bbox_inches='tight', facecolor='white')
        logger.info(f"Saved: {output_path}")

    return fig


def plot_similarity_by_atc(
    df_similarity: pd.DataFrame,
    atc_mapping: pd.DataFrame,
    min_compounds: int = 10,
    figsize: Tuple[int, int] = (14, 8),
    output_path: Optional[Path] = None
) -> plt.Figure:
    """
    Create similarity distribution by ATC class (box/violin plot).

    Parameters
    ----------
    df_similarity : pd.DataFrame
        Similarity data
    atc_mapping : pd.DataFrame
        ATC mapping (inchi_key, atc_class)
    min_compounds : int
        Minimum number of compounds (classes with fewer are grouped into Other)
    figsize : Tuple[int, int]
        Figure size
    output_path : Optional[Path]
        Output path

    Returns
    -------
    plt.Figure
        Created figure
    """
    from scipy.stats import kruskal, mannwhitneyu
    from statsmodels.stats.multitest import multipletests

    # Join ATC mapping
    df_plot = df_similarity.merge(atc_mapping[['inchi_key', 'atc_class']], on='inchi_key', how='left')
    df_plot['atc_class'] = df_plot['atc_class'].fillna('Unknown')

    # Filter by class size
    class_counts = df_plot['atc_class'].value_counts()
    small_classes = class_counts[class_counts < min_compounds].index.tolist()
    df_plot['atc_class_grouped'] = df_plot['atc_class'].apply(lambda x: 'Other' if x in small_classes else x)

    # Sort by median
    order = df_plot.groupby('atc_class_grouped')['similarity_spearman'].median().sort_values(ascending=False).index.tolist()

    fig, ax = plt.subplots(figsize=figsize)

    # Box plot
    sns.boxplot(data=df_plot, x='atc_class_grouped', y='similarity_spearman', order=order, ax=ax, palette='Set2')

    ax.set_xlabel('ATC Class', fontweight='bold')
    ax.set_ylabel(r'Spearman $\rho$ (Similarity to PHH)', fontweight='bold')
    ax.set_title('PHH Similarity by ATC Drug Class', fontweight='bold')
    ax.axhline(0, color='red', linestyle='--', linewidth=1, alpha=0.7)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.grid(axis='y', alpha=0.3, linestyle='--')

    # Kruskal-Wallis test
    groups = [group['similarity_spearman'].values for name, group in df_plot.groupby('atc_class_grouped')]
    if len(groups) >= 2:
        stat, pvalue = kruskal(*groups)
        ax.text(0.98, 0.98, f'Kruskal-Wallis p = {pvalue:.2e}',
                transform=ax.transAxes, ha='right', va='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
        fig.savefig(output_path.with_suffix('.pdf'), bbox_inches='tight', facecolor='white')
        fig.savefig(output_path.with_suffix('.svg'), format='svg', bbox_inches='tight', facecolor='white')
        logger.info(f"Saved: {output_path}")

    return fig


def plot_umap_similarity(
    df_compound_sig: pd.DataFrame,
    df_similarity: pd.DataFrame,
    gene_columns: List[str],
    figsize: Tuple[int, int] = (12, 10),
    output_path: Optional[Path] = None
) -> plt.Figure:
    """
    UMAP visualization of compound rank signatures (highlighting PHH-similar compounds).

    Parameters
    ----------
    df_compound_sig : pd.DataFrame
        Compound rank signature
    df_similarity : pd.DataFrame
        Similarity data
    gene_columns : List[str]
        Glycogene list
    figsize : Tuple[int, int]
        Figure size
    output_path : Optional[Path]
        Output path

    Returns
    -------
    plt.Figure
        Created figure
    """
    # Convert rank signature to wide format
    df_wide = df_compound_sig.pivot(index='compound_id', columns='gene', values='rank_compound')

    # Fill missing values with median
    df_wide = df_wide.fillna(df_wide.median())

    # Join similarity
    compound_ids = df_wide.index.tolist()
    sim_dict = df_similarity.set_index('compound_id')['similarity_spearman'].to_dict()
    similarities = [sim_dict.get(c, np.nan) for c in compound_ids]

    # UMAP
    X = df_wide.values
    reducer = UMAP(n_neighbors=15, min_dist=0.1, n_components=2, random_state=RANDOM_SEED)
    embedding = reducer.fit_transform(X)

    fig, ax = plt.subplots(figsize=figsize)

    # Dynamically set color scale range (excluding outliers)
    valid_sims = [s for s in similarities if np.isfinite(s)]
    if valid_sims:
        vmin = np.percentile(valid_sims, 5)
        vmax = np.percentile(valid_sims, 95)
        # Make symmetric around 0
        abs_max = max(abs(vmin), abs(vmax))
        vmin, vmax = -abs_max, abs_max
    else:
        vmin, vmax = -0.5, 0.5

    # Scatter plot (colored by similarity)
    scatter = ax.scatter(
        embedding[:, 0], embedding[:, 1],
        c=similarities, cmap='RdYlBu', s=30, alpha=0.7,
        vmin=vmin, vmax=vmax, edgecolors='none'
    )

    # Highlight Top 10
    top10_ids = df_similarity.nlargest(10, 'similarity_spearman')['compound_id'].tolist()
    for i, cid in enumerate(compound_ids):
        if cid in top10_ids:
            ax.scatter(embedding[i, 0], embedding[i, 1], c='gold', s=150, edgecolors='black', linewidth=2, marker='*', zorder=10)
            ax.annotate(cid, (embedding[i, 0], embedding[i, 1]), fontsize=7, ha='left', va='bottom')

    plt.colorbar(scatter, ax=ax, label='Similarity to PHH (Spearman ρ)')
    ax.set_xlabel('UMAP 1', fontweight='bold')
    ax.set_ylabel('UMAP 2', fontweight='bold')
    ax.set_title('UMAP of Compound Rank Signatures\n(Stars = Top 10 PHH-similar)', fontweight='bold')

    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
        fig.savefig(output_path.with_suffix('.pdf'), bbox_inches='tight', facecolor='white')
        fig.savefig(output_path.with_suffix('.svg'), format='svg', bbox_inches='tight', facecolor='white')
        logger.info(f"Saved: {output_path}")

    return fig


# =============================================================================
# Report generation
# =============================================================================

def generate_report(
    df_similarity: pd.DataFrame,
    df_phh_sig: pd.DataFrame,
    df_compound_sig: pd.DataFrame,
    atc_mapping: Optional[pd.DataFrame] = None,
    params: dict = None,
    output_path: Path = None
) -> str:
    """
    Generate report Markdown.

    Parameters
    ----------
    df_similarity : pd.DataFrame
        Similarity data
    df_phh_sig : pd.DataFrame
        PHH rank signature
    df_compound_sig : pd.DataFrame
        Compound rank signature
    atc_mapping : Optional[pd.DataFrame]
        ATC mapping
    params : dict
        Parameter dictionary
    output_path : Path
        Output path

    Returns
    -------
    str
        Report Markdown
    """
    if params is None:
        params = {}

    min_genes = params.get('min_genes', 300)
    bootstrap_n = params.get('bootstrap_n', 100)
    bootstrap_fraction = params.get('bootstrap_fraction', 0.8)

    n_compounds = len(df_similarity)
    n_genes = df_phh_sig['n_samples_used'].iloc[0] if len(df_phh_sig) > 0 else 0
    mean_sim = df_similarity['similarity_spearman'].mean()
    median_sim = df_similarity['similarity_spearman'].median()

    # Top 10
    top10 = df_similarity.nlargest(10, 'similarity_spearman')

    if atc_mapping is not None:
        top10 = top10.merge(atc_mapping[['inchi_key', 'atc_class']], on='inchi_key', how='left')
        top10_str = '\n'.join([
            f"  {i+1}. {row['compound_id']}: ρ={row['similarity_spearman']:.3f} (CI: {row['ci_low']:.3f}-{row['ci_high']:.3f}) [{row.get('atc_class', 'N/A')}]"
            for i, row in top10.iterrows()
        ])
    else:
        top10_str = '\n'.join([
            f"  {i+1}. {row['compound_id']}: ρ={row['similarity_spearman']:.3f} (CI: {row['ci_low']:.3f}-{row['ci_high']:.3f})"
            for i, row in top10.iterrows()
        ])

    report = f"""# PHH Rank Similarity Analysis Report

## 1. Objective

Using PHH (Primary Human Hepatocytes) glycogene rank signature as reference,
calculate similarity (Spearman correlation) between each LINCS compound's glycogene rank signature and PHH,
and rank compounds by their similarity to PHH.

## 2. Input Data

- Number of PHH samples: {df_phh_sig['n_samples_used'].max() if len(df_phh_sig) > 0 else 'N/A'}
- Number of LINCS compounds (analyzed): {n_compounds:,}
- Number of glycogenes used: {len(df_phh_sig)}
- Minimum gene filter: {min_genes}

## 3. Methods

### Rank Signature Creation
- For each sample, rank glycogene z-scores across genes
- Rank definition: lower z-score = rank 1 (ascending), ties use average rank
- PHH: when replicates exist, take median(rank) per gene for representative signature
- Compounds: rank per condition, then take median(rank) per gene per compound

### Similarity Calculation
- Spearman correlation coefficient: similarity(d) = SpearmanCorr(rank_compound(d, G), rank_phh(G))
- Pearson correlation coefficient also calculated for reference

### Stability Assessment
- Bootstrap (gene resampling with replacement) for CI estimation
- Sampling rate: {bootstrap_fraction*100:.0f}%
- Number of iterations: {bootstrap_n}
- CI: 2.5%-97.5% percentile

## 4. Results

### Similarity Distribution
- Mean Spearman ρ: {mean_sim:.4f}
- Median Spearman ρ: {median_sim:.4f}
- Maximum: {df_similarity['similarity_spearman'].max():.4f}
- Minimum: {df_similarity['similarity_spearman'].min():.4f}

### Top 10 Compounds (Most Similar to PHH)
{top10_str}

## 5. Output Files

- `phh_rank_signature.parquet`: PHH rank signature
- `compound_rank_signatures.parquet`: Rank signature for each compound
- `phh_similarity_scores.parquet`: Similarity scores
- `top20_genes_closest_to_phh.csv`: Genes most similar to PHH in top compounds
- `top20_genes_farthest_from_phh.csv`: Genes least similar to PHH in top compounds
- `fig_top30_similarity.png`: Top 30 compounds bar chart
- `fig_similarity_distribution.png`: Similarity distribution
- `fig_similarity_by_atc.png`: Similarity by ATC class
- `fig_umap_similarity.png`: UMAP visualization

## 6. Notes

- Random seed: {RANDOM_SEED}
- Rank definition: ascending (lower z-score = smaller rank)
- Analysis date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}
"""

    if output_path:
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(report)
        logger.info(f"Saved: {output_path}")

    return report
