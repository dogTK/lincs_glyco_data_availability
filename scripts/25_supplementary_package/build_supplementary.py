"""Collect the Supplementary Materials into the files uploaded with the submission.

1. Copies the supplementary figure panels kept in results/ and manuscript/figure/ into
   manuscript/Supplementary_Materials/Supplementary Figures/ (width <= 3,000 px, panel letter on multi-panel figures)
   and writes Supplementary Tables S7-S14 from the outputs of notebook 03 and scripts 20, 23, 24 and 27
   (TABLES), so that SUPPLEMENTARY_MATERIALS{,_EN}.md only point to tracked files. Tables S1-S6 come from
   notebooks 00-02 and are kept as they are.
2. Writes the submission files to manuscript/submission/ (Glycobiology: each file <= 2 MB), split into
   consecutive groups that stay under the limit:
   Supplementary_Figures_Sx-Sy.pdf  legends and figures from SUPPLEMENTARY_MATERIALS_EN.md
   Supplementary_Tables_Sx-Sy.xlsx  a Contents sheet with the table legends and one sheet per table

usage: python scripts/25_supplementary_package/build_supplementary.py
"""
import re
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
from PIL import Image, ImageDraw, ImageFont

Image.MAX_IMAGE_PIXELS = None  # the network figure is 16,632 x 7,319 px

ROOT = Path(__file__).resolve().parents[2]
MS = ROOT / 'manuscript'
SUPP = MS / 'Supplementary_Materials'
FIGDIR = SUPP / 'Supplementary Figures'
TABDIR = SUPP / 'Supplementary Tables'
OUT = MS / 'submission'
ATC = ROOT / 'results/03_atc_pathway_gsea'
MAX_WIDTH = 3000
LIMIT = 2 * 1024 * 1024
FONT = ImageFont.truetype('/System/Library/Fonts/Supplemental/Arial Bold.ttf', 90)

# panel name -> source image (S1-S3 are already in FIGDIR)
PANELS = {
    'S4': ROOT / 'results/02_wgcna/wgcna_dendrogram_modules.png',
    'S5': ROOT / 'results/02_wgcna/module_eigengene_correlation_heatmap.png',
    'S6': ROOT / 'results/02_wgcna/module_specificity_ranking_top10.png',
    'S7a': ATC / 'atc_L1_pathway_biclustered_heatmap.png',
    'S7b': ATC / 'atc_L1_atc_correlation_heatmap.png',
    'S7c': ATC / 'atc_L1_pathway_correlation_heatmap.png',
    'S8a': ATC / 'atc_L2_pathway_biclustered_heatmap.png',
    'S8b': ATC / 'atc_L2_atc_correlation_heatmap.png',
    'S8c': ATC / 'atc_L2_pathway_correlation_heatmap.png',
    'S9a': ATC / 'atc_L3_pathway_biclustered_heatmap.png',
    'S9b': ATC / 'atc_L3_atc_correlation_heatmap.png',
    'S9c': ATC / 'atc_L3_pathway_correlation_heatmap.png',
    'S10': MS / 'figure/Supplementary_Figure_S10.png',
    'S11': ROOT / 'results/02_wgcna/module_intra_gene_network_publication.png',
}


def copy_panels():
    for name, src in PANELS.items():
        img = Image.open(src).convert('RGBA')
        img = Image.alpha_composite(Image.new('RGBA', img.size, 'white'), img).convert('RGB')  # transparent -> white
        if img.width > MAX_WIDTH:
            img = img.resize((MAX_WIDTH, round(img.height * MAX_WIDTH / img.width)), Image.LANCZOS)
        if name[-1] in 'abc':
            canvas = Image.new('RGB', (img.width, img.height + 130), 'white')
            canvas.paste(img, (0, 130))
            ImageDraw.Draw(canvas).text((10, 10), f'({name[-1]})', fill='black', font=FONT)
            img = canvas
        img.save(FIGDIR / f'Supplementary_Figure_{name}.png', optimize=True)
        print('figure', name, img.size)


RES = ROOT / 'results'
# Supplementary Table -> source file in results/ (S7 is written by write_table_s7)
TABLES = {
    'S8': RES / '20_moa_glyco_programs/moa_module_shift.csv',
    'S9a': RES / '20_moa_glyco_programs/hdaci_module_profile.csv',
    'S9b': RES / '20_moa_glyco_programs/hdaci_per_condition_lysosomal.csv',
    'S9c': RES / '20_moa_glyco_programs/hdaci_lysosomal_gene_level.csv',
    'S10': RES / '20_moa_glyco_programs/moa_side_shift.csv',
    'S11': RES / '23_hdaci_cross_cell_line/hdaci_cross_cell_line_program_shift.csv',
    'S12a': RES / '24_program_branch_analysis/moa_branch_shift.csv',
    'S12b': RES / '24_program_branch_analysis/moa_gene_shift.csv',
    'S13': RES / '24_program_branch_analysis/hdaci_all_sialyltransferases.csv',
    'S14a': RES / '27_hdaci_like_compounds/hdaci_like_top_compounds.csv',
    'S14b': RES / '27_hdaci_like_compounds/paroxetine_ssri_cross_cell_line.csv',
}


def write_tables():
    for old in TABDIR.glob('Supplementary_Table_S*.csv'):
        if int(re.match(r'Supplementary_Table_S(\d+)', old.name).group(1)) >= 7:
            old.unlink()
    write_table_s7()
    for name, src in TABLES.items():
        pd.read_csv(src).to_csv(TABDIR / f'Supplementary_Table_{name}.csv', index=False)


def write_table_s7():
    """(a) summary per ATC level, (b) module GSEA for Level 1 (Figure 2), (c) pathway GSEA for Levels 1-3"""
    summary = pd.read_csv(ATC / 'atc_level_comparison_summary.csv')
    summary['Level'] = summary['Level'].str.extract(r'^(L\d)', expand=False)
    summary.to_csv(TABDIR / 'Supplementary_Table_S7a.csv', index=False)
    module = pd.read_csv(ATC / 'atc_L1_module_gsea.csv').rename(columns={'pathway': 'module'})
    module.insert(0, 'atc_class', module.pop('atc_class'))
    module.to_csv(TABDIR / 'Supplementary_Table_S7b.csv', index=False)
    frames = []
    for level in ('L1', 'L2', 'L3'):
        d = pd.read_csv(ATC / f'atc_{level}_pathway_gsea.csv')
        d['atc_level'] = level
        frames.append(d)
    pathway = pd.concat(frames, ignore_index=True)
    lead = ['atc_level', 'atc_class', 'n_compounds', 'n_samples']
    pathway[lead + [c for c in pathway.columns if c not in lead]].to_csv(
        TABDIR / 'Supplementary_Table_S7c.csv', index=False)


def sections(md):
    figs = md.split('## Supplementary Figures', 1)[1].split('## Supplementary Tables', 1)[0]
    tabs = md.split('## Supplementary Tables', 1)[1]
    return figs, tabs


def group_name(prefix, numbers, ext):
    return f'{prefix}_{numbers[0]}-{numbers[-1]}.{ext}' if len(numbers) > 1 else f'{prefix}_{numbers[0]}.{ext}'


def build_pdfs(md, title, tmp):
    """Figures as JPEG (width 1,800 px, quality 85); consecutive figures are added to one PDF while it stays under 2 MB"""
    figs, _ = sections(md)

    def to_jpeg(m):
        src = (MS / m.group(2).replace('%20', ' ')).resolve()
        img = Image.open(src).convert('RGB')
        if img.width > 1800:
            img = img.resize((1800, round(img.height * 1800 / img.width)), Image.LANCZOS)
        dst = tmp / (src.stem + '.jpg')
        img.save(dst, quality=85, optimize=True)
        return f'![]({dst}){{height=65%}}'  # leaves room for the legend on the same page

    body = re.sub(r'!\[([^\]]*)\]\(([^)]+)\)\n?', lambda m: to_jpeg(m) + '\n\n', figs)  # one panel per paragraph.replace('log₁₀', 'log~10~')
    blocks = re.split(r'\n(?=Supplementary Figure S\d+\. )', body.strip())
    figures = [(re.match(r'Supplementary Figure (S\d+)', b).group(1),
                re.sub(r'^(Supplementary Figure S\d+\. )([^\n]+)', r'**\1\2**', b)) for b in blocks]

    def render(group, out):
        first, last = group[0][0], group[-1][0]
        heading = f'Supplementary Figures {first}–{last}' if first != last else f'Supplementary Figure {first}'
        doc = (f'---\ntitle: "{heading}"\nsubtitle: "{title}"\n---\n\n'
               + '\n\n\\newpage\n\n'.join(b for _, b in group))
        (tmp / 'figs.md').write_text(doc)
        subprocess.run(['pandoc', str(tmp / 'figs.md'), '--from', 'markdown-implicit_figures',
                        '--pdf-engine=xelatex', '-V', 'mainfont=Arial', '-V', 'geometry:margin=2cm',
                        '-V', 'fontsize=10pt', '-o', str(out)], check=True, capture_output=True)
        return out.stat().st_size

    outs, group = [], []
    trial = tmp / 'trial.pdf'
    for fig in figures:
        if group and render(group + [fig], trial) > LIMIT:
            outs.append(group)
            group = []
        group.append(fig)
    outs.append(group)
    paths = []
    for g in outs:
        out = OUT / group_name('Supplementary_Figures', [n for n, _ in g], 'pdf')
        render(g, out)
        paths.append(out)
    return paths


def build_xlsxs(md, title, tmp):
    """One sheet per table; consecutive tables are added to one workbook while it stays under 2 MB"""
    _, tabs = sections(md)
    entries = re.findall(r'^(Supplementary Table (S\d+)\..*?)\s*\(((?:`[^`]+`(?:, )?)+)\)\s*$', tabs, flags=re.M)
    tables = []
    for legend, number, files in entries:
        legend = re.sub(r'\*([^*]+)\*', r'\1', legend)
        sheets = [(f.removeprefix('Supplementary_Table_').removesuffix('.csv'), pd.read_csv(TABDIR / f))
                  for f in re.findall(r'`([^`]+)`', files)]
        tables.append((number, legend, sheets))

    def write(group, out):
        contents = [{'Sheet': '', 'Legend': f'Supplementary Tables for: {title}'}]
        contents += [{'Sheet': s, 'Legend': legend} for _, legend, sheets in group for s, _ in sheets]
        with pd.ExcelWriter(out, engine='openpyxl') as xw:
            pd.DataFrame(contents).to_excel(xw, sheet_name='Contents', index=False)
            for _, _, sheets in group:
                for sheet, df in sheets:
                    df.to_excel(xw, sheet_name=sheet, index=False)
        return out.stat().st_size

    outs, group = [], []
    trial = tmp / 'trial.xlsx'
    for table in tables:
        if group and write(group + [table], trial) > LIMIT:
            outs.append(group)
            group = []
        group.append(table)
    outs.append(group)
    paths = []
    for g in outs:
        out = OUT / group_name('Supplementary_Tables', [n for n, _, _ in g], 'xlsx')
        write(g, out)
        paths.append(out)
    return paths


def main():
    copy_panels()
    write_tables()
    OUT.mkdir(exist_ok=True)
    for old in OUT.glob('Supplementary_*'):
        old.unlink()
    md = (MS / 'SUPPLEMENTARY_MATERIALS_EN.md').read_text()
    title = re.search(r'^# (.+)$', (MS / 'FULL_MANUSCRIPT_DRAFT_EN.md').read_text(), flags=re.M).group(1)
    with tempfile.TemporaryDirectory() as tmp:
        for out in build_pdfs(md, title, Path(tmp)) + build_xlsxs(md, title, Path(tmp)):
            size = out.stat().st_size
            print(out.relative_to(ROOT), f'{size / 1024 / 1024:.2f} MB', 'OK' if size <= LIMIT else 'OVER 2 MB')


if __name__ == '__main__':
    main()
