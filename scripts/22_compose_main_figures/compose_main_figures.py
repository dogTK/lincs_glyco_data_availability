"""Assemble main Figures 1-3 from panel images and write them to manuscript/figure/.

Figure 1: module pathway enrichment dot plot (results/02_wgcna/dtc_module_enrichment_dotplot.png, notebook 02).
Figure 2: fig2_atc_module_gsea (script 20).
Figure 3: A = fig3a_moa_module_shift, B = fig3b_moa_program_shift, D = fig3d_hdaci_gene_level (script 20),
          C = fig3c_hdaci_cross_cell_line (script 23), E and F = fig3e_hdaci_like_compounds and
          fig3f_paroxetine_cross_cell_line (script 27).
"""
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / 'results/figures'
OUT = ROOT / 'manuscript/figure'
DPI = 600
WIDTH_PX = round(178 / 25.4 * DPI)  # 178 mm double-column width
FONT = ImageFont.truetype('/System/Library/Fonts/Supplemental/Arial Bold.ttf', 240)
PAD = 160


def load(name):
    return Image.open(FIG / name / f'{name}.png').convert('RGB')


def trim(img, thresh=245):
    """余白（白）を切り詰める"""
    gray = img.convert('L').point(lambda v: 255 if v < thresh else 0)
    return img.crop(gray.getbbox())


def scale_to(img, width=None, height=None):
    if width:
        return img.resize((width, round(img.height * width / img.width)), Image.LANCZOS)
    return img.resize((round(img.width * height / img.height), height), Image.LANCZOS)


def labeled(img, letter):
    canvas = Image.new('RGB', (img.width + PAD + 300, img.height + PAD + 300), 'white')
    canvas.paste(img, (300, 300))
    ImageDraw.Draw(canvas).text((10, 0), letter, fill='black', font=FONT)
    return canvas


def hstack(imgs):
    h = max(i.height for i in imgs)
    canvas = Image.new('RGB', (sum(i.width for i in imgs), h), 'white')
    x = 0
    for i in imgs:
        canvas.paste(i, (x, 0))
        x += i.width
    return canvas


def vstack(imgs):
    w = max(i.width for i in imgs)
    canvas = Image.new('RGB', (w, sum(i.height for i in imgs)), 'white')
    y = 0
    for i in imgs:
        canvas.paste(i, (0, y))
        y += i.height
    return canvas


def save(img, name):
    """PNG for the docx draft and LZW-compressed TIFF for submission (Glycobiology: line art >= 600 dpi),
    scaled to the double-column width"""
    if img.width > WIDTH_PX:
        img = scale_to(img, width=WIDTH_PX)
    img.save(OUT / f'{name}.png', dpi=(DPI, DPI))
    img.save(OUT / f'{name}.tif', dpi=(DPI, DPI), compression='tiff_lzw')
    print(name, img.size, f'{img.width / DPI * 25.4:.0f} x {img.height / DPI * 25.4:.0f} mm at {DPI} dpi')


def figure1():
    save(trim(Image.open(ROOT / 'results/02_wgcna/dtc_module_enrichment_dotplot.png').convert('RGB')), 'Figure.1')


def figure2():
    save(trim(load('fig2_atc_module_gsea')), 'Figure.2')


def figure3():
    a = labeled(trim(load('fig3a_moa_module_shift')), 'A')
    b = labeled(trim(load('fig3b_moa_program_shift')), 'B')
    c = labeled(trim(load('fig3c_hdaci_cross_cell_line')), 'C')
    d = labeled(trim(load('fig3d_hdaci_gene_level')), 'D')
    e = labeled(trim(load('fig3e_hdaci_like_compounds')), 'E')
    f = labeled(trim(load('fig3f_paroxetine_cross_cell_line')), 'F')
    top = hstack([a, scale_to(b, height=a.height)])
    rows = [top] + [scale_to(hstack([x, scale_to(y, height=x.height)]), width=top.width) for x, y in ((c, d), (e, f))]
    save(vstack(rows), 'Figure.3')


if __name__ == '__main__':
    figure1()
    figure2()
    figure3()
