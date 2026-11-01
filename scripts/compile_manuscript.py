"""Compile the manuscript markdown to Word.

draft:      legends stay under each figure (for reading and the thesis draft); alt text is omitted
            -> manuscript/FULL_MANUSCRIPT_DRAFT_{JA,EN}.docx
submission: Glycobiology section order (... References, Legends to figures with alt text under each legend, Tables);
            figure images are removed
            because figures are uploaded as separate files (manuscript/figure/Figure.N.tif)
            -> manuscript/FULL_MANUSCRIPT_SUBMISSION_{JA,EN}.docx

usage: python scripts/compile_manuscript.py [draft|submission|both] [JA|EN|both]
"""
import re
import subprocess
import sys
import tempfile
from pathlib import Path

MS = Path(__file__).resolve().parents[1] / 'manuscript'
PANDOC = ['pandoc', '--from', 'markdown-implicit_figures', '--to', 'docx', '--citeproc',
          '--bibliography=refs/lincs_glyco_2nd_paper.bib', '--csl=cse-name-year-8th-edition.csl']


ALT = re.compile(r'\*Alt text \(Figure \d+\)\.\*')


def to_draft(md):
    """alt text is only for the submission file"""
    return '\n\n'.join(b for b in md.split('\n\n') if not ALT.match(b.strip()))


def to_submission(md):
    blocks = md.split('\n\n')
    legends, tables, kept = [], [], []
    for b in blocks:
        s = b.strip()
        if re.fullmatch(r'!\[Figure \d+\]\([^)]*\)', s):
            continue  # figures are uploaded separately
        if re.match(r'\*\*Figure \d+\.\*\*', s) or ALT.match(s):
            legends.append(s)  # Glycobiology: alt text directly under each legend
        elif s.startswith('|') or s.startswith(': Table'):
            tables.append(s)
        else:
            kept.append(b)
    out = '\n\n'.join(kept).rstrip()
    assert out.endswith('## References'), 'References must be the last section of the source'
    return (out + '\n\n::: {#refs}\n:::\n\n## Legends to figures\n\n' + '\n\n'.join(legends)
            + '\n\n## Tables\n\n' + '\n\n'.join(tables) + '\n')


def compile_one(mode, lang):
    src = MS / f'FULL_MANUSCRIPT_DRAFT_{lang}.md'
    if mode == 'draft':
        md, out = to_draft(src.read_text()), MS / f'FULL_MANUSCRIPT_DRAFT_{lang}.docx'
    else:
        md, out = to_submission(src.read_text()), MS / f'FULL_MANUSCRIPT_SUBMISSION_{lang}.docx'
    with tempfile.NamedTemporaryFile('w', suffix='.md', dir=MS, delete=False) as tmp:
        tmp.write(md)
    try:
        subprocess.run(PANDOC + [tmp.name, '-o', str(out)], cwd=MS, check=True)
    finally:
        Path(tmp.name).unlink()
    print(out.relative_to(MS.parent))


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'draft'
    lang = sys.argv[2] if len(sys.argv) > 2 else 'both'
    for m in (['draft', 'submission'] if mode == 'both' else [mode]):
        for l in (['JA', 'EN'] if lang == 'both' else [lang]):
            compile_one(m, l)
