"""Pathway name abbreviation utilities for GlycoEnzOnto pathway names.

Standard glycobiology abbreviations applied in order of priority.
Used for figure axis labels and manuscript text.
"""

PATHWAY_ABBREVIATIONS = [
    ("glycosaminoglycan", "GAG"),
    ("glycosphingolipid", "GSL"),
    ("glycosylphosphatidylinositol", "GPI"),
    ("uridine-5'-diphosphate", "UDP"),
    ("guanosine-5'-diphosphate", "GDP"),
    ("cytidine-5'-monophosphate", "CMP"),
    ("phosphoadenosine phosphosulfate", "PAPS"),
    ("n-acetylgalactosamine", "GalNAc"),
    ("n-acetylglucosamine", "GlcNAc"),
    ("n-acetylneuraminic acid", "Neu5Ac"),
    ("biosynthetic pathway", "biosynthesis"),
    ("degradation pathway", "degradation"),
    ("synthetic pathway", "synthesis"),
    ("biosynthetic process", "biosynthesis"),
    (" pathway", ""),  # catch-all, applied last
]


def abbreviate_pathway(name: str) -> str:
    """Abbreviate a GlycoEnzOnto pathway name using standard glycobiology abbreviations.

    Parameters
    ----------
    name : str
        Original pathway name (may include surrounding quotes from GMT file).

    Returns
    -------
    str
        Abbreviated pathway name.
    """
    import re
    result = name.strip('"')
    for old, new in PATHWAY_ABBREVIATIONS:
        result = re.sub(re.escape(old), new, result, flags=re.IGNORECASE)
    return result.strip()
