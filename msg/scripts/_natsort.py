"""Natural string sort key — port of the in-script function from
msg/pull_thin_tsv.py (Connelly Barnes / Seo Sanghyeon).

Case-sensitive (matches legacy behavior).
"""
import re


def natsort_key(s):
    return [int(t) if t.isdigit() else t
            for t in re.findall(r'\d+|\D+', s)]
