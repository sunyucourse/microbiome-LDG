"""Minimal NCBI taxdump parser: map tax_id -> class-level lineage info."""
import numpy as np
import pandas as pd


def load_taxonomy(taxdir='/workspace/data/external/ncbi_taxonomy'):
    """Parse nodes.dmp + names.dmp. Returns (parent, rank, name) dicts keyed by tax_id."""
    parent, rank = {}, {}
    with open(f'{taxdir}/nodes.dmp') as fh:
        for line in fh:
            parts = [p.strip() for p in line.split('|')]
            tid = int(parts[0])
            parent[tid] = int(parts[1])
            rank[tid] = parts[2]
    name = {}
    with open(f'{taxdir}/names.dmp') as fh:
        for line in fh:
            parts = [p.strip() for p in line.split('|')]
            if parts[3] == 'scientific name':
                name[int(parts[0])] = parts[1]
    return parent, rank, name


def classify_to_class(tax_id, parent, rank, name, max_depth=60):
    """Walk up the tree; return (class_name, phylum_name, kingdom_name, scientific_name)."""
    cls = phy = king = None
    tid = tax_id
    seen = set()
    depth = 0
    while tid in parent and tid != 1 and tid not in seen and depth < max_depth:
        seen.add(tid)
        r = rank.get(tid, '')
        if r == 'class' and cls is None:
            cls = name.get(tid)
        elif r == 'phylum' and phy is None:
            phy = name.get(tid)
        elif r in ('kingdom', 'superkingdom') and king is None:
            king = name.get(tid)
        tid = parent[tid]
        depth += 1
    return cls, phy, king, name.get(tax_id)


ENDOTHERM_CLASSES = {'Mammalia', 'Aves'}


def thermo_group(class_name):
    if class_name in ENDOTHERM_CLASSES:
        return 'endotherm'
    if class_name is None:
        return None
    return 'ectotherm'
