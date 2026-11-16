#!/usr/bin/env python3
"""Generate compact mapping-order cases from the pinned Python implementation."""
import json
from pathlib import Path

from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from DockQ import DockQ as dq

ROOT = Path(__file__).resolve().parent
CASES = [
    {"name": "model_surplus_homomers", "model": ["A", "A", "C"], "native": ["A", "A"]},
    {"name": "native_surplus_reverse", "model": ["A", "A"], "native": ["A", "A", "A"]},
    {"name": "allowed_mismatch", "model": ["AC", "AA"], "native": ["AD", "AA"], "allowed_mismatches": 1},
    {"name": "fixed_pair_and_subset", "model": ["A", "A", "C", "X"], "native": ["A", "A", "C"], "fixed": [[2, 2]], "model_selected": [0, 1, 2], "native_selected": [0, 1, 2]},
    {"name": "hetero_identity", "model": [{"sequence": "A"}, {"sequence": "HEM", "hetero": "HEM"}, {"sequence": "LIG", "hetero": "LIG"}], "native": [{"sequence": "A"}, {"sequence": "HEM", "hetero": "HEM"}]},
    {"name": "incompatible_chain", "model": [{"sequence": "HEM", "hetero": "HEM"}], "native": [{"sequence": "LIG", "hetero": "LIG"}]},
    {"name": "no_injective_assignment", "model": ["A", "X"], "native": ["A", "A"]},
]


def make_structure(sequences, prefix):
    structure = {}
    for i, item in enumerate(sequences):
        spec = {"sequence": item} if isinstance(item, str) else item
        chain_id = f"{prefix}{i}"
        chain = Chain(chain_id)
        chain.sequence = spec["sequence"]
        chain.is_het = spec.get("hetero")
        for j, _ in enumerate(chain.sequence):
            chain.add(Residue((" ", j + 1, " "), "ALA", ""))
        structure[chain_id] = chain
    return structure


def run(case_index, case):
    model = make_structure(case["model"], f"M{case_index}C")
    native = make_structure(case["native"], f"N{case_index}C")
    model_ids = list(model)
    native_ids = list(native)
    selected_model = case.get("model_selected", list(range(len(model_ids))))
    selected_native = case.get("native_selected", list(range(len(native_ids))))
    fixed = [(native_ids[n], model_ids[m]) for n, m in case.get("fixed", [])]
    fixed_native = {n for n, _ in fixed}
    fixed_model = {m for _, m in fixed}
    model_combo = [model_ids[i] for i in selected_model if model_ids[i] not in fixed_model]
    native_combo = [native_ids[i] for i in selected_native if native_ids[i] not in fixed_native]
    dq.align_chains.cache_clear()
    try:
        clusters, reverse = dq.group_chains(
            model, native, model_combo, native_combo,
            case.get("allowed_mismatches", 0),
        )
        mappings = list(dq.get_all_chain_maps(
            clusters, dict(fixed), reverse, model_combo, native_combo
        ))
        rows = [[[native_ids.index(n), model_ids.index(m)] for n, m in row.items()] for row in mappings]
        result = {"ok": True, "mappings": rows}
    except SystemExit:
        result = {"ok": False, "mappings": []}
    return {**case, "expected": result}


output = {"schema_version": 1, "cases": [run(i, case) for i, case in enumerate(CASES)]}
(ROOT / "mapping_cases.json").write_text(json.dumps(output, indent=2) + "\n")
print(f"Wrote {len(output['cases'])} mapping cases")
