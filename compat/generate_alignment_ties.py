#!/usr/bin/env python3
"""Generate small exhaustive alignment tie fixtures from Biopython."""
import itertools
import json
from pathlib import Path

from Bio import Align

ROOT = Path(__file__).resolve().parent
sequences = [
    "".join(chars)
    for size in range(1, 5)
    for chars in itertools.product("AC", repeat=size)
]
aligner = Align.PairwiseAligner()
aligner.match_score = 5
aligner.mismatch_score = 0
aligner.open_gap_score = -4
aligner.extend_gap_score = -0.5

records = []
for model in sequences:
    for native in sequences:
        alignment = aligner.align(model, native)[0]
        aligned_model = alignment[0, :]
        aligned_native = alignment[1, :]
        matches = "".join(
            "|" if a == b else " " if "-" in (a, b) else "."
            for a, b in zip(aligned_model, aligned_native)
        )
        records.append({
            "model": model,
            "native": native,
            "model_aligned": aligned_model,
            "native_aligned": aligned_native,
            "matches": matches,
        })

(ROOT / "alignment_ties.json").write_text(
    json.dumps({"schema_version": 1, "alignments": records}, indent=2) + "\n"
)
print(f"Wrote {len(records)} exhaustive alignment cases")
