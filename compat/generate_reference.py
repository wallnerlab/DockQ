#!/usr/bin/env python3
"""Generate versioned scientific fixtures from the unmodified Python reference.

Run in an environment with this checkout's compiled DockQ installed. Deliberately
does not import the local source tree, whose Cython extension may not be built.
"""
import argparse
import importlib.metadata
import json
import platform
import subprocess
from pathlib import Path

import numpy as np
from DockQ import DockQ as dq
from DockQ import operations, operations_nocy

ROOT = Path(__file__).resolve().parents[1]


def clean(value):
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    if isinstance(value, np.ndarray):
        return clean(value.tolist())
    if isinstance(value, np.generic):
        return value.item()
    return value


def structure_summary(structure):
    return [{
        'id': chain.id, 'sequence': chain.sequence, 'hetero': chain.is_het,
        'residues': [
            {'name': r.resname, 'number': r.id[1], 'insertion_code': r.id[2].strip(),
             'hetero': r.id[0] != ' ', 'atom_names': [a.id for a in r],
             'first_coord': clean(next(r.get_atoms()).coord)}
            for r in chain
        ],
    } for chain in structure]


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'compat/reference.json')
    args = parser.parse_args()
    cases = [
        ('1a2k', '1A2K_r_l_b.model.pdb', '1A2K_r_l_b.pdb', 'BAC:ABC', {}),
        ('1a2k_numbering', '1A2K_r_l_b.model.pdb', '1A2K_r_l_b.pdb', 'BAC:ABC', {'no_align': True}),
        ('dimer', 'dimer_dimer.model.pdb', 'dimer_dimer.pdb', 'ABLH:ABLH', {}),
        ('mismatch', 'model.pdb', 'native.pdb', 'AB:AB', {}),
        ('octamer', '1EXB_r_l_b.model.pdb', '1EXB_r_l_b.pdb', 'BACDFHEG:ABDCEGFH', {}),
        ('cif', '1EXB_r_l_b.model.pdb', '1EXB.cif.gz', 'DH:AE', {}),
        ('peptide', '6qwn-assembly1.cif.gz', '6qwn-assembly2.cif.gz', 'AF:BG', {'capri_peptide': True}),
        ('ligand', '1HHO_hem.cif', '2HHB_hem.cif', 'ABDCF:ABEFG', {'small_molecule': True}),
        ('nucleic', '../rust/tests/data/nucleic.model.pdb', '../rust/tests/data/nucleic.native.pdb', 'AB:AB', {}),
    ]
    output = {'schema_version': 1, 'baseline_revision': subprocess.check_output(
        ['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        'python': platform.python_version(),
        'dependencies': {p: importlib.metadata.version(p) for p in
                         ['dockq', 'numpy', 'biopython', 'networkx', 'parallelbar', 'Cython']},
        'operations_module': str(operations.__file__), 'structures': {}, 'cases': []}
    for name, mf, nf, mapping, opts in cases:
        small = opts.get('small_molecule', False)
        model, native = [dq.load_PDB(str(ROOT/'examples'/f), small_molecule=small) for f in (mf,nf)]
        for f, s in [(mf,model),(nf,native)]:
            output['structures'][f] = structure_summary(s)
        if mapping is None:
            # Characterize a fixed assignment without factorial search. Select matching
            # chain IDs where possible; ligand author/label semantics stay observable.
            ids = {c.id for c in model}
            cm = {c.id:c.id for c in native if c.id in ids}
        else:
            mm, nn = mapping.split(':')
            cm = dict(zip(nn, mm))
        score_opts = {k:v for k,v in opts.items() if k != 'small_molecule'}
        record = {'name': name, 'model': mf, 'native': nf,
                  'mapping': list(cm.items()), 'options': opts}
        try:
            results, total = dq.run_on_all_native_interfaces(model,native,chain_map=cm,**score_opts)
            record.update(results=clean(results), total_dockq=total,
                          global_dockq=total/len(results) if results else None)
        except Exception as exc:
            record['reference_error'] = f'{type(exc).__name__}: {exc}'
        output['cases'].append(record)
    # A controlled example exposes Cython's early-exit behavior independently of
    # floating rounding: the second atom would be a contact under exhaustive search.
    a = np.array([[100.,0.,0.],[0.,0.,0.]],dtype=np.float32)
    b = np.array([[0.,0.,0.]],dtype=np.float32)
    ga, gb = np.array([2],dtype=int), np.array([1],dtype=int)
    output['kernel_discrepancy'] = {
        'a': clean(a), 'b': clean(b), 'groups_a': [2], 'groups_b': [1],
        'compiled': clean(np.asarray(operations.residue_distances(a,b,ga,gb))),
        'fallback': clean(operations_nocy.residue_distances(a,b,ga,gb))}
    # Tie-breaking and numbering are tested with actual Bio chains, not a reimplementation.
    from Bio.PDB.Chain import Chain
    from Bio.PDB.Residue import Residue
    output['alignments'] = []
    for sa,sb in [('AAAA','AAA'),('ACDE','ACE'),('ACX','ACY'),('ABAB','BABA'),('M','M')]:
        ca, cb = Chain('M'), Chain('N')
        ca.sequence, cb.sequence = sa,sb
        for c,s in [(ca,sa),(cb,sb)]:
            for i in range(len(s)):
                c.add(Residue((' ',i-2,' '),'ALA',''))
        for numbering in (False,True):
            # Standalone Bio chains with the same IDs compare equal, even when
            # sequences differ; do not reuse align_chains' production LRU here.
            dq.align_chains.cache_clear()
            aln = dq.format_alignment(dq.align_chains(ca,cb,use_numbering=numbering))
            output['alignments'].append({'model':sa,'native':sb,'numbering':numbering,**aln})
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(clean(output),indent=2,allow_nan=False)+'\n')
    print(f'Wrote {args.output}: {len(output["cases"])} cases, {len(output["structures"])} structures')


if __name__ == '__main__':
    main()
