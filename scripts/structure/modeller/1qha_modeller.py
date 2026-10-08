"""Homology modeling workflow for Hexokinase-1 (HK1, PDB 1QHA).

Purpose
-------
Extracts a sequence from the HK1 PDB, aligns it with the supplied
target FASTA using Clustal Omega, creates a MODELLER PIR alignment while
accounting for missing residues, and builds a homology model.

Inputs
------
- ``1qha.pdb``
- ``1qha.fasta``
- Clustal Omega executable in ``PATH``
- MODELLER installation

Outputs
-------
- ``1qha_pdb_seq.fasta``
- ``1qha_clustalo.aln``
- ``1qha.alignment.ali``
- ``1qha_model/`` containing MODELLER output

Run from
--------
Run from the directory containing ``1qha.pdb`` and ``1qha.fasta``.

Usage
-----
``python 1qha_modeller.py``

Dependencies
------------
Biopython, Clustal Omega, MODELLER, and Python.
"""

import os
import subprocess
import re
from Bio.PDB import PDBParser, PPBuilder
from modeller import Environ
from modeller.automodel import automodel, assess

def extract_sequence_from_pdb(pdb_file, output_fasta):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('', pdb_file)
    ppb = PPBuilder()
    seq = ''
    for pp in ppb.build_peptides(structure):
        seq += str(pp.get_sequence())
    with open(output_fasta, 'w') as f:
        f.write(f'>{os.path.basename(pdb_file).split(".")[0]}_pdb_seq\n')
        f.write(seq + '\n')
    return output_fasta

def run_clustalo(seq1_file, seq2_file, output_file):
    combined_fasta = seq1_file.replace('_pdb_seq.fasta', '_combined.fasta')
    with open(combined_fasta, 'w') as outf:
        for f in [seq1_file, seq2_file]:
            with open(f, 'r') as inf:
                outf.write(inf.read())
    cmd = ['clustalo', '-i', combined_fasta, '--outfmt=clu', '--force', '-o', output_file]
    subprocess.run(cmd, check=True)
    os.remove(combined_fasta)

def parse_clustal_alignment(aln_file):
    with open(aln_file) as f:
        lines = f.readlines()
    seqs = {}
    for line in lines:
        line = line.strip()
        if line and not line.startswith(('CLUSTAL', 'MUSCLE')) and not line.startswith(' '):
            parts = line.split()
            if len(parts) >= 2:
                label, seqpart = parts[0], parts[1]
                seqs.setdefault(label, '')
                seqs[label] += seqpart
    return seqs

def format_sequence(seq, line_length=60):
    return '\n'.join(seq[i:i+line_length] for i in range(0, len(seq), line_length))

def get_first_residue_number_and_chain(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('X', pdb_file)
    for model in structure:
        for chain in model:
            if chain.id == 'A':
                for res in chain:
                    if res.id[0] == ' ':
                        return res.id[1], chain.id
    raise ValueError(f"No suitable residues found in chain A of {pdb_file}")

def parse_remark_465_missing_residues(pdb_file):
    missing_residues = {}
    with open(pdb_file) as f:
        for line in f:
            if line.startswith('REMARK 465') and len(line) > 30 and line[10:15].strip() != 'M' and line[15:18].strip() != '':
                parts = line[11:].split()
                if len(parts) >= 4:
                    resname, chain_id, resnum = parts[1], parts[2], parts[3]
                    try:
                        resnum = int(resnum)
                        if chain_id not in missing_residues:
                            missing_residues[chain_id] = []
                        missing_residues[chain_id].append(resnum)
                    except ValueError:
                        continue
    for c in missing_residues:
        missing_residues[c].sort()
    return missing_residues

def write_pir_alignment(struct_label, pdb_file, pdb_seq, fasta_label, fasta_seq, output_file, missing_residues):
    first_resnum, chain_id = get_first_residue_number_and_chain(pdb_file)
    min_missing = min(missing_residues.get('A', []), default=None)
    if min_missing is not None and min_missing < first_resnum:
        start_resnum = min_missing
    else:
        start_resnum = first_resnum
    n_gaps = max(0, first_resnum - start_resnum)
    pdb_seq_with_gaps = '-' * n_gaps + pdb_seq
    total_residues = len(pdb_seq_with_gaps.replace('-', ''))
    pdb_name = os.path.basename(pdb_file).split('.')[0]
    with open(output_file, 'w') as f:
        f.write(f'>P1;{struct_label}\n')
        f.write(f'structureX:{pdb_name}.pdb:{start_resnum}:{chain_id}:+{total_residues}:{chain_id}:::-1.00:-1.00\n')
        f.write(format_sequence(pdb_seq_with_gaps) + '*\n\n')
        f.write(f'>P1;{fasta_label}\n')
        f.write(f'sequence:{fasta_label}:::::::0.00:0.00\n')
        f.write(format_sequence(fasta_seq) + '*\n')

def build_model_with_modeller(alignment_file, knowns_label, sequence_label):
    print(f"Building model from alignment file: {alignment_file}")
    pdb_name = os.path.basename(alignment_file).split('.')[0].replace('_alignment', '')
    output_dir = f"{pdb_name}_model"
    os.makedirs(output_dir, exist_ok=True)
    env = Environ()
    env.io.output_directory = output_dir
    a = automodel(env,
                  alnfile=os.path.abspath(alignment_file),
                  knowns=knowns_label,
                  sequence=sequence_label,
                  assess_methods=(assess.DOPE, assess.GA341))
    a.starting_model = 1
    a.ending_model = 1
    a.make()

def process_pair(pdb_file, fasta_file):
    print(f"Processing: {pdb_file} and {fasta_file}")
    pdb_seq_file = pdb_file.replace('.pdb', '_pdb_seq.fasta')
    extract_sequence_from_pdb(pdb_file, pdb_seq_file)
    alignment_file = pdb_file.replace('.pdb', '_clustalo.aln')
    try:
        run_clustalo(pdb_seq_file, fasta_file, alignment_file)
    except subprocess.CalledProcessError:
        print(f"Clustal Omega failed for {pdb_file}")
        return
    seqs = parse_clustal_alignment(alignment_file)
    labels = list(seqs.keys())
    if len(labels) < 2:
        print(f"Alignment parsing failed for {alignment_file}")
        return
    pdb_label, fasta_label = labels[0], labels[1]
    pdb_aln_seq = seqs[pdb_label]
    fasta_aln_seq = seqs[fasta_label]
    missing_residues = parse_remark_465_missing_residues(pdb_file)
    pir_out = pdb_file.replace('.pdb', '.alignment.ali')
    write_pir_alignment(pdb_label, pdb_file, pdb_aln_seq, fasta_label, fasta_aln_seq, pir_out, missing_residues)
    print(f"Generated PIR alignment: {pir_out}")
    build_model_with_modeller(pir_out, pdb_label, fasta_label)

def main():
    files = [f for f in os.listdir('.') if f.endswith('.pdb')]
    for pdb_file in files:
        fasta_file = pdb_file.replace('.pdb', '.fasta')
        if not os.path.exists(fasta_file):
            print(f"Skipping {pdb_file}: Missing {fasta_file}")
            continue
        process_pair(pdb_file, fasta_file)

if __name__ == '__main__':
    main()

