"""Homology modeling workflow for Hexokinase-2 (HK2, PDB 2NZT).

Purpose
-------
Extracts the sequence from the HK2 PDB, aligns it with the supplied
target FASTA using Clustal Omega, creates a MODELLER PIR alignment while
accounting for missing residues, and builds a homology model.

Inputs
------
- ``2nzt.pdb``
- ``2nzt.fasta``
- Clustal Omega
- MODELLER

Outputs
-------
- ``2nzt_pdb_seq.fasta``
- ``2nzt_clustalo.aln``
- ``2nzt.alignment.ali``
- ``2nzt_model/``

Run from
--------
Run from the directory containing the HK2 PDB and FASTA inputs.

Usage
-----
``python 2nzt_modeller.py``

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
    subprocess.run(['clustalo', '-i', combined_fasta, '--outfmt=clu', '--force', '-o', output_file], check=True)
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
            for res in chain:
                if res.id[0] == ' ':
                    return res.id[1], chain.id
    raise ValueError(f"No suitable residues found in {pdb_file}")

def parse_remark_465_missing_residues(pdb_file, target_chain='A'):
    missing_residues = []
    with open(pdb_file) as f:
        for line in f:
            if line.startswith('REMARK 465') and 'MISSING RESIDUES' not in line:
                m = re.match(r'REMARK 465\s+\w+\s+([A-Z]+)\s+([A-Z])\s+(\d+)', line)
                if m:
                    res_name, chain_id, res_num = m.groups()
                    if chain_id == target_chain:
                        missing_residues.append(int(res_num))
    return sorted(missing_residues)

def write_pir_alignment(struct_label, pdb_file, pdb_seq, fasta_label, fasta_seq, output_file, start_resnum, chain_id='A'):
    pdb_name = os.path.basename(pdb_file).split('.')[0]
    first_resnum, _ = get_first_residue_number_and_chain(pdb_file)
    gap_count = first_resnum - start_resnum

    padded_pdb_seq = '-' * gap_count + pdb_seq
    padded_fasta_seq = '-' * gap_count + fasta_seq
    total_residues = len(padded_pdb_seq.replace('-', '')) + gap_count

    with open(output_file, 'w') as f:
        f.write(f'>P1;{struct_label}\n')
        f.write(f'structureX:{pdb_name}.pdb:{start_resnum}:{chain_id}:+{total_residues}:{chain_id}:::-1.00:-1.00\n')
        f.write(format_sequence(padded_pdb_seq) + '*\n\n')
        f.write(f'>P1;{fasta_label}\n')
        f.write(f'sequence:{fasta_label}:::::::0.00:0.00\n')
        f.write(format_sequence(padded_fasta_seq) + '*\n')

def build_model_with_modeller(alignment_file, knowns_label, sequence_label, renumber_start):
    env = Environ()

    class MyModel(automodel):
        def special_patches(self, aln):
            if renumber_start:
                self.rename_segments(segment_ids=['A'], renumber_residues=[renumber_start -2])

    a = MyModel(env,
                alnfile=os.path.abspath(alignment_file),
                knowns=knowns_label,
                sequence=sequence_label,
                assess_methods=(assess.DOPE, assess.GA341))
    a.starting_model = 1
    a.ending_model = 1
    # Removed a.renumber_residues = renumber_start to avoid conflict with rename_segments()
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
    if len(seqs) < 2:
        print(f"Alignment parsing failed for {alignment_file}")
        return

    labels = list(seqs)
    pdb_label, fasta_label = labels[0], labels[1]
    pdb_aln_seq = seqs[pdb_label]
    fasta_aln_seq = seqs[fasta_label]

    first_resnum, chain_id = get_first_residue_number_and_chain(pdb_file)
    missing_residues = parse_remark_465_missing_residues(pdb_file, target_chain=chain_id)

    if missing_residues and min(missing_residues) < first_resnum:
        start_resnum = min(missing_residues)
        print(f"Using missing residue renumbering start: {start_resnum}")
    else:
        start_resnum = first_resnum
        print(f"Using default start from PDB: {start_resnum}")

    pir_out = pdb_file.replace('.pdb', '.alignment.ali')
    write_pir_alignment(pdb_label, pdb_file, pdb_aln_seq, fasta_label, fasta_aln_seq, pir_out, start_resnum, chain_id)
    print(f"Generated PIR alignment: {pir_out}")

    build_model_with_modeller(pir_out, pdb_label, fasta_label, start_resnum)

def main():
    pdb_file = '2nzt.pdb'
    fasta_file = '2nzt.fasta'
    if not os.path.exists(pdb_file) or not os.path.exists(fasta_file):
        print("Required files 2nzt.pdb or 2nzt.fasta not found.")
        return
    process_pair(pdb_file, fasta_file)

if __name__ == '__main__':
    main()

