"""HTMD/ACEMD workflow for G6P-bound mutant and modeled Hexokinase systems.

Purpose
-------
Prepares G6P-bound systems by removing/repositioning crystallographic
G6P, inserting the pH-specific parameterized G6P ligand, solvating the
system, building Amber topology, equilibrating for 10 ns, and running
300 ns production.

Inputs
------
- Protein PDB files under the configured bound-data directory.
- pH-specific ACPYPE G6P files under
  ``$HEXOKINASE_DATA_ROOT/ligand_regulatory_site``.
- HTMD/ACEMD installation.

Outputs
-------
For each system/pH:
- ``build-amber/``
- ``equil/``
- ``prod/``
- simulation log/error files and ACEMD outputs

Run from
--------
Run from the repository root with the HTMD/ACEMD environment active.

Usage
-----
``python scripts/structure/mutations/mutation_g6p.py``

Dependencies
------------
HTMD, ACEMD, NumPy, Amber parameter files generated for G6P, and
``scripts.utils.paths``.
"""

import os
import numpy as np
from datetime import datetime

from scripts.utils.paths import BOUND_DIR, LIGAND_DIR
os.environ["HTMD_LICENSE"] = "skip"

from htmd.ui import *
from htmd.protocols.equilibration_v3 import Equilibration
from htmd.protocols.production_v6 import Production


class MolecularSimulation:
    def __init__(self, pdb_id, base_dir, log_file):
        self.pdb_id = pdb_id
        self.base_dir = base_dir
        self.log_file = log_file

    def log(self, message):
        with open(self.log_file, 'a') as log:
            log.write(f"{datetime.now()} : {message}\n")
        print(message)

    # ----------------------------------------------------
    # SYSTEM PREPARATION
    # ----------------------------------------------------
    def prepare_system(self, pH):
        self.log(f"Starting system preparation for {self.pdb_id} at pH {pH}")

        protein_path = os.path.join(self.base_dir, f"{self.pdb_id}.pdb")
        protein = Molecule(protein_path)

        # FIX 1: Assign segid to protein
        protein.set("segid", "PROT")

        # Extract crystallographic G6P coordinates
        crystal = protein.copy()
        sel_crystal = crystal.atomselect("resname G6P")

        if np.sum(sel_crystal) == 0:
            self.log("WARNING: No crystal G6P found; ligand will be placed approximately")
            crystal_coords = None
        else:
            crystal_coords = crystal.get("coords")[sel_crystal]

        # Remove crystallographic ligand
        protein.remove("resname G6P")
        self.log("Removed crystallographic G6P ligand")

        # Format pH string
        pH_str = str(int(pH)) if pH == int(pH) else str(int(pH * 10))

        ligand_dir = LIGAND_DIR / f"g6p_pH{pH_str}.acpype"

        mol2_files = [f for f in os.listdir(ligand_dir) if f.endswith(".mol2")]
        frcmod_files = [f for f in os.listdir(ligand_dir) if f.endswith(".frcmod")]

        if not mol2_files or not frcmod_files:
            raise FileNotFoundError(f"Missing ACPYPE files in {ligand_dir}")

        ligand_mol2 = os.path.join(ligand_dir, mol2_files[0])
        ligand_frcmod = os.path.join(ligand_dir, frcmod_files[0])

        self.log(f"Using ligand mol2: {ligand_mol2}")
        self.log(f"Using ligand frcmod: {ligand_frcmod}")

        ligand = Molecule(ligand_mol2)

        # FIX 2: Assign segid to ligand
        ligand.set("segid", "LIG")

        # Place ligand at crystal position if possible
        if crystal_coords is not None and ligand.numAtoms == crystal_coords.shape[0]:
            ligand.set("coords", crystal_coords)
        elif crystal_coords is not None:
            ligand.moveBy(
                crystal_coords.mean(axis=0) - ligand.get("coords").mean(axis=0)
            )

        self.log("Placed parameterized G6P in binding site")

        # Combine protein + ligand
        system = protein.copy()
        system.append(ligand)

        # Solvate system
        system = solvate(system, pad=10)

        # Build Amber topology
        topos = amber.defaultTopo() + [ligand_mol2]
        params = amber.defaultParam() + [ligand_frcmod]

        build_dir = os.path.join(
            self.base_dir,
            f"{self.pdb_id}_pH{pH}",
            "build-amber"
        )
        os.makedirs(build_dir, exist_ok=True)

        amber.build(system, topo=topos, param=params, outdir=build_dir)

        self.log(f"Amber build completed: {build_dir}")

        return build_dir

    # ----------------------------------------------------
    # EQUILIBRATION
    # ----------------------------------------------------
    def equilibrate(self, build_dir, pH):
        self.log(f"Starting equilibration for {self.pdb_id} at pH {pH}")

        equil_dir = os.path.join(
            self.base_dir,
            f"{self.pdb_id}_pH{pH}",
            "equil"
        )

        md = Equilibration()
        md.runtime = 10
        md.timeunits = "ns"
        md.temperature = 300
        md.useconstantratio = False

        md.write(build_dir, equil_dir)

        queue = LocalGPUQueue()
        queue.submit(equil_dir)
        queue.wait()

        self.log("Equilibration completed")

        return equil_dir

    # ----------------------------------------------------
    # PRODUCTION
    # ----------------------------------------------------
    def production(self, equil_dir, pH):
        self.log(f"Starting production for {self.pdb_id} at pH {pH}")

        prod_dir = os.path.join(
            self.base_dir,
            f"{self.pdb_id}_pH{pH}",
            "prod"
        )

        md = Production()
        md.runtime = 300
        md.timeunits = "ns"
        md.temperature = 300
        md.acemd.bincoordinates = "output.coor"
        md.acemd.extendedsystem = "output.xsc"

        md.write(equil_dir, prod_dir)

        queue = LocalGPUQueue()
        queue.submit(prod_dir)
        queue.wait()

        self.log("Production completed")

    # ----------------------------------------------------
    # RUN FULL PIPELINE
    # ----------------------------------------------------
    def run(self):
        pH_levels = [3.0, 5.0, 7.0]

        for pH in pH_levels:
            run_dir = os.path.join(self.base_dir, f"{self.pdb_id}_pH{pH}")

            if os.path.exists(run_dir):
                self.log(f"Skipping {self.pdb_id} at pH {pH} (already exists)")
                continue

            try:
                os.makedirs(run_dir, exist_ok=True)
                build_dir = self.prepare_system(pH)
                equil_dir = self.equilibrate(build_dir, pH)
                self.production(equil_dir, pH)

            except Exception as e:
                self.log(f"ERROR at pH {pH}: {e}")


# ----------------------------------------------------
# MAIN EXECUTION
# ----------------------------------------------------

def get_pdb_ids(directory):
    return [
        os.path.splitext(f)[0]
        for f in os.listdir(directory)
        if f.endswith(".pdb")
    ]


folder_path = BOUND_DIR
log_file = os.path.join(folder_path, "simulation_log.txt")

pdb_ids = get_pdb_ids(folder_path)

for pdb_id in pdb_ids:
    sim = MolecularSimulation(pdb_id, folder_path, log_file)
    sim.run()

