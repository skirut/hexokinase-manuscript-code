"""HTMD-based G6P-bound Hexokinase system preparation/simulation workflow.

Purpose
-------
Provides the HTMD workflow used to prepare and run Hexokinase systems
containing G6P, including equilibration and production stages. The
module contains helper logic for identifying PDB inputs and uses the
configured study data directories.

Inputs
------
- Protein PDB files from the configured bound-data directory.
- Parameterized G6P ligand files under the configured ligand directory.
- HTMD/ACEMD installation and license/environment configuration.

Outputs
-------
Creates system-specific build, equilibration, and production
directories under the configured data tree, together with simulation
logs and ACEMD output files.

Run from
--------
Run from the repository root with the HTMD/ACEMD environment active and
``HEXOKINASE_DATA_ROOT`` configured.

Usage
-----
Run the script according to the input PDB layout expected by the module;
inspect the main block in the script before launching production jobs.

Dependencies
------------
HTMD, ACEMD, NumPy, and the repository path utilities.
"""

import os
import time
import numpy as np
from datetime import datetime, timedelta

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

    def prepare_system(self, pH):
        self.log(f"Starting system preparation for {self.pdb_id} at pH {pH}")

        # Load protein
        protein_path = os.path.join(self.base_dir, f"{self.pdb_id}.pdb")
        protein = Molecule(protein_path)

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

        # pH formatting
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
        ligand.set("segid", "LIG")

        # Place ligand
        if crystal_coords is not None and ligand.numAtoms == crystal_coords.shape[0]:
            ligand.set("coords", crystal_coords)
        elif crystal_coords is not None:
            ligand.moveBy(crystal_coords.mean(axis=0) - ligand.get("coords").mean(axis=0))

        self.log("Placed parameterized G6P in binding site")

        # Combine protein + ligand
        system = protein.copy()
        system.append(ligand)

        # Solvate
        system = solvate(system, pad=10)

        # Amber build
        topos = amber.defaultTopo() + [ligand_mol2]
        params = amber.defaultParam() + [ligand_frcmod]

        build_dir = os.path.join(self.base_dir, f"{self.pdb_id}_pH{pH}", "build-amber")
        os.makedirs(build_dir, exist_ok=True)

        amber.build(system, topo=topos, param=params, outdir=build_dir)
        self.log(f"Amber build completed: {build_dir}")

        return build_dir

    def equilibrate(self, build_dir, pH):
        self.log(f"Starting equilibration for {self.pdb_id} at pH {pH}")

        equil_dir = os.path.join(self.base_dir, f"{self.pdb_id}_pH{pH}", "equil")

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

    def production(self, equil_dir, pH):
        self.log(f"Starting production for {self.pdb_id} at pH {pH}")

        prod_dir = os.path.join(self.base_dir, f"{self.pdb_id}_pH{pH}", "prod")

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


# ------------------------
# Main execution
# ------------------------

def get_pdb_ids(directory):
    return [os.path.splitext(f)[0] for f in os.listdir(directory) if f.endswith(".pdb")]

folder_path = BOUND_DIR
log_file = os.path.join(folder_path, "simulation_log.txt")

pdb_ids = get_pdb_ids(folder_path)

for pdb_id in pdb_ids:
    sim = MolecularSimulation(pdb_id, folder_path, log_file)
    sim.run()


