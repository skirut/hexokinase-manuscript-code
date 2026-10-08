#!/usr/bin/env python3
"""HTMD/ACEMD production workflow for Hexokinase simulations.

Purpose
-------
Prepares, equilibrates, and runs molecular-dynamics simulations for
Hexokinase PDB inputs across the pH conditions defined in the script.
The workflow uses HTMD system preparation, 10 ns equilibration, and
300 ns production runs at 300 K.

Inputs
------
- Protein PDB files discovered from the configured apo-data directory.
- pH values defined in the script (7.0, 5.0, and 3.0 in the current
  workflow).
- HTMD/ACEMD environment and GPU execution support.

Outputs
-------
For each system/pH:
- build directory
- equilibration directory
- production directory
- ACEMD outputs and simulation logs
- completion/error markers used by the workflow

Run from
--------
Run from the repository root with the HTMD/ACEMD environment active
and ``HEXOKINASE_DATA_ROOT`` configured.

Usage
-----
``python scripts/simulation/hk_sim.py``

Dependencies
------------
HTMD, ACEMD, and ``scripts.utils.paths``.
"""

import os
import time
from htmd.ui import *
from htmd.protocols.equilibration_v3 import Equilibration
from htmd.protocols.production_v6 import Production
from datetime import datetime, timedelta
from pathlib import Path

from scripts.utils.paths import APO_DIR
class MolecularSimulation:
    def __init__(self, pdb_id, base_dir, log_file):
        self.pdb_id = pdb_id
        self.base_dir = Path(base_dir)
        self.log_file = Path(log_file)
        self.pdb_base = self.base_dir / f'{self.pdb_id}'

    def log(self, message):
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        line = f"{timestamp}: {message}\n"
        with open(self.log_file, 'a') as log:
            log.write(line)
        print(line, end='')

    def _done_flag_path(self, pH):
        # Use a dotfile inside the pdb+pH folder to mark completion
        folder = self.base_dir / f'{self.pdb_id}_pH{pH}'
        return folder / '.simulation_done'

    def is_done(self, pH):
        return self._done_flag_path(pH).exists()

    def mark_done(self, pH):
        flag = self._done_flag_path(pH)
        flag.parent.mkdir(parents=True, exist_ok=True)
        # atomic write
        tmp = flag.with_suffix('.tmp')
        with open(tmp, 'w') as f:
            f.write(f"Done: {datetime.now().isoformat()}\n")
        tmp.replace(flag)

    def prepare_system(self, pH):
        self.log(f"Starting system preparation for {self.pdb_id} at pH {pH}..")
        pdb_path = self.base_dir / f'{self.pdb_id}.pdb'
        if not pdb_path.exists():
            raise FileNotFoundError(f"PDB file not found: {pdb_path}")
        protein = Molecule(str(pdb_path))
        protein.filter('protein')
        protein_op = systemPrepare(protein, pH=pH)
        protein_seg = autoSegment(protein_op)
        protein_solv = solvate(protein_seg, pad=10)

        amber_dir = self.base_dir / f'{self.pdb_id}_pH{pH}' / 'build-amber'
        amber_dir = str(amber_dir)
        amber.build(protein_solv, outdir=amber_dir)
        return amber_dir

    def equilibrate(self, build_dir, pH):
        self.log(f"Starting equilibration for {self.pdb_id} at pH {pH}...")
        equil_dir = os.path.join(self.base_dir, f'{self.pdb_id}_pH{pH}', 'equil')
        md = Equilibration()
        md.runtime = 20  # ns
        md.timeunits = 'ns'
        # md.timestep = 2
        md.temperature = 300
        md.useconstantratio = False
        md.write(build_dir, equil_dir)
        local = LocalGPUQueue()
        local.submit(equil_dir)
        local.wait()
        return equil_dir

    def production(self, equil_dir, pH):
        self.log(f"Starting production run for {self.pdb_id} at pH {pH}...")
        prod_dir = os.path.join(self.base_dir, f'{self.pdb_id}_pH{pH}', 'prod')
        md = Production()
        md.runtime = 300  # ns
        md.timeunits = 'ns'
        # md.timestep = 2
        md.temperature = 300
        md.acemd.bincoordinates = 'output.coor'
        md.acemd.extendedsystem = 'output.xsc'
        md.write(equil_dir, prod_dir)
        local = LocalGPUQueue()
        local.submit(prod_dir)
        local.wait()
        return prod_dir

    def run(self):
        pH_levels = [7.0, 5.0, 3.0]
        total_start_time = time.time()

        for i, pH in enumerate(pH_levels):
            # Skip if already completed
            if self.is_done(pH):
                self.log(f"SKIPPING {self.pdb_id} at pH {pH} — already completed (flag present). ({i+1}/{len(pH_levels)})")
                continue

            self.log(f"\nStarting simulations for {self.pdb_id} at pH {pH} ({i+1}/{len(pH_levels)})...")
            step_start_time = time.time()

            try:
                # ensure pdb-specific dir exists
                os.makedirs(os.path.join(self.base_dir, f'{self.pdb_id}_pH{pH}'), exist_ok=True)

                build_dir = self.prepare_system(pH)
                equil_dir = self.equilibrate(build_dir, pH)
                self.production(equil_dir, pH)

                # if we reach here production finished without exception -> mark done
                self.mark_done(pH)
                self.log(f"Marked {self.pdb_id} at pH {pH} as completed.")

            except Exception as e:
                # log error, continue with next pH / pdb
                self.log(f"ERROR while running {self.pdb_id} at pH {pH}: {repr(e)}")
                # optionally: write a small error file for debugging
                errfile = self.base_dir / f'{self.pdb_id}_pH{pH}' / 'error.txt'
                errfile.parent.mkdir(parents=True, exist_ok=True)
                with open(errfile, 'a') as ef:
                    ef.write(f"{datetime.now().isoformat()} - {repr(e)}\n")
                # do not mark done; next run will retry this pH
                continue

            step_duration = time.time() - step_start_time
            steps_remaining = len(pH_levels) - (i + 1)
            estimated_remaining_time = steps_remaining * step_duration
            estimated_completion_time = datetime.now() + timedelta(seconds=estimated_remaining_time)
            self.log(f"Completed simulations for {self.pdb_id} at pH {pH} in {step_duration:.2f} seconds.")
            self.log(f"Estimated time remaining: {timedelta(seconds=estimated_remaining_time)}. "
                     f"Estimated completion time: {estimated_completion_time.strftime('%Y-%m-%d %H:%M:%S')}\n")

        total_duration = time.time() - total_start_time
        self.log(f"All simulations (attempted) for {self.pdb_id} completed in {total_duration:.2f} seconds.\n")


def get_pdb_files_from_directory(directory):
    return [f for f in os.listdir(directory) if f.endswith('.pdb')]

if __name__ == '__main__':
    folder_path = APO_DIR
    log_file = os.path.join(folder_path, 'simulation_log.txt')
    pdb_files = get_pdb_files_from_directory(folder_path)
    pdb_ids = [os.path.splitext(f)[0] for f in pdb_files]

    for pdb_id in pdb_ids:
        simulation = MolecularSimulation(pdb_id, folder_path, log_file)
        simulation.run()


