#!/usr/bin/env python3
"""HTMD/ACEMD workflow for Hexokinase wild-type and mutant simulations.

Purpose
-------
Runs or resumes Hexokinase molecular-dynamics simulations for explicit
PDB targets. The workflow prepares systems, equilibrates them, runs
300 ns production, and records completion in both per-run flags and a
JSON ledger.

Inputs
------
- PDB files under the configured apo-data directory.
- Optional command-line targets and pH values.
- HTMD/ACEMD execution environment.

Outputs
-------
- ``build-amber/``, ``equil/``, and ``prod/`` directories per system/pH
- ``.simulation_done`` completion flags
- ``completed_runs.json`` ledger
- simulation logs and error files

Run from
--------
Run from the repository root with the HTMD/ACEMD environment active.

Usage
-----
``python scripts/simulation/mutants_sim.py --help``

The script supports ``--base``, ``--log``, ``--ph``, and ``--targets``.

Dependencies
------------
HTMD, ACEMD, Python standard-library modules, and
``scripts.utils.paths``.
"""

import os
import time
import json
import argparse
from datetime import datetime, timedelta
from pathlib import Path

from htmd.ui import *                     # brings in Molecule, systemPrepare, autoSegment, solvate, amber, LocalGPUQueue
from htmd.protocols.equilibration_v3 import Equilibration
from htmd.protocols.production_v6 import Production

from scripts.utils.paths import APO_DIR
LEDGER = 'completed_runs.json'

# ---------- Utilities ----------
def norm_pdb_id_to_file(base_dir: Path, pdb_id_or_file: str) -> tuple[str, Path]:
    """
    Return the run prefix (no .pdb) and absolute PDB path.
    Accepts '1qha', '1qha.pdb', '3hm8_T546M', '3hm8_T546M.pdb', etc.
    """
    p = Path(pdb_id_or_file)
    runprefix = p.stem if p.suffix.lower() == '.pdb' else pdb_id_or_file
    pdbpath = base_dir / f"{runprefix}.pdb"
    return runprefix, pdbpath

def good_file(p: Path, min_bytes=100) -> bool:
    try:
        return p.exists() and p.stat().st_size > min_bytes
    except Exception:
        return False

# ---------- Simulation Driver ----------
class MolecularSimulation:
    def __init__(self, runprefix: str, base_dir: Path, log_file: Path):
        self.runprefix = runprefix          # e.g., '1qha', '3hm8_T546M'
        self.base_dir = base_dir
        self.log_file = log_file

    def log(self, msg: str):
        ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        line = f"{ts}: {msg}\n"
        with open(self.log_file, 'a') as f:
            f.write(line)
        print(line, end='')

    def _root(self, pH):            # e.g., /.../1qha_pH7.0
        return self.base_dir / f"{self.runprefix}_pH{pH}"

    def _build_dir(self, pH):
        return self._root(pH) / "build-amber"

    def _equil_dir(self, pH):
        return self._root(pH) / "equil"

    def _prod_dir(self, pH):
        return self._root(pH) / "prod"

    def _prod_paths(self, pH):
        d = self._prod_dir(pH)
        return {
            "dir": d,
            "coor": d / "output.coor",
            "xsc":  d / "output.xsc",
            "log":  d / "acemd.log",
        }

    def _done_flag(self, pH):
        return self._root(pH) / ".simulation_done"

    def _read_ledger(self) -> dict:
        path = self.base_dir / LEDGER
        try:
            if path.exists():
                return json.load(open(path))
        except Exception:
            return {}
        return {}

    def _write_ledger_entry(self, pH):
        path = self.base_dir / LEDGER
        data = self._read_ledger()
        key = f"{self.runprefix}_pH{pH}"
        data[key] = {"done_at": datetime.now().isoformat()}
        with open(path, "w") as f:
            json.dump(data, f, indent=2)

    def is_done(self, pH) -> bool:
        """A run is 'done' if .simulation_done exists OR ledger has it OR prod outputs look valid."""
        flag = self._done_flag(pH)
        ledger = self._read_ledger()
        key = f"{self.runprefix}_pH{pH}"
        if flag.exists() or key in ledger:
            return True
        paths = self._prod_paths(pH)
        return any(good_file(paths[k]) for k in ["coor", "xsc", "log"])

    def mark_done(self, pH):
        flag = self._done_flag(pH)
        flag.parent.mkdir(parents=True, exist_ok=True)
        tmp = flag.with_suffix('.tmp')
        tmp.write_text(f"Done: {datetime.now().isoformat()}\n")
        tmp.replace(flag)
        self._write_ledger_entry(pH)

    # ---------- Pipeline steps (idempotent where possible) ----------
    def prepare_system(self, pdb_path: Path, pH: float) -> str:
        self.log(f"[{self.runprefix} pH {pH}] Preparing system...")
        build_dir = self._build_dir(pH)
        if build_dir.exists() and any(build_dir.iterdir()):
            self.log(f"[{self.runprefix} pH {pH}] Reusing existing build: {build_dir}")
            return str(build_dir)

        mol = Molecule(str(pdb_path))
        mol.filter('protein')
        mol_op = systemPrepare(mol, pH=pH)
        mol_seg = autoSegment(mol_op)
        mol_solv = solvate(mol_seg, pad=10)

        build_dir.parent.mkdir(parents=True, exist_ok=True)
        amber.build(mol_solv, outdir=str(build_dir))
        return str(build_dir)

    def equilibrate(self, build_dir: str, pH: float) -> str:
        self.log(f"[{self.runprefix} pH {pH}] Equilibration...")
        equil_dir = self._equil_dir(pH)
        if equil_dir.exists() and any(equil_dir.iterdir()):
            self.log(f"[{self.runprefix} pH {pH}] Reusing existing equil: {equil_dir}")
            return str(equil_dir)

        md = Equilibration()
        md.runtime = 20           # ns
        md.timeunits = 'ns'
        md.temperature = 300
        md.useconstantratio = False
        md.write(build_dir, str(equil_dir))

        q = LocalGPUQueue()
        q.submit(str(equil_dir))
        q.wait()
        return str(equil_dir)

    def production(self, equil_dir: str, pH: float) -> str:
        self.log(f"[{self.runprefix} pH {pH}] Production...")
        prod_dir = self._prod_dir(pH)
        # If prod already has outputs, we don't overwrite; ACEMD can resume if configured, but here we treat it as done.
        if prod_dir.exists() and any(prod_dir.iterdir()):
            self.log(f"[{self.runprefix} pH {pH}] Found existing prod dir: {prod_dir} (will continue/validate)")
        md = Production()
        md.runtime = 300          # ns
        md.timeunits = 'ns'
        md.temperature = 300
        # seed prod from equil outputs
        md.acemd.bincoordinates = 'output.coor'
        md.acemd.extendedsystem = 'output.xsc'
        md.write(equil_dir, str(prod_dir))

        q = LocalGPUQueue()
        q.submit(str(prod_dir))
        q.wait()
        return str(prod_dir)

    def run_pH(self, pdb_path: Path, pH: float):
        if self.is_done(pH):
            self.log(f"[{self.runprefix} pH {pH}] SKIP — already completed.")
            return

        start = time.time()
        try:
            self._root(pH).mkdir(parents=True, exist_ok=True)

            bdir = self.prepare_system(pdb_path, pH)
            edir = self.equilibrate(bdir, pH)
            pdir = self.production(edir, pH)

            # Validate prod outputs
            paths = self._prod_paths(pH)
            if any(good_file(paths[k]) for k in ["coor", "xsc", "log"]):
                self.mark_done(pH)
                self.log(f"[{self.runprefix} pH {pH}]  Marked as completed.")
            else:
                self.log(f"[{self.runprefix} pH {pH}]  Finished but outputs look incomplete.")
        except Exception as e:
            self.log(f"[{self.runprefix} pH {pH}]  ERROR: {repr(e)}")
            errfile = self._root(pH) / 'error.txt'
            errfile.parent.mkdir(parents=True, exist_ok=True)
            with open(errfile, 'a') as ef:
                ef.write(f"{datetime.now().isoformat()} - {repr(e)}\n")
        finally:
            dt = time.time() - start
            self.log(f"[{self.runprefix} pH {pH}] Elapsed: {dt:.1f}s")

# ---------- Main ----------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Rerun ACEMD only for incomplete hexokinase simulations.")
    parser.add_argument("--base", default=str(APO_DIR),
                        help="Base directory containing PDBs and pH subfolders")
    parser.add_argument("--log", default="simulationrerun_log.txt", help="Log file under base directory")
    parser.add_argument("--ph", nargs="+", type=float, default=[3.0, 5.0, 7.0],
                        help="pH values to run (default: 2.0 3.0 5.0 7.0)")
    parser.add_argument("--targets", nargs="+",
                        help="Explicit targets (e.g., 1qha 1v4s_D205N 3hm8_T546M 2nzt_D84Y). "
                             "If omitted, scans for all *.pdb in base and attempts reruns where incomplete.")
    parser.add_argument("--dry", action="store_true", help="Dry-run: print what would be executed and exit.")
    args = parser.parse_args()

    base_dir = Path(args.base).resolve()
    log_file = base_dir / args.log

    # Build target list
    if args.targets:
        wanted = []
        for t in args.targets:
            runprefix, pdbpath = norm_pdb_id_to_file(base_dir, t)
            if not pdbpath.exists():
                print(f"WARNING: PDB not found: {pdbpath} — skipping")
                continue
            wanted.append((runprefix, pdbpath))
    else:
        # auto-scan: every *.pdb in base_dir
        wanted = []
        for pdbpath in sorted(base_dir.glob("*.pdb")):
            wanted.append((pdbpath.stem, pdbpath))

    # DRY plan print
    if args.dry:
        print("Dry-run: planned (target × pH):")
        for runprefix, pdbpath in wanted:
            print(f"  - {runprefix}: file={pdbpath.name}  pH={', '.join(str(p) for p in args.ph)}")
        raise SystemExit(0)

    # Run without any lock
    for runprefix, pdbpath in wanted:
        sim = MolecularSimulation(runprefix, base_dir, log_file)
        sim.log(f"\n=== {runprefix} ({pdbpath.name}) ===")
        for pH in args.ph:
            sim.run_pH(pdbpath, pH)


