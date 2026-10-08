#!/usr/bin/env python

"""
Extend selected Hexokinase production simulations from 300 ns to 500 ns.

Purpose
-------
Prepares and sequentially runs 500 ns ACEMD extensions for selected
Hexokinase molecular dynamics systems. The script copies an existing
300 ns production directory, modifies the ACEMD input to use the
500 ns run length, and launches ACEMD from the extended directory.

The script also maintains a log so that systems already completed
successfully in a previous run can be skipped.

Inputs
------
- HEXOKINASE_DATA_ROOT environment variable.
- Existing Hexokinase simulation directories under the data root.
- A ``prod`` directory containing:
    - ``restart.chk``
    - ``input``

The original production configuration is expected to correspond to
75,000,000 steps (300 ns at a 4 fs timestep).

Outputs
-------
For each selected system:

- ``prod_extended_500ns/`` directory containing the copied production
  setup and modified ACEMD input.
- ``extended.log`` containing ACEMD output.
- ``hk_extension_500ns_log.txt`` containing script-level progress and
  completion messages.

Run from
--------
Repository root:

    /path/to/hexokinase_simulation

Set the external data root before running:

    export HEXOKINASE_DATA_ROOT=/path/to/hexokinase

Usage
-----
Run:

    python scripts/simulation/extended_simulation.py

The script processes the selected systems sequentially.

Dependencies
------------
- Python standard library
- ACEMD
- Existing Hexokinase simulation directories and restart checkpoints
"""
import os
import shutil
import subprocess
import time
import sys
from datetime import datetime

class HexokinaseSimulationExtender:
    def __init__(self, base_dir):
        self.base_dir = base_dir
        self.log_file = os.path.join(base_dir, 'hk_extension_500ns_log.txt')
        self.extension_steps = 125000000  # 500 ns with 4 fs timestep
        self.original_steps = 75000000    # 300 ns with 4 fs timestep
        
    def log(self, message):
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        with open(self.log_file, 'a') as f:
            f.write(f"{timestamp}: {message}\n")
        print(f"{timestamp}: {message}")

    def extend_system(self, system_path):
        self.log(f"\n{'='*60}")
        self.log(f"Processing: {system_path}")
        
        # Original prod directory
        orig_dir = os.path.join(self.base_dir, system_path, 'prod')
        if not os.path.exists(orig_dir):
            self.log(f"ERROR: Original prod directory not found at {orig_dir}")
            return False
        
        # Check if checkpoint exists
        checkpoint = os.path.join(orig_dir, 'restart.chk')
        if not os.path.exists(checkpoint):
            self.log(f"ERROR: Checkpoint file not found at {checkpoint}")
            return False
        
        # Extended directory structure
        extended_dir = os.path.join(self.base_dir, system_path, 'prod_extended_500ns')
        
        # Remove if exists to prevent appending to old test files
        if os.path.exists(extended_dir):
            self.log(f"Removing existing directory: {extended_dir}")
            shutil.rmtree(extended_dir)
        
        # Copy prod directory
        self.log(f"Copying prod directory to: {extended_dir}")
        shutil.copytree(orig_dir, extended_dir)
        
        # Modify input file
        input_file = os.path.join(extended_dir, 'input')
        if not os.path.exists(input_file):
            self.log(f"ERROR: Input file not found at {input_file}")
            return False
        
        # Read the input file
        with open(input_file, 'r') as f:
            content = f.read()
        
        self.log(f"Original run value: {self.original_steps}")
        self.log(f"New run value: {self.extension_steps} (500 ns)")
        
        # Replace run value
        content = content.replace(f'run                  \t{self.original_steps}', 
                                   f'run                  \t{self.extension_steps}')
        
        # Write the modified input file
        with open(input_file, 'w') as f:
            f.write(content)
        
        self.log(f"Modified input file successfully")
        
        # Verify the change
        with open(input_file, 'r') as f:
            for line in f:
                if 'run' in line:
                    self.log(f"Verified: {line.strip()}")
        
        return extended_dir

    def run_simulation(self, extended_dir, system_path):
        """Run ACEMD and wait for completion"""
        self.log(f"Starting ACEMD for {system_path}")
        
        original_dir = os.getcwd()
        os.chdir(extended_dir)
        
        cmd = "acemd input > extended.log 2>&1"
        self.log(f"Running: {cmd}")
        
        try:
            process = subprocess.Popen(cmd, shell=True)
            pid = process.pid
            
            self.log(f"ACEMD started with PID: {pid}")
            self.log(f"Monitor progress: tail -f {extended_dir}/extended.log")
            
            process.wait()
            
            if process.returncode == 0:
                self.log(f"Simulation completed successfully for {system_path}")
                os.chdir(original_dir)
                return True
            else:
                self.log(f"Simulation failed for {system_path} (return code: {process.returncode})")
                os.chdir(original_dir)
                return False
                
        except Exception as e:
            self.log(f" Exception occurred for {system_path}: {e}")
            os.chdir(original_dir)
            return False

    def get_target_systems(self):
        """Returns the strategically selected Hexokinase systems."""
        return [
            "G6P_structure/2nzt_model_with_G6P_pH7.0",
            "G6P_structure/2nzt_model_with_G6P_pH5.0",
            "G6P_structure/2nzt_D84Y_G6P_pH7.0",
            "G6P_structure/2nzt_D84Y_G6P_pH5.0",
            "G6P_structure/3hm8_D867Y_G6P_pH5.0",
            "G6P_structure/3hm8_D867Y_G6P_pH3.0",
            "G6P_structure/1QHA_G87E_G6P_pH7.0",
            "structures/1v4s_pH5.0"
        ]

    def run_all_sequential(self):
        systems = self.get_target_systems()
        total_start = time.time()
        
        # Parse the log file to identify systems that completed in previous runs
        completed_systems = set()
        if os.path.exists(self.log_file):
            with open(self.log_file, 'r') as f:
                for line in f:
                    if "Simulation completed successfully for" in line:
                        system_name = line.strip().split("for ")[-1]
                        completed_systems.add(system_name)
        
        for i, system in enumerate(systems, 1):
            self.log(f"\n{'='*60}")
            self.log(f"Processing system {i}/{len(systems)}: {system}")
            
            # Skip system if marked as successfully finished in the log
            if system in completed_systems:
                self.log(f"SKIPPING: {system} (Already completed in a previous run)")
                continue
            
            extended_dir = self.extend_system(system)
            if not extended_dir:
                self.log(f" Failed to prepare {system}, skipping...")
                continue
            
            success = self.run_simulation(extended_dir, system)
            
            if success:
                self.log(f"Completed {i}/{len(systems)}: {system}")
            else:
                self.log(f" Failed {i}/{len(systems)}: {system}")
            
            elapsed = time.time() - total_start
            avg_time = elapsed / i
            remaining = (len(systems) - i) * avg_time
            self.log(f"Elapsed: {elapsed/60:.1f} min, Est. remaining: {remaining/60:.1f} min")
        
        total_time = time.time() - total_start
        self.log(f"\n{'='*60}")
        self.log(f"ALL EXTENSIONS COMPLETED in {total_time/60:.2f} minutes ({total_time/3600:.2f} hours)")

def main():
    base_dir = '/media/simulations/hexokinase'
    
    if not os.path.exists(base_dir):
        print(f"ERROR: Directory not found: {base_dir}")
        sys.exit(1)
    
    extender = HexokinaseSimulationExtender(base_dir)
    
    print(f"\n{'='*60}")
    print("ACEMD Hexokinase Simulation Extender")
    print(f"Base directory: {base_dir}")
    print(f"Extending to 500 ns")
    print(f"{'='*60}\n")
    
    extender.run_all_sequential()

if __name__ == '__main__':
    main()
