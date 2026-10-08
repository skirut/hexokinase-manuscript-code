"""Comparative binding-energy plotting for selected Hexokinase systems.

Purpose
-------
Reads selected ``binding_energy.dat`` files from the study data tree,
calculates total interaction energy as Elec + VdW, applies a 20-frame
moving average, and plots the selected systems together.

Inputs
------
- ``$HEXOKINASE_DATA_ROOT/G6P_structure/.../prod_extended_500ns/binding_energy.dat``
- ``$HEXOKINASE_DATA_ROOT/structures/1v4s_pH5.0/prod_extended_500ns/binding_energy.dat``
- The exact relative system list is defined in ``systems`` in the script.

Outputs
-------
- ``Combined_Total_Binding_Energy_500ns.jpg`` in the current working
  directory.

Run from
--------
Run from the repository root. The output file is written to the current
working directory rather than ``results/``.

Usage
-----
``python scripts/binding_energy/energy.py``

Dependencies
------------
Python, pandas, Matplotlib, and ``scripts.utils.paths``.
"""

import pandas as pd
import matplotlib.pyplot as plt
import os

from scripts.utils.paths import DATA_ROOT

def plot_combined_binding_energy():
    base_dir = DATA_ROOT
    
    systems = [
        "G6P_structure/2nzt_model_with_G6P_pH7.0",
        "G6P_structure/2nzt_model_with_G6P_pH5.0",
        "G6P_structure/2nzt_D84Y_G6P_pH7.0",
        "G6P_structure/2nzt_D84Y_G6P_pH5.0",
        "G6P_structure/3hm8_D867Y_G6P_pH5.0",
        "G6P_structure/3hm8_D867Y_G6P_pH3.0",
        "G6P_structure/1QHA_G87E_G6P_pH7.0",
        "structures/1v4s_pH5.0"
    ]
    
    plt.figure(figsize=(14, 7))
    found_any = False
    
    for system in systems:
        data_file = os.path.join(base_dir, system, "prod_extended_500ns", "binding_energy.dat")
        system_name = system.split('/')[-1]
        
        if not os.path.exists(data_file):
            print(f"Skipping {system_name}: File not found at {data_file}")
            continue
        
        print(f"Processing data for {system_name}...")
        df = pd.read_csv(data_file, sep=r'\s+', comment="#", engine='python')
        
        # 1. Filter out Frame 0 artifact (Pyruvate Kinase logic)
        if "Frame" in df.columns:
            df = df[df["Frame"] != 0]
        else:
            df = df.iloc[1:]  # Drop first row if no header column exists
            
        # 2. Calculate Elec + VdW interaction energy
        if "Elec" in df.columns and "VdW" in df.columns:
            total_energy = df["Elec"] + df["VdW"]
        elif "Total" in df.columns:
            total_energy = df["Total"]
        else:
            total_energy = df.iloc[:, 3]
            
        time_val = df["Time"] if "Time" in df.columns else df.iloc[:, 0]
        
        # 3. Calculate 20-frame moving average for smoothing
        energy_ma = total_energy.rolling(window=20, min_periods=1).mean()
        
        # Plot raw trajectory (faint background) and smoothed trendline
        plt.plot(time_val, total_energy, color="gray", alpha=0.15, linewidth=0.5)
        plt.plot(time_val, energy_ma, label=system_name, linewidth=1.8)
        found_any = True
        
    if not found_any:
        print("No valid binding_energy.dat files were found.")
        return

    # Chart formatting
    plt.title('Binding Energy Dynamics (500 ns)', fontsize=14, fontweight='bold')
    plt.xlabel('Time / Frame', fontsize=12)
    plt.ylabel('Binding Energy (kcal/mol)', fontsize=12)
    plt.axhline(0, color='gray', linestyle='--', linewidth=1)
    plt.grid(True, linestyle=':', alpha=0.7)
    
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=9.5, frameon=True)
    plt.tight_layout()
    
    output_filename = "Combined_Total_Binding_Energy_500ns.jpg"
    plt.savefig(output_filename, dpi=300, format='jpg', bbox_inches='tight')
    print(f"\nSaved clean comparative plot to: {output_filename}")
    plt.close()

if __name__ == '__main__':
    plot_combined_binding_energy()

