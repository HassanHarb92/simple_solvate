import os
import subprocess
import numpy as np
from math import pi

def count_atoms(filename):
    """Counts the number of atoms in an XYZ file."""
    with open(filename, 'r') as f:
        return int(f.readline().strip())

def create_frozen_input(frozen_count):
    """
    Creates an xTB input file that correctly freezes the first 'frozen_count' atoms.
    """
    with open("freeze.inp", "w") as f:
        f.write("$fix\n")
        f.write("  atoms: ")
        f.write(" ".join(str(i) for i in range(1, frozen_count + 1)))  # xTB uses 1-based index
        f.write("\n$end\n")
    print(f"[INFO] Created freeze.inp (freezing atoms 1 to {frozen_count})")

def run_xtb_opt_frozen(xyz_file, charge=0, solvent="water"):
    """
    Runs xTB geometry optimization with frozen main molecule and user-defined solvent.
    """
    cmd = [
        "xtb",
        xyz_file,
        "--opt", "tight",
        "--gfn2",
        "--alpb", solvent,
        "--input", "freeze.inp",
        "--chrg", str(charge)
    ]
    print(f"[INFO] Running xTB with solvent '{solvent}'...")
    completed_process = subprocess.run(cmd, capture_output=True, text=True)
    
    if completed_process.returncode != 0:
        print("[ERROR] xTB optimization failed:")
        print(completed_process.stderr)
    else:
        print("[INFO] xTB optimization completed successfully!")
        if os.path.exists("xtbopt.xyz"):
            os.rename("xtbopt.xyz", "solvated_optimized.xyz")
            print("[INFO] Renamed 'xtbopt.xyz' to 'solvated_optimized.xyz'")

def main():
    main_xyz = "Citrate.xyz"  # The frozen molecule
    solvated_xyz = "solvated_system.xyz"  # The system with solvents

    # Count the number of atoms in the main molecule to freeze them
    frozen_count = count_atoms(main_xyz)
    create_frozen_input(frozen_count)

    # Ask for system charge and solvent
    total_charge = input("Enter total net charge for the solvated system [0]: ").strip()
    if not total_charge:
        total_charge = "0"
    total_charge = int(total_charge)

    solvent_model = input("Enter solvent model for xTB (default = 'water'): ").strip()
    if not solvent_model:
        solvent_model = "water"

    # Run xTB with the frozen atoms
    run_xtb_opt_frozen(solvated_xyz, charge=total_charge, solvent=solvent_model)

    print("[INFO] Done. Check 'solvated_optimized.xyz' for the final structure.")

if __name__ == "__main__":
    main()

