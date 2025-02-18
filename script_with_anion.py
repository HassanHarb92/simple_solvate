import os
import subprocess
import numpy as np
import argparse
from math import pi

# Approximate solvent molecular volumes (Å³)
SOLVENT_VOLUMES = {
    "Water": 30.0,
    "Methanol": 40.0,
    "Ethanol": 60.0,
    "Acetone": 70.0,
    "DMSO": 90.0,
    "Hexane": 130.0,
}

def parse_arguments():
    parser = argparse.ArgumentParser(description="Solvated system generation and xTB optimization.")
    parser.add_argument("--structure", required=True, help="Path to the input XYZ structure.")
    parser.add_argument("--charge", type=int, required=True, help="Charge of the system.")
    parser.add_argument("--solvent", required=True, help="Name of the solvent.")
    parser.add_argument("--solvent_SMILES", required=True, help="SMILES string of the solvent.")
    parser.add_argument("--counter_ion", default="NONE", help="Name of the counter ion.")
    parser.add_argument("--counter_ion_SMILES", help="SMILES string of the counter ion.")
    parser.add_argument("--counter_ion_count", type=int, default=0, help="Number of counter ions.")
    return parser.parse_args()

def read_xyz(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    return [(line.split()[0], list(map(float, line.split()[1:]))) for line in lines[2:] if len(line.split()) == 4]

def generate_xyz_from_smiles(smiles, name, script_path="/Users/hharb/Desktop/Codes/simple_solvate/SMILES2XYZ.py"):
    output_file = f"{name}.xyz"
    if os.path.exists(output_file):
        return output_file
    command = ["python", script_path, smiles, name]
    print(f"[INFO] Generating XYZ for {name} ({smiles})...")
    subprocess.run(command, check=True)
    if not os.path.exists(output_file):
        raise RuntimeError(f"[ERROR] Failed to generate XYZ file for {name}")
    return output_file

def calculate_center_of_mass(xyz_data):
    coords = np.array([coord for _, coord in xyz_data])
    return np.mean(coords, axis=0)

def calculate_max_distance(center, xyz_data):
    return max(np.linalg.norm(np.array(coord) - center) for _, coord in xyz_data)

def place_counter_ions(main_coords, counter_ion_xyz, count):
    center = calculate_center_of_mass(main_coords)
    max_distance = calculate_max_distance(center, main_coords) + 3.0
    counter_ions = []
    for i in range(count):
        sign = 1 if i % 2 == 0 else -1
        ion_coords = read_xyz(counter_ion_xyz)
        ion_com = calculate_center_of_mass(ion_coords)
        translation = np.array([center[0] + sign * max_distance, center[1], center[2]]) - ion_com
        counter_ions.extend([(atom, tuple(np.array(coord) + translation)) for atom, coord in ion_coords])
    return counter_ions

def distribute_solvents_on_sphere(center, radius, n_solvents, solvent_xyz):
    solvents = []
    golden_ratio = (1 + 5**0.5) / 2
    solvent_coords = read_xyz(solvent_xyz)
    solvent_com = calculate_center_of_mass(solvent_coords)
    for i in range(n_solvents):
        theta = 2 * pi * i / golden_ratio
        phi = np.arccos(1 - 2 * (i + 0.5) / n_solvents)
        new_pos = np.array([center[0] + radius * np.sin(phi) * np.cos(theta),
                            center[1] + radius * np.sin(phi) * np.sin(theta),
                            center[2] + radius * np.cos(phi)])
        translation = new_pos - solvent_com
        solvents.extend([(atom, tuple(np.array(coord) + translation)) for atom, coord in solvent_coords])
    return solvents

def write_combined_xyz(filename, all_coords):
    with open(filename, 'w') as f:
        f.write(f"{len(all_coords)}\n")
        f.write("Solvated system with counter ions\n")
        for atom, (x, y, z) in all_coords:
            f.write(f"{atom}  {x:.6f}  {y:.6f}  {z:.6f}\n")

def create_frozen_input(main_xyz_file):
    frozen_count = sum(1 for _ in open(main_xyz_file)) - 2  # Exclude header
    with open("freeze.inp", "w") as f:
        f.write("$fix\n  atoms: " + " ".join(map(str, range(1, frozen_count + 1))) + "\n$end\n")

def run_xtb_optimization(xyz_file, charge, solvent):
    if not os.path.exists(xyz_file):
        raise RuntimeError(f"[ERROR] No solvated system found: {xyz_file}")
    cmd = ["xtb", xyz_file, "--opt", "tight", "--gfn2", "--alpb", solvent, "--input", "freeze.inp", "--chrg", str(charge)]
#    cmd = ["xtb", xyz_file, "--opt", "tight", "--gfn2", "--cbonds", "--alpb", solvent, "--input", "freeze.inp", "--chrg", str(charge)]
    print(f"[INFO] Running xTB with solvent '{solvent}'...")
    subprocess.run(cmd, check=True)
    os.rename("xtbopt.xyz", "solvated_optimized.xyz")
    print("[INFO] Optimization complete! Output saved as 'solvated_optimized.xyz'.")

def main():
    args = parse_arguments()
    main_coords = read_xyz(args.structure)
    all_coords = main_coords
    if args.counter_ion.upper() != "NONE":
        counter_ion_xyz = generate_xyz_from_smiles(args.counter_ion_SMILES, args.counter_ion)
        all_coords += place_counter_ions(main_coords, counter_ion_xyz, args.counter_ion_count)
    if args.solvent.upper() != "NONE":
        solvent_xyz = generate_xyz_from_smiles(args.solvent_SMILES, args.solvent)
        radius = (3 * SOLVENT_VOLUMES.get(args.solvent, 100) / (4 * pi))**(1/3) + 5.0
        solvents = distribute_solvents_on_sphere(calculate_center_of_mass(main_coords), radius, 20, solvent_xyz)
        all_coords += solvents
    write_combined_xyz("solvated_system.xyz", all_coords)
    create_frozen_input(args.structure)
    run_xtb_optimization("solvated_system.xyz", args.charge, args.solvent)

if __name__ == "__main__":
    main()

