import os
import subprocess
import numpy as np
from math import pi
import sys

# Approximate solvent molecular volumes (Å³)
SOLVENT_VOLUMES = {
    "Water": 30.0,
    "Methanol": 40.0,
    "Ethanol": 60.0,
    "Acetone": 70.0,
    "DMSO": 90.0,
    "Hexane": 130.0,
}

def read_xyz(filename):
    """Reads an XYZ file and returns a list of (atom_label, (x, y, z))."""
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f]

    if len(lines) < 3:
        raise ValueError(f"[ERROR] '{filename}' does not contain enough lines to be a valid XYZ file.")

    declared_atom_count = int(lines[0])
    actual_atom_count = sum(1 for line in lines[2:] if len(line.split()) == 4)

    xyz_data = []
    for line in lines[2:2 + actual_atom_count]:
        parts = line.split()
        atom_label = parts[0]
        x, y, z = map(float, parts[1:4])
        xyz_data.append((atom_label, (x, y, z)))
    
    return xyz_data

def count_atoms(filename):
    """Counts the number of atoms in an XYZ file."""
    with open(filename, 'r') as f:
        return int(f.readline().strip())

def oldparse_input_file(filename):
    """Parses the input file and extracts required parameters."""
    params = {}
    with open(filename, 'r') as f:
        for line in f:
            if line.strip():
                key, value = line.strip().split(maxsplit=1)
                params[key.strip()] = value.strip()

    # Validate required fields
    required_fields = ["--structure", "--charge"]
    for field in required_fields:
        if field not in params:
            raise ValueError(f"[ERROR] Missing required field: {field}")

    return params

def parse_input_file(filename):
    """Parses the input file and extracts required parameters."""
    params = {}

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):  # Ignore empty lines or comments
                continue

            parts = line.split(maxsplit=1)
            if len(parts) == 2:  # Proper key-value pair
                key, value = parts
                params[key.strip()] = value.strip()
            else:
                print(f"[WARNING] Skipping malformed line: {line}")

    # Validate required fields
    required_fields = ["--structure", "--charge"]
    for field in required_fields:
        if field not in params:
            raise ValueError(f"[ERROR] Missing required field: {field}")

    return params



def generate_xyz_from_smiles(smiles, name, script_path="/Users/hharb/Desktop/Codes/simple_solvate/SMILES2XYZ.py"):
    """Generates an XYZ file from a SMILES string using an external script."""
    output_file = f"{name}.xyz"
    command = ["python", script_path, smiles, name]

    print(f"[INFO] Generating XYZ for {name} ({smiles})...")
    subprocess.run(command, capture_output=True, text=True)

    if not os.path.exists(output_file):
        raise RuntimeError(f"[ERROR] Failed to generate XYZ file for {name}")

    return output_file

def calculate_center_of_mass(xyz_data):
    """Calculates the center of mass of a molecule."""
    coords = np.array([coord for _, coord in xyz_data])
    return np.mean(coords, axis=0)

def calculate_max_distance(center, xyz_data):
    """Finds the maximum distance from the center of mass to any atom."""
    return max(np.linalg.norm(np.array(coord) - center) for _, coord in xyz_data)

def place_counter_ions(main_coords, counter_ion_xyz, count):
    """Places counter ions at opposite ends of the main molecule."""
    center = calculate_center_of_mass(main_coords)
    max_distance = calculate_max_distance(center, main_coords) + 3.0

    counter_ions = []
    for i in range(count):
        sign = 1 if i % 2 == 0 else -1
        ion_coords = read_xyz(counter_ion_xyz)
        ion_com = calculate_center_of_mass(ion_coords)

        translation = np.array([center[0] + sign * max_distance, center[1], center[2]]) - ion_com
        translated_ions = [(atom, tuple(np.array(coord) + translation)) for atom, coord in ion_coords]
        counter_ions.extend(translated_ions)

    return counter_ions

def distribute_solvents_on_sphere(center, radius, n_solvents, solvent_xyz):
    """Distributes full solvent molecules evenly on the surface of a sphere."""
    solvents = []
    golden_ratio = (1 + 5**0.5) / 2
    solvent_coords = read_xyz(solvent_xyz)
    solvent_com = calculate_center_of_mass(solvent_coords)

    for i in range(n_solvents):
        theta = 2 * pi * i / golden_ratio
        phi = np.arccos(1 - 2 * (i + 0.5) / n_solvents)

        new_x = center[0] + radius * np.sin(phi) * np.cos(theta)
        new_y = center[1] + radius * np.sin(phi) * np.sin(theta)
        new_z = center[2] + radius * np.cos(phi)

        translation = np.array([new_x, new_y, new_z]) - solvent_com
        translated_solvent = [(atom, tuple(np.array(coord) + translation)) for atom, coord in solvent_coords]
        solvents.extend(translated_solvent)

    return solvents

def write_combined_xyz(filename, all_coords):
    """Writes an XYZ file with all molecules (main, counter ions, solvent)."""
    with open(filename, 'w') as f:
        f.write(f"{len(all_coords)}\n")
        f.write("Solvated system with counter ions\n")
        for atom_label, (x, y, z) in all_coords:
            f.write(f"{atom_label}  {x:.6f}  {y:.6f}  {z:.6f}\n")

def create_frozen_input(main_xyz_file):
    """Creates an xTB input file that correctly freezes ONLY the main molecule's atoms."""
    frozen_count = count_atoms(main_xyz_file)
    
    with open("freeze.inp", "w") as f:
        f.write("$fix\n")
        f.write("  atoms: ")
        f.write(" ".join(str(i) for i in range(1, frozen_count + 1)))
        f.write("\n$end\n")
    
    print(f"[INFO] Created freeze.inp (freezing atoms 1 to {frozen_count})")

def run_xtb_optimization(xyz_file, charge, solvent, output_name):
    """Runs xTB geometry optimization while freezing the main structure."""
    if not os.path.exists(xyz_file):
        raise RuntimeError(f"[ERROR] No solvated system found: {xyz_file}")

    cmd = ["xtb", xyz_file, "--opt", "tight", "--gfn2", "--alpb", solvent, "--input", "freeze.inp", "--chrg", str(charge)]
    print(f"[INFO] Running xTB with solvent '{solvent}'...")
    subprocess.run(cmd, check=True)
    os.rename("xtbopt.xyz", output_name)
    print(f"[INFO] Optimization complete! Output saved as '{output_name}'.")

def main():
    input_file = sys.argv[1]
    params = parse_input_file(input_file)

    main_coords = read_xyz(params["--structure"])
    all_coords = main_coords

    if params["--counter_ion"].upper() != "NONE":
        counter_ion_xyz = generate_xyz_from_smiles(params["--counter_ion_SMILES"], params["--counter_ion"])
        all_coords += place_counter_ions(main_coords, counter_ion_xyz, int(params["--counter_ion_count"]))

        write_combined_xyz("counter_ion_added.xyz", all_coords)
        create_frozen_input(params["--structure"])
        run_xtb_optimization("counter_ion_added.xyz", params["--charge"], params["--solvent"], "counter_ion_optimized.xyz")

    if params["--solvent"].upper() != "NONE":
        solvent_xyz = generate_xyz_from_smiles(params["--solvent_SMILES"], params["--solvent"])
        solvents = distribute_solvents_on_sphere(calculate_center_of_mass(main_coords), 10.0, 20, solvent_xyz)
        all_coords += solvents

        write_combined_xyz("solvated_system.xyz", all_coords)
        run_xtb_optimization("solvated_system.xyz", params["--charge"], params["--solvent"], "solvated_optimized.xyz")

if __name__ == "__main__":
    main()

