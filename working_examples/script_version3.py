import os
import subprocess
import numpy as np
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

def read_xyz(filename):
    """Reads an XYZ file and returns a list of (atom_label, (x, y, z))."""
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f ] #if line.strip()]

    if len(lines) < 3:
        raise ValueError(f"[ERROR] '{filename}' does not contain enough lines to be a valid XYZ file.")

    declared_atom_count = int(lines[0])
    actual_atom_count = sum(1 for line in lines[2:] if len(line.split()) == 4)

    if declared_atom_count != actual_atom_count:
        print(f"[WARNING] Atom count mismatch in '{filename}': Declared={declared_atom_count}, Found={actual_atom_count}. Using {actual_atom_count}.")

    xyz_data = []
    for line in lines[2:2 + actual_atom_count]:
        parts = line.split()
        if len(parts) != 4:
            raise ValueError(f"[ERROR] Invalid line format in '{filename}': '{line}'. Ensure it follows 'Element X Y Z'.")
        
        atom_label = parts[0]
        x, y, z = map(float, parts[1:4])
        xyz_data.append((atom_label, (x, y, z)))
    
    return xyz_data

def generate_solvent_xyz(smiles, solvent_name, script_path="/Users/hharb/Desktop/Codes/simple_solvate/SMILES2XYZ.py"):
    """
    Uses SMILES2XYZ.py to generate an XYZ file for the solvent molecule.
    """
    output_file = f"{solvent_name}.xyz"
    command = ["python", script_path, smiles, solvent_name]

    print(f"[INFO] Generating XYZ for solvent: {solvent_name} ({smiles})...")

    result = subprocess.run(command, capture_output=True, text=True)

    if not os.path.exists(output_file):
        raise RuntimeError(f"[ERROR] Failed to generate XYZ file for {solvent_name}")

    return output_file

def calculate_center_of_mass(xyz_data):
    """Calculates the center of mass of a molecule."""
    coords = np.array([coord for _, coord in xyz_data])
    return np.mean(coords, axis=0)

def calculate_max_distance(center, xyz_data):
    """Finds the maximum distance from the center of mass to any atom."""
    return max(np.linalg.norm(np.array(coord) - center) for _, coord in xyz_data)

def calculate_sphere_volume(radius):
    """Calculates the volume of a sphere given its radius."""
    return (4/3) * pi * (radius ** 3)

def distribute_solvents_on_sphere(center, radius, n_solvents, solvent_xyz):
    """
    Distributes full solvent molecules evenly on the surface of a sphere.
    Moves the solvent molecule to each position.
    """
    solvents = []
    golden_ratio = (1 + 5**0.5) / 2  # Approximate spherical distribution
    solvent_coords = read_xyz(solvent_xyz)  # Read solvent molecule structure
    solvent_com = calculate_center_of_mass(solvent_coords)  # Get solvent center

    for i in range(n_solvents):
        theta = 2 * pi * i / golden_ratio  # Azimuthal angle
        phi = np.arccos(1 - 2 * (i + 0.5) / n_solvents)  # Polar angle

        # Determine new solvent center position on the sphere surface
        new_x = center[0] + radius * np.sin(phi) * np.cos(theta)
        new_y = center[1] + radius * np.sin(phi) * np.sin(theta)
        new_z = center[2] + radius * np.cos(phi)

        # Translation vector to move solvent from its COM to the new position
        translation = np.array([new_x, new_y, new_z]) - solvent_com

        # Move each solvent atom to the new position
        translated_solvent = [(atom, tuple(np.array(coord) + translation)) for atom, coord in solvent_coords]
        solvents.extend(translated_solvent)

    return solvents

def write_combined_xyz(filename, main_coords, solvent_coords):
    """Writes the combined XYZ file with main molecule and solvent molecules."""
    combined_coords = main_coords + solvent_coords
    with open(filename, 'w') as f:
        f.write(f"{len(combined_coords)}\n")
        f.write("Molecule with solvated sphere\n")
        for atom_label, (x, y, z) in combined_coords:
            f.write(f"{atom_label}  {x:.6f}  {y:.6f}  {z:.6f}\n")

def main():
    main_xyz = input("Enter the main molecule (frozen) XYZ file: ")
    main_coords = read_xyz(main_xyz)

    # Compute solvation sphere
    com = calculate_center_of_mass(main_coords)
    max_distance = calculate_max_distance(com, main_coords)
    solvation_radius = max_distance + 2.5
    sphere_volume = calculate_sphere_volume(solvation_radius)

    print(f"[INFO] Solvation sphere radius: {solvation_radius:.2f} Å")
    print(f"[INFO] Solvation sphere volume: {sphere_volume:.2f} Å³")

    # Get solvent choice
    solvent_name = input("Enter the solvent name (e.g., 'Water'): ").strip()
    solvent_smiles = input("Enter the SMILES string of the solvent: ").strip()

    # Generate solvent XYZ file
    solvent_xyz = generate_solvent_xyz(solvent_smiles, solvent_name)

    # Estimate solvent molecule volume
    solvent_volume = SOLVENT_VOLUMES.get(solvent_name, 50.0)  # Default if not found
    print(f"[INFO] Estimated solvent volume: {solvent_volume:.2f} Å³")

    # Determine number of solvent molecules
    n_solvents = min(int(sphere_volume / solvent_volume), 20)
    print(f"[INFO] Placing approximately {n_solvents} solvent molecules.")

    # Distribute full solvent molecules on the sphere
    solvent_coords = distribute_solvents_on_sphere(com, solvation_radius, n_solvents, solvent_xyz)

    # Write the final solvated system
    output_file = "solvated_system.xyz"
    write_combined_xyz(output_file, main_coords, solvent_coords)

    print(f"[INFO] Solvated system written to {output_file}")

if __name__ == "__main__":
    main()

