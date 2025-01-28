import os
import subprocess

# Define the solvent dictionary
solvent_dict = {
    "Water": "O",
    "Methanol": "CO",
    "Ethanol": "CCO",
    "Acetone": "CC(=O)C",
    "Chloroform": "C(Cl)(Cl)Cl",
    "Benzene": "c1ccccc1",
    "Toluene": "CC1=CC=CC=C1",
    "Acetonitrile": "CC#N",
    "Dimethylformamide (DMF)": "CN(C)C=O",
    "Dimethyl sulfoxide (DMSO)": "CS(=O)C",
    "Diethyl ether": "CCOCC",
    "Ethyl acetate": "CC(=O)OCC",
    "Dichloromethane (DCM)": "C(Cl)Cl",
    "Hexane": "CCCCCC",
    "Heptane": "CCCCCCC",
    "Tetrahydrofuran (THF)": "C1CCOC1",
    "1,4-Dioxane": "C1COCCO1",
    "Pyridine": "C1=CC=NC=C1",
    "Formic acid": "OC=O",
    "Acetic acid": "CC(=O)O",
    "Propanol": "CCCO",
    "Butanol": "CCCCO",
    "Isopropanol": "CC(O)C",
    "Anisole": "COc1ccccc1"
}

# Define path to the script
script_path = "/Users/hharb/Desktop/Codes/simple_solvate/SMILES2XYZ.py"

# Store invalid SMILES
invalid_smiles = []

# Run test for each solvent
for solvent, smiles in solvent_dict.items():
    output_file = f"{solvent}.xyz"
    
    # Construct command
    command = ["python", script_path, smiles, solvent]

    print(f"[INFO] Testing: {solvent} ({smiles})...")

    try:
        # Run subprocess
        result = subprocess.run(command, capture_output=True, text=True)

        # Check if .xyz file was created
        if not os.path.exists(output_file):
            print(f"[ERROR] Conversion failed for: {solvent} ({smiles})")
            invalid_smiles.append((solvent, smiles))
        else:
            print(f"[SUCCESS] {solvent} converted successfully to {output_file}")

    except Exception as e:
        print(f"[EXCEPTION] Error processing {solvent}: {e}")
        invalid_smiles.append((solvent, smiles))

# Print summary
print("\n========== Test Summary ==========")
if invalid_smiles:
    print("[WARNING] The following SMILES failed conversion:")
    for solvent, smiles in invalid_smiles:
        print(f" - {solvent}: {smiles}")
else:
    print("[SUCCESS] All solvents converted successfully!")


