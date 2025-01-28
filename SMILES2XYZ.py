import sys

def Smiles2XYZ(molecule_name,molecule_smiles):
    from rdkit import Chem
    from rdkit.Chem import Draw
    from rdkit.Chem.Draw import IPythonConsole
    from rdkit.Chem import Descriptors
    from rdkit.Chem import AllChem
    from rdkit import DataStructs
    import numpy as np

    functional = 'uwb97xd'
    Basis_set = '6-31G(2df,p)'
    route_line = '#p ' + functional + ' ' + Basis_set + ' opt=tight freq '
    molecule_label = molecule_smiles
    molecule_smiles = Chem.MolFromSmiles(molecule_smiles)
    molecule_smiles_1 = molecule_smiles
    molecule_smiles = Chem.AddHs(molecule_smiles)
    nAtoms = molecule_smiles.GetNumAtoms()
    print (nAtoms)

    AllChem.Compute2DCoords(molecule_smiles)
    AllChem.EmbedMolecule(molecule_smiles)
    AllChem.MMFFOptimizeMolecule(molecule_smiles)
    charge = Chem.GetFormalCharge(molecule_smiles)
    multiplicity = '1'
   # xyz_string = "\n{} {}\n".format(charge, multiplicity)
    xyz_string = ''
    #print ('Charge = ', charge)
    for atom in molecule_smiles.GetAtoms():
        pos = molecule_smiles.GetConformer().GetAtomPosition(atom.GetIdx())
        xyz_string += "{} {} {} {}\n".format(atom.GetSymbol(),pos.x, pos.y, pos.z)

    #print (xyz_string)

    inputfile = molecule_name+'.xyz'
    checkpoint = molecule_name+'.chk'

    with open(inputfile, 'w') as inpt:
 #       inpt.write("%chk=")
 #       inpt.write(checkpoint)
 #       inpt.write('\n')
 #       inpt.write(route_line)
 #       inpt.write('\n\n')
 #       inpt.write('Optimization for ')
        inpt.write(str(nAtoms))
        inpt.write('\n\n')
        inpt.write(xyz_string)


    molecule_smiles_1

test = 'c1ccccc1'
test_name = 'benzene'
test = sys.argv[1]
test_name = sys.argv[2]

Smiles2XYZ(test_name,test)

