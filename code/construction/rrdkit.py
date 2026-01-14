import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def process_excel(input_file, output_file):
    df = pd.read_excel(input_file)
    
    if 'MajorMicrospecies' not in df.columns:
        print("The 'MajorMicrospecies' column is missing in the input file.")
        return

    df['Chemical_Formula'] = ''
    df['Net_Charge'] = ''

    for index, row in df.iterrows():
        inchi = row['MajorMicrospecies']

        if pd.notnull(inchi) and isinstance(inchi, str):
            chemical_formula, net_charge = get_molecular_details_from_inchi(inchi)
        else:
            chemical_formula, net_charge = 'Invalid InChI', 0
        df.at[index, 'Chemical_Formula'] = chemical_formula
        df.at[index, 'Net_Charge'] = net_charge

    df.to_excel(output_file, index=False, engine='openpyxl')

def get_molecular_details_from_inchi(inchi):
    mol = Chem.MolFromInchi(inchi)
    if mol is None:
        return 'Invalid InChI', 0
    chemical_formula = rdMolDescriptors.CalcMolFormula(mol)
    net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    return chemical_formula, net_charge

input_excel_path = '/mnt/NFS/fengch/aftercxcalcadd.xlsx'  
output_excel_path = '/mnt/NFS/fengch/afterrdkitadd.xlsx' 

process_excel(input_excel_path, output_excel_path)

