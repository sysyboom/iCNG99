import pandas as pd
import subprocess

input_excel_path = '/mnt/NFS/fengch/beforecxcalc.xlsx'
output_excel_path = '/mnt/NFS/fengch/aftercxcalc.xlsx'

df = pd.read_excel(input_excel_path)

df['MajorMicrospecies'] = ''

for index, row in df.iterrows():
    smiles = row['InChI']  
    if smiles:
        command = f'cxcalc majormicrospecies -H 7.2 "{smiles}"'

        try:
            result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            major_microspecies = result.stdout.strip()
            
            df.at[index, 'MajorMicrospecies'] = major_microspecies
        except subprocess.CalledProcessError as e:
            error_message = e.stderr.strip()
            print(f"Error executing command for SMILES {smiles}: {error_message}")

df.to_excel(output_excel_path, index=False)
