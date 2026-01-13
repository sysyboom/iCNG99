import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

# 读取Excel文件并处理InChI字符串
def process_excel(input_file, output_file):
    # 读取Excel文件
    df = pd.read_excel(input_file)

    # 确保'MajorMicrospecies'列存在
    if 'MajorMicrospecies' not in df.columns:
        print("The 'MajorMicrospecies' column is missing in the input file.")
        return

    # 初始化新列
    df['Chemical_Formula'] = ''
    df['Net_Charge'] = ''

    # 遍历DataFrame中的每一行
    for index, row in df.iterrows():
        inchi = row['MajorMicrospecies']
        # 检查InChI值是否为字符串
        if pd.notnull(inchi) and isinstance(inchi, str):
            chemical_formula, net_charge = get_molecular_details_from_inchi(inchi)
        else:
            chemical_formula, net_charge = 'Invalid InChI', 0
        df.at[index, 'Chemical_Formula'] = chemical_formula
        df.at[index, 'Net_Charge'] = net_charge

    # 将更新后的DataFrame保存到一个新的Excel文件
    df.to_excel(output_file, index=False, engine='openpyxl')

def get_molecular_details_from_inchi(inchi):
    mol = Chem.MolFromInchi(inchi)
    if mol is None:
        return 'Invalid InChI', 0
    chemical_formula = rdMolDescriptors.CalcMolFormula(mol)
    net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    return chemical_formula, net_charge

# 定义输入和输出文件路径
input_excel_path = '/mnt/NFS/fengch/aftercxcalcadd.xlsx'  # 更新为您的输入Excel文件路径
output_excel_path = '/mnt/NFS/fengch/afterrdkitadd.xlsx'  # 更新为您希望保存结果的Excel文件路径

# 调用函数处理Excel文件
process_excel(input_excel_path, output_excel_path)

