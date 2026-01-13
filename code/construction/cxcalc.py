import pandas as pd
import subprocess

# 输入和输出 Excel 文件路径
input_excel_path = '/mnt/NFS/fengch/beforecxcalc.xlsx'
output_excel_path = '/mnt/NFS/fengch/aftercxcalc.xlsx'

# 读取 Excel 文件
df = pd.read_excel(input_excel_path)

# 添加一列来存储计算结果
df['MajorMicrospecies'] = ''

# 循环处理每一行
for index, row in df.iterrows():
    smiles = row['InChI']  # 假设第二列的列名是'InChI'
    if smiles:
        # 构建 cxcalc 命令
        command = f'cxcalc majormicrospecies -H 7.2 "{smiles}"'

        try:
            # 执行命令
            result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            major_microspecies = result.stdout.strip()  # 获取主要微种

            # 将计算结果写入 DataFrame
            df.at[index, 'MajorMicrospecies'] = major_microspecies
        except subprocess.CalledProcessError as e:
            error_message = e.stderr.strip()
            print(f"Error executing command for SMILES {smiles}: {error_message}")

# 将结果保存到新的 Excel 文件
df.to_excel(output_excel_path, index=False)

print("计算完成并将结果写入 Excel 文件。")
