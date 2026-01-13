import pandas as pd

# 读取基因关联表格
gene_associations = pd.read_excel('/mnt/NFS/fengch/TPM/gene_ass.xlsx')

# 读取基因读数表格
gene_readings = pd.read_excel('/mnt/NFS/fengch/TPM/drug/dr.xlsx', sheet_name='TPM' )
#gene_readings = pd.read_csv('/mnt/NFS/fengch/new/transcriptome/new/drug_rawdata.csv')

# 计算每一行基因关联的读数
gene_averages = []
for association in gene_associations.iloc[:, 0]:
    # 先分割 'or'，并对每组 'and' 进行处理
    or_groups = association.split(' or ')
    row_values = []

    for group in or_groups:
        # 对每组 'and' 连接的基因取最小值
        associated_genes = group.split(' and ')
        gene_values = []

        for gene in associated_genes:
            # 获取基因对应的读数
            gene_values_in_row = gene_readings[gene_readings['gene'] == gene].iloc[:, 1:].values
            if len(gene_values_in_row) > 0:
                gene_values.append(gene_values_in_row[0])  # 取该基因对应行的值

        # 如果存在值，取基因组内的最小值，否则为 0
        if len(gene_values) > 0:
            # 取每列的最小值作为这一组的值
            row_values.append([min(values) for values in zip(*gene_values)])
        else:
            row_values.append([0] * (gene_readings.shape[1] - 1))  # 如果没有值，填充0

    # 对每个 'or' 组取平均
    if len(row_values) > 0:
        column_averages = [sum(column) / len(column) for column in zip(*row_values)]
    else:
        column_averages = [0] * (gene_readings.shape[1] - 1)
    gene_averages.append(column_averages)

# 将结果添加到基因关联表格中，使用原有的列名
original_column_names = gene_readings.columns[1:]  # 获取原有的列名，跳过第一列“Gene”
for i, column_name in enumerate(original_column_names):
    gene_associations[column_name] = [row[i] for row in gene_averages]

# 保存结果到新的 Excel 文件中
gene_associations.to_excel('/mnt/NFS/fengch/TPM/drug/drug1_reactions.xlsx', index=False)
