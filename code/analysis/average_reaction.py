import pandas as pd

gene_associations = pd.read_excel('/mnt/NFS/fengch/TPM/gene_ass.xlsx')
gene_readings = pd.read_excel('/mnt/NFS/fengch/TPM/drug/dr.xlsx', sheet_name='TPM' )
#gene_readings = pd.read_csv('/mnt/NFS/fengch/new/transcriptome/new/drug_rawdata.csv')

gene_averages = []
for association in gene_associations.iloc[:, 0]:
    or_groups = association.split(' or ')
    row_values = []

    for group in or_groups:
        associated_genes = group.split(' and ')
        gene_values = []

        for gene in associated_genes:
            gene_values_in_row = gene_readings[gene_readings['gene'] == gene].iloc[:, 1:].values
            if len(gene_values_in_row) > 0:
                gene_values.append(gene_values_in_row[0])  

        if len(gene_values) > 0:
            row_values.append([min(values) for values in zip(*gene_values)])
        else:
            row_values.append([0] * (gene_readings.shape[1] - 1))  

    if len(row_values) > 0:
        column_averages = [sum(column) / len(column) for column in zip(*row_values)]
    else:
        column_averages = [0] * (gene_readings.shape[1] - 1)
    gene_averages.append(column_averages)

original_column_names = gene_readings.columns[1:]  
for i, column_name in enumerate(original_column_names):
    gene_associations[column_name] = [row[i] for row in gene_averages]

gene_associations.to_excel('/mnt/NFS/fengch/TPM/drug/drug1_reactions.xlsx', index=False)
