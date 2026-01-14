import pandas as pd
from matplotlib import pyplot as plt
from matplotlib_venn import venn3

def read_excel_gene_column(file_path):
    df = pd.read_excel(file_path, engine='openpyxl') 
    return set(df['gene'].dropna())  

genes_set1 = read_excel_gene_column('/mnt/NFS/fengch/new/draft/annotation_kegg.xlsx')
genes_set2 = read_excel_gene_column('/mnt/NFS/fengch/new/draft/homology_kegg.xlsx')
genes_set3 = read_excel_gene_column('/mnt/NFS/fengch/new/draft/homology_SamPler.xlsx')

data = {
    "File 1": pd.Series(list(genes_set1 - genes_set2 - genes_set3), dtype=object),
    "File 2": pd.Series(list(genes_set2 - genes_set1 - genes_set3), dtype=object),
    "File 3": pd.Series(list(genes_set3 - genes_set1 - genes_set2), dtype=object),
    "File 1 and File 2": pd.Series(list(genes_set1 & genes_set2 - genes_set3), dtype=object),
    "File 1 and File 3": pd.Series(list(genes_set1 & genes_set3 - genes_set2), dtype=object),
    "File 2 and File 3": pd.Series(list(genes_set2 & genes_set3 - genes_set1), dtype=object),
    "all three": pd.Series(list(genes_set1 & genes_set2 & genes_set3), dtype=object)
}

df = pd.DataFrame(data)

output_file = '/mnt/NFS/fengch/new/draft/results_genes.xlsx'
df.to_excel(output_file, index=False)

rgb_color_set1 = (67, 170, 139)   
rgb_color_set2 = (39, 125, 161)   
rgb_color_set3 = (144, 190, 109)   

color_set1 = tuple(x / 255.0 for x in rgb_color_set1)
color_set2 = tuple(x / 255.0 for x in rgb_color_set2)
color_set3 = tuple(x / 255.0 for x in rgb_color_set3)

plt.figure(figsize=(10, 8))
v = venn3(subsets=(
    len(genes_set1 - genes_set2 - genes_set3),
    len(genes_set2 - genes_set1 - genes_set3),
    len(genes_set1 & genes_set2 - genes_set3),
    len(genes_set3 - genes_set1 - genes_set2),
    len(genes_set1 & genes_set3 - genes_set2),
    len(genes_set2 & genes_set3 - genes_set1),
    len(genes_set1 & genes_set2 & genes_set3),
), set_labels=(' ' * len('annotation_kegg'), ' ' * len('homology_kegg'), ' ' * len('homology_SamPler')),
    set_colors=(color_set1, color_set2, color_set3))
for patch in v.patches:
    if patch:
        patch.set_alpha(0.6)

for text in v.subset_labels:
    if text:  
        text.set_fontsize(14)
plt.text(-0.8, 0.8, 'annotation_kegg', fontsize=12)
plt.text(0.8, 0.8, 'homology_kegg', fontsize=12)
plt.text(0.8, 0.8, 'homology_SamPler', fontsize=12)

plt.show()
