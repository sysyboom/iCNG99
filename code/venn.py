import pandas as pd
from matplotlib import pyplot as plt
from matplotlib_venn import venn3

# 读取Excel文件中的gene列
def read_excel_gene_column(file_path):
    df = pd.read_excel(file_path, engine='openpyxl')  # 确保指定engine='openpyxl'以读取xlsx文件
    return set(df['gene'].dropna())  # 删除任何NaN值

# 读取三个文件
genes_set1 = read_excel_gene_column('/mnt/NFS/fengch/new/draft/annotation_kegg.xlsx')
genes_set2 = read_excel_gene_column('/mnt/NFS/fengch/new/draft/homology_kegg.xlsx')
genes_set3 = read_excel_gene_column('/mnt/NFS/fengch/new/draft/homology_SamPler.xlsx')

# 准备数据集合
data = {
    "只在File 1中的genes": pd.Series(list(genes_set1 - genes_set2 - genes_set3), dtype=object),
    "只在File 2中的genes": pd.Series(list(genes_set2 - genes_set1 - genes_set3), dtype=object),
    "只在File 3中的genes": pd.Series(list(genes_set3 - genes_set1 - genes_set2), dtype=object),
    "在File 1和File 2中的genes": pd.Series(list(genes_set1 & genes_set2 - genes_set3), dtype=object),
    "在File 1和File 3中的genes": pd.Series(list(genes_set1 & genes_set3 - genes_set2), dtype=object),
    "在File 2和File 3中的genes": pd.Series(list(genes_set2 & genes_set3 - genes_set1), dtype=object),
    "在所有三个文件中的genes": pd.Series(list(genes_set1 & genes_set2 & genes_set3), dtype=object)
}

# 转换为DataFrame
df = pd.DataFrame(data)

# 写入Excel文件
output_file = '/mnt/NFS/fengch/new/draft/results_genes.xlsx'
df.to_excel(output_file, index=False)

print(f"集合的内容已经写入到 {output_file}。")

# 定义RGB值
rgb_color_set1 = (67, 170, 139)   # 红色
rgb_color_set2 = (39, 125, 161)   # 绿色
rgb_color_set3 = (144, 190, 109)   # 蓝色

# 将RGB值归一化到0到1之间
color_set1 = tuple(x / 255.0 for x in rgb_color_set1)
color_set2 = tuple(x / 255.0 for x in rgb_color_set2)
color_set3 = tuple(x / 255.0 for x in rgb_color_set3)
# 绘制韦恩图
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

#
for text in v.subset_labels:
    if text:  # 检查文本对象是否存在
        text.set_fontsize(14)
# 手动添加标签到图外
plt.text(-0.8, 0.8, 'annotation_kegg', fontsize=12)
plt.text(0.8, 0.8, 'homology_kegg', fontsize=12)
plt.text(0.8, 0.8, 'homology_SamPler', fontsize=12)

plt.show()
