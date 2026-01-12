import pandas as pd
from collections import Counter
import numpy as np

# 读取四个Excel文件
df1 = pd.read_excel('/mnt/NFS/fengch/new/location/deeploc.xlsx')
df2 = pd.read_excel('/mnt/NFS/fengch/new/location/yloc.xlsx')
df3 = pd.read_excel('/mnt/NFS/fengch/new/location/wolf.xlsx')
df4 = pd.read_excel('/mnt/NFS/fengch/new/location/yloc2.xlsx')  # 新添加的Excel文件

# 初始化一个空的字典来存储最终的投票结果
final_results = {}

# 对每一个Excel文件进行处理
for df in [df1, df2, df3]:
    # 对每一行进行处理
    for index, row in df.iterrows():
        # 获取蛋白质名称和位置信息
        protein_name = row[0]
        locations = row[1:]

        # 如果蛋白质名称还没有在最终结果中，则添加它
        if protein_name not in final_results:
            final_results[protein_name] = Counter()

        # 对每一个位置进行投票
        for location in locations:
            # 检查位置信息是否存在
            if pd.isnull(location) or location == '':
                continue
            # 位置的每次出现投两票
            final_results[protein_name][location] += 1

# 对新添加的Excel文件中的位置进行投票，每个位置加0.5分
for index, row in df4.iterrows():
    protein_name = row[0]
    locations = row[1:]
    if protein_name in final_results:
        for location in locations:
            if pd.isnull(location) or location == '':
                continue
            final_results[protein_name][location] += 0.5

# 将投票结果转换为列表并进行排序
sorted_final_list = []
for protein_name, votes in final_results.items():
    # 按投票数对位置进行排序
    sorted_votes = sorted(votes.items(), key=lambda x: x[1], reverse=True)
    row_data = [protein_name] + [item for sublist in sorted_votes for item in sublist]
    sorted_final_list.append(row_data)

# 找出最大的列数
max_columns = max(len(r) for r in sorted_final_list)

# 创建列名
columns = ['Protein'] + [f'Location{i // 2 + 1}' if i % 2 == 0 else f'Votes{i // 2 + 1}' for i in range(1, max_columns)]

# 将列表转换为数据框
sorted_results_df = pd.DataFrame(sorted_final_list, columns=columns)

# 保存到Excel文件
sorted_results_df.to_excel('/mnt/NFS/fengch/new/location/locationcombined_sorted.xlsx', index=False)
