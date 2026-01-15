import pandas as pd
from collections import Counter
import numpy as np

df1 = pd.read_excel('/mnt/NFS/fengch/new/location/deeploc.xlsx')
df2 = pd.read_excel('/mnt/NFS/fengch/new/location/yloc.xlsx')
df3 = pd.read_excel('/mnt/NFS/fengch/new/location/wolf.xlsx')
df4 = pd.read_excel('/mnt/NFS/fengch/new/location/yloc2.xlsx')  # 0.5

final_results = {}

for df in [df1, df2, df3]:
    for index, row in df.iterrows():
        protein_name = row[0]
        locations = row[1:]

        if protein_name not in final_results:
            final_results[protein_name] = Counter()

        for location in locations:
            if pd.isnull(location) or location == '':
                continue
            final_results[protein_name][location] += 1

for index, row in df4.iterrows():
    protein_name = row[0]
    locations = row[1:]
    if protein_name in final_results:
        for location in locations:
            if pd.isnull(location) or location == '':
                continue
            final_results[protein_name][location] += 0.5

sorted_final_list = []
for protein_name, votes in final_results.items():
    sorted_votes = sorted(votes.items(), key=lambda x: x[1], reverse=True)
    row_data = [protein_name] + [item for sublist in sorted_votes for item in sublist]
    sorted_final_list.append(row_data)

max_columns = max(len(r) for r in sorted_final_list)

columns = ['Protein'] + [f'Location{i // 2 + 1}' if i % 2 == 0 else f'Votes{i // 2 + 1}' for i in range(1, max_columns)]

sorted_results_df = pd.DataFrame(sorted_final_list, columns=columns)

sorted_results_df.to_excel('/mnt/NFS/fengch/new/location/locationcombined_sorted.xlsx', index=False)
