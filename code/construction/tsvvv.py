import pandas as pd


def extract_and_save(filename, output_file):
    # 读取数据
    data = pd.read_csv(filename, sep='\t')

    # 定义一个函数来从 'xrefs' 列中提取所需的信息
    def extract_all_references(info, keyword):
        # 使用 ' ... ' 分割字符串
        items = info.split(' ... ')
        # 收集所有包含关键词的项
        matches = [item for item in items if keyword in item]
        if matches:
            return ' ... '.join(matches)  # 连接所有找到的匹配项
        return "Not found"  # 如果没有找到，返回 'Not found'

    # 应用函数并创建新列
    data['bigg.reaction'] = data['xrefs'].apply(lambda x: extract_all_references(x, 'bigg.reaction'))
    data['seed.reaction'] = data['xrefs'].apply(lambda x: extract_all_references(x, 'seed.reaction'))
    data['metacyc.reaction'] = data['xrefs'].apply(lambda x: extract_all_references(x, 'metacyc.reaction'))
    #data['seed.compound'] = data['xrefs'].apply(lambda x: extract_all_references(x, 'seed.compound'))
    #data['reactome:'] = data['xrefs'].apply(lambda x: extract_all_references(x, 'reactome:'))

    # 创建新的DataFrame，包含所需的列
    output_data = data[['#query', 'mnx_id', 'bigg.reaction', 'seed.reaction','metacyc.reaction']]
    # 假设第一列和第二列的列名是 'Column1' 和 'Column2'

    # 输出到Excel文件
    output_data.to_excel(output_file, index=False)


# 使用示例
filename = '/mnt/NFS/fengch/new/annotationall/id-mapper43'  # 替换为你的文件路径
output_file = '/mnt/NFS/fengch/new/annotationall/id-mapper43.xlsx'  # 输出Excel文件的路径
extract_and_save(filename, output_file)
