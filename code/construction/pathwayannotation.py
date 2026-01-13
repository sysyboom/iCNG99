import xlrd
import xlwt
import requests


def get_pathway_name(pathway_id):
    url = f"https://rest.kegg.jp/get/{pathway_id}"
    response = requests.get(url)

    if response.status_code == 200:
        lines = response.text.strip().split("\n")
        for line in lines:
            if line.startswith("NAME"):
                _, pathway_name = line.split(" ", 1)
                return pathway_name.strip()

    return "无法找到通路名称"


def get_pathway_by_reaction(reaction_id):
    url = f"https://rest.kegg.jp/link/pathway/{reaction_id}"
    response = requests.get(url)

    if response.status_code == 200:
        lines = response.text.strip().split("\n")
        pathways = []
        for line in lines:
            fields = line.split("\t")
            if len(fields) >= 2:
                _, pathway_id = fields[:2]
                pathway_name = get_pathway_name(pathway_id)
                pathways.append(pathway_name)
        return pathways
    else:
        return None


# 读取输入的 .xls 文件
input_file = "/mnt/NFS/fengch/new/gapfill/pathwayanno.xls"  # 替换为你的输入文件路径
workbook = xlrd.open_workbook(input_file)
sheet = workbook.sheet_by_index(0)

# 创建输出的 .xls 文件
output_file = "/mnt/NFS/fengch/new/gapfill/res.xls"  # 替换为你的输出文件路径
output_workbook = xlwt.Workbook()
output_sheet = output_workbook.add_sheet("Output")

# 遍历每一行，获取反应ID并查询通路信息
for row in range(sheet.nrows):
    reaction_id = sheet.cell_value(row, 0)
    pathways = get_pathway_by_reaction(reaction_id)
    if pathways:
        pathway_str = ', '.join(pathways)
        output_sheet.write(row, 0, reaction_id)
        output_sheet.write(row, 1, pathway_str)
    else:
        output_sheet.write(row, 0, reaction_id)
        output_sheet.write(row, 1, "无法找到通路信息")

# 保存输出的 .xls 文件
output_workbook.save(output_file)

print("输出文件已保存成功。")