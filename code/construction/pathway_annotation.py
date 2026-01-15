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

    return "can't find"


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


input_file = "/mnt/NFS/fengch/new/gapfill/pathwayanno.xls" 
workbook = xlrd.open_workbook(input_file)
sheet = workbook.sheet_by_index(0)

output_file = "/mnt/NFS/fengch/new/gapfill/res.xls"  
output_workbook = xlwt.Workbook()
output_sheet = output_workbook.add_sheet("Output")

for row in range(sheet.nrows):
    reaction_id = sheet.cell_value(row, 0)
    pathways = get_pathway_by_reaction(reaction_id)
    if pathways:
        pathway_str = ', '.join(pathways)
        output_sheet.write(row, 0, reaction_id)
        output_sheet.write(row, 1, pathway_str)
    else:
        output_sheet.write(row, 0, reaction_id)
        output_sheet.write(row, 1, "can't find")

output_workbook.save(output_file)
print("done")
