import cobra
import pandas as pd

# 加载模型
model = cobra.io.read_sbml_model("/mnt/NFS/fengch/TPM/drug/Gal_GIMME.xml")

# 要修改的反应列表及其新的上下限
# reactions_to_modify = {
    #"EX_003": {"lower_bound": -1000, "upper_bound": 1000},  # 示例：葡萄糖交换反应
    #"EX_111": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_112": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_113": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_114": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_116": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_117": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_118": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_119": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_120": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_121": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_122": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_123": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_124": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_125": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_126": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_127": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_128": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_129": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_130": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_162": {"lower_bound": -1000, "upper_bound": 1000},
    #"biomass1": {"lower_bound": 5, "upper_bound": 1000},
#}

# 修改指定反应的上下限
# for reaction_id, bounds in reactions_to_modify.items():
    # reaction = model.reactions.get_by_id(reaction_id)
    # reaction.lower_bound = bounds["lower_bound"]
    # reaction.upper_bound = bounds["upper_bound"]

# 定义新的目标函数 - 使用占位符
new_objective_reaction_id = "biomass1"  # 示例：生物量生产反应
model.objective = new_objective_reaction_id

# 创建一个空的DataFrame来存储结果
gene_knockout_results = pd.DataFrame(columns=["Gene", "Wild-type Objective Value", "Knockout Objective Value"])

# 遍历模型中的每个基因
for gene in model.genes:
    gene_id = gene.id

    # 创建一个新的模型，不包括当前基因
    knockout_model = model.copy()
    knockout_model.genes.get_by_id(gene_id).knock_out()

    # 执行 FBA 分析以获得单基因敲除效应
    wildtype_solution = model.slim_optimize()
    knockout_solution = knockout_model.slim_optimize()

    # 将结果添加到DataFrame中
    gene_knockout_results = gene_knockout_results.append({"Gene": gene_id,
                                                          "Wild-type Objective Value": wildtype_solution,
                                                          "Knockout Objective Value": knockout_solution},
                                                         ignore_index=True)

# 根据敲除效应排序结果
gene_knockout_results["Objective Value Change"] = gene_knockout_results["Wild-type Objective Value"] - \
                                                  gene_knockout_results["Knockout Objective Value"]
gene_knockout_results = gene_knockout_results.sort_values(by="Objective Value Change", ascending=False)

# 将结果保存为CSV文件
gene_knockout_results.to_csv("/mnt/NFS/fengch/TPM/drug/knock_Gal.csv", index=False)


