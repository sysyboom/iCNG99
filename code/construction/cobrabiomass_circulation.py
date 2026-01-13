import cobra
import pandas as pd

# 需要依次修改的交换反应列表
# exchange_reactions_group1 = ['EX_022', 'EX_032', 'EX_033', 'EX_035', 'EX_037','EX_038','EX_039']
# exchange_reactions_group2 = ['EX_023', 'EX_024', 'EX_025', 'EX_026', 'EX_027','EX_028','EX_029','EX_030','EX_031','EX_034','EX_036']
# Alternative nitrogen sources
exchange_reactions_group1 = []
exchange_reactions_group2 = ['EX_040', 'EX_111', 'EX_112', 'EX_113', 'EX_114','EX_115','EX_116','EX_117','EX_118','EX_119','EX_120','EX_121','EX_122','EX_123','EX_124','EX_125','EX_126','EX_127','EX_128','EX_129','EX_130']
# 固定要修改的反应ID
fixed_reaction_id = 'EX_015'  # 例如 EX_020 为固定的交换反应

# 目标反应的ID对应
biomass1_reaction_id = 'biomass2_c' # with capsule
biomass2_reaction_id = 'biomass1_c' # no capsule

# 用于保存结果的列表
results = []

# 进行循环，处理第一组交换反应，目标函数为 biomass1
for i in range(len(exchange_reactions_group1)):
    print(f"Starting iteration {i + 1} for Group 1 (biomass1)...")

    # 每次循环重新加载模型
    model = cobra.io.read_sbml_model('/mnt/NFS/fengch/new/models/mergemodeltest40001.xml')
    model.solver = 'cplex'

    # 获取当前要修改的交换反应ID
    current_reaction_id = exchange_reactions_group1[i]

    # 修改当前交换反应的上下限
    reaction = model.reactions.get_by_id(current_reaction_id)
    reaction.lower_bound = -1000
    reaction.upper_bound = 1000

    # 修改固定的交换反应的上下限
    fixed_reaction = model.reactions.get_by_id(fixed_reaction_id)
    fixed_reaction.lower_bound = 0
    fixed_reaction.upper_bound = 0

    # 设置目标函数为 biomass1
    target_reaction = model.reactions.get_by_id(biomass1_reaction_id)
    model.objective = target_reaction
    model.objective_direction = 'max'

    # 优化模型
    solution = model.optimize()

    # 将结果保存到列表中
    results.append({
        'Iteration': i + 1,
        'Group': 'Group 1',
        'Modified_Reaction': current_reaction_id,
        'Fixed_Reaction': fixed_reaction_id,
        'Objective_Function': biomass1_reaction_id,
        'Optimized_Value': solution.objective_value
    })

    print(f'Iteration {i + 1} - {current_reaction_id} and {fixed_reaction_id} modified:')
    print('Optimized value:', solution.objective_value)

# 进行循环，处理第二组交换反应，目标函数为 biomass2
for i in range(len(exchange_reactions_group2)):
    print(f"Starting iteration {i + 1} for Group 2 (biomass2)...")

    # 每次循环重新加载模型
    model = cobra.io.read_sbml_model('/mnt/NFS/fengch/new/models/mergemodeltest40001.xml')
    model.solver = 'cplex'

    # 获取当前要修改的交换反应ID
    current_reaction_id = exchange_reactions_group2[i]

    # 修改当前交换反应的上下限
    reaction = model.reactions.get_by_id(current_reaction_id)
    reaction.lower_bound = -1000
    reaction.upper_bound = 1000

    # 修改固定的交换反应的上下限
    fixed_reaction = model.reactions.get_by_id(fixed_reaction_id)
    fixed_reaction.lower_bound = 0
    fixed_reaction.upper_bound = 0

    # 设置目标函数为 biomass2
    target_reaction = model.reactions.get_by_id(biomass2_reaction_id)
    model.objective = target_reaction
    model.objective_direction = 'max'

    # 优化模型
    solution = model.optimize()

    # 将结果保存到列表中
    results.append({
        'Iteration': i + 1,
        'Group': 'Group 2',
        'Modified_Reaction': current_reaction_id,
        'Fixed_Reaction': fixed_reaction_id,
        'Objective_Function': biomass2_reaction_id,
        'Optimized_Value': solution.objective_value
    })

    print(f'Iteration {i + 1} - {current_reaction_id} and {fixed_reaction_id} modified:')
    print('Optimized value:', solution.objective_value)

# 将结果转换为pandas DataFrame
df_results = pd.DataFrame(results)

# 将DataFrame保存到Excel文件
df_results.to_excel('/mnt/NFS/fengch/new/transport/nitrogen_biomass_results.xlsx', index=False)

print("Results have been saved to 'results.xlsx'.")
