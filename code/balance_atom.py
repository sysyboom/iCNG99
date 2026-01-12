from cobra.io import read_sbml_model
import pandas as pd
from collections import Counter
import re


def parse_formula(formula):
    """解析化学式，返回一个包含原子计数的字典"""
    elements = Counter()
    matches = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    for (element, count) in matches:
        elements[element] += int(count) if count else 1
    return elements


def calculate_element_difference(reactant_elements, product_elements):
    """计算反应物和产物元素计数的差异"""
    difference = Counter(reactant_elements)
    for element, count in product_elements.items():
        difference[element] -= count
    return difference


def check_reaction_balance(reaction):
    """检查反应的原子守恒性，返回不平衡的详细信息"""
    reactant_elements = Counter()
    product_elements = Counter()

    for met, coeff in reaction.metabolites.items():
        formula = met.formula
        if formula:
            elements_count = parse_formula(formula)
            if coeff < 0:  # 反应物
                reactant_elements += Counter({el: count * -coeff for el, count in elements_count.items()})
            else:  # 产物
                product_elements += Counter({el: count * coeff for el, count in elements_count.items()})

    if reactant_elements != product_elements:
        element_difference = calculate_element_difference(reactant_elements, product_elements)
        return {
            'Reaction ID': reaction.id,
            'Reactant Elements': dict(reactant_elements),
            'Product Elements': dict(product_elements),
            'Element Difference': dict(element_difference)
        }
    return None  # 平衡时返回None


# 读取模型
model = read_sbml_model('/mnt/NFS/fengch/new/models/merge40001_noc.xml')

# 检查并收集不平衡的反应详细信息
unbalanced_reactions_info = []
for reaction in model.reactions:
    balance_info = check_reaction_balance(reaction)
    if balance_info:
        unbalanced_reactions_info.append(balance_info)

# 转换为DataFrame
df = pd.DataFrame(unbalanced_reactions_info)
excel_path = '/mnt/NFS/fengch/atom_merge40001_noc.xlsx'  # 确保路径正确
df.to_excel(excel_path, index=False)
