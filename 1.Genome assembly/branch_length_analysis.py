import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from Bio import Phylo
from io import StringIO

# 设置 Seaborn 样式
sns.set_theme(style="whitegrid", font_scale=1.2)
sns.set_palette("colorblind") 

# 1. 定义输入文件和目标物种
tree_file = "all_trees_ufboot.tre"
target_taxa = ["HnothoanatolicumXam", "HnothoanatolicumXag", "Hmarinum", "Hgussoneanum"]

# 用于存储数据的列表
data = []

# 2. 读取树文件并提取端点分支长度
print("正在使用 Biopython 解析树文件提取分支长度...")
with open(tree_file, "r") as f:
    for tree_idx, line in enumerate(f):
        line = line.strip()
        if not line:
            continue
        
        try:
            tree = Phylo.read(StringIO(line), "newick")
        except Exception as e:
            print(f"树 {tree_idx} 解析错误: {e}")
            continue
            
        tree_lengths = {"Gene_Tree": f"Tree_{tree_idx+1}"}
        all_taxa_found = True
        
        for taxon in target_taxa:
            clade = tree.find_any(name=taxon)
            if clade and clade.branch_length is not None:
                tree_lengths[taxon] = clade.branch_length
            else:
                all_taxa_found = False
                break
        
        if all_taxa_found:
            data.append(tree_lengths)

# 3. 转换为 DataFrame
df = pd.DataFrame(data)
print(f"成功提取了 {len(df)} 棵包含所有目标分类群的基因树。")

# 保存原始数据（可选）
df.to_csv("branch_lengths_raw.csv", index=False)

# 4. 统计检验 (完整数据)
print("\n--- 统计检验结果（完整数据）---")
stat_xam_xag, p_xam_xag = stats.wilcoxon(df["HnothoanatolicumXam"], df["HnothoanatolicumXag"])
print(f"Xam vs Xag -> Wilcoxon stat: {stat_xam_xag:.2f}, p-value: {p_xam_xag:.2e}")

stat_xam_mar, p_xam_mar = stats.wilcoxon(df["HnothoanatolicumXam"], df["Hmarinum"])
print(f"Xam vs Hmarinum -> Wilcoxon stat: {stat_xam_mar:.2f}, p-value: {p_xam_mar:.2e}")

stat_xag_gus, p_xag_gus = stats.wilcoxon(df["HnothoanatolicumXag"], df["Hgussoneanum"])
print(f"Xag vs Hgussoneanum -> Wilcoxon stat: {stat_xag_gus:.2f}, p-value: {p_xag_gus:.2e}")

# 描述性统计（完整数据）
desc_stats = df[target_taxa].agg(['mean', 'median', 'std']).T
desc_stats.columns = ['Mean', 'Median', 'Std']
print("\n--- 完整数据描述性统计 ---")
print(desc_stats.round(6))

# ================== 绘图部分 ==================
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# ================== 1. 数据准备 ==================
threshold = 0.1
df_melt = pd.melt(df, id_vars=["Gene_Tree"], value_vars=target_taxa, 
                  var_name="Taxon", value_name="Branch_Length")
df_filtered = df_melt[df_melt["Branch_Length"] <= threshold].copy()

# ================== 2. 风格与配色配置 ==================
custom_palette = {
    "HnothoanatolicumXam": "#a14d61",  # Xam (实心深红)
    "HnothoanatolicumXag": "#8cb1ca",  # Xag (实心蓝灰)
    "Hmarinum": "#f5c765",             # H. marinum (实心淡黄)
    "Hgussoneanum": "#989898"          # H. gussoneanum (实心深灰)
}

# 解决 Linux 字体兼容
plt.rcParams['font.sans-serif'] = ['Arial', 'Liberation Sans', 'DejaVu Sans', 'sans-serif']
plt.rcParams['pdf.fonttype'] = 42

order_list = ["HnothoanatolicumXam", "HnothoanatolicumXag", "Hmarinum", "Hgussoneanum"]

# ================== 3. “硬核三明治”绘图逻辑 ==================
fig, ax = plt.subplots(figsize=(6, 8)) 

# --- 底层：颜色填充 (无边框) ---
sns.boxplot(
    x="Taxon", y="Branch_Length", data=df_filtered,
    hue="Taxon", palette=custom_palette, order=order_list,
    width=0.4, showfliers=False, legend=False,
    linewidth=0,  
    ax=ax
)

# --- 中层：实心散点 (无透明度，黑边加粗) ---
sns.stripplot(
    x="Taxon", y="Branch_Length", data=df_filtered,
    hue="Taxon", palette=custom_palette, order=order_list,
    size=5,             # 稍微增大，更显饱满
    alpha=1.0,          # 关键修改：取消透明度，完全实心
    jitter=0.12,        
    edgecolor="black",  # 黑色描边
    linewidth=0.8,      # 关键修改：描边再粗一号 (由0.3增至0.8)
    legend=False, 
    zorder=2,           # 确保在填充色之上
    ax=ax
)

# --- 顶层：超粗黑色外框架 (压在最上面) ---
sns.boxplot(
    x="Taxon", y="Branch_Length", data=df_filtered,
    hue="Taxon", order=order_list,
    width=0.4, showfliers=False, legend=False,
    fill=False,         
    color="black",      
    zorder=10,          # 确保框架压死散点边缘
    boxprops=dict(color="black", linewidth=1.5), 
    whiskerprops=dict(color="black", linewidth=2.5),
    capprops=dict(color="black", linewidth=2.5),
    medianprops=dict(color="black", linewidth=3.2), # 中位数线加厚，压在最顶
    ax=ax
)

# ================== 4. 细节修饰与保存 ==================
ax.set_ylabel("Terminal Branch Length (substitutions/site)", 
              fontsize=14, fontweight='bold', labelpad=15)
ax.set_xlabel("") 

# 修复 Tick 警告
ax.set_xticks(range(len(order_list)))
nice_labels = ["Xam\n(Sub)", "Xag\n(Sub)", "H. mar\n(Anc)", "H. guss\n(Anc)"]
ax.set_xticklabels(nice_labels, fontsize=11, fontweight='bold')

# 去掉上右边框，加厚刻度线
sns.despine(offset=10, trim=True)
ax.tick_params(width=2)

# Y轴自适应留白
ymax = df_filtered["Branch_Length"].max()
ax.set_ylim(-0.001, ymax * 1.15)

plt.tight_layout()
plt.savefig("branch_lengths_solid_style.png", dpi=600, bbox_inches="tight")
plt.savefig("branch_lengths_solid_style.pdf", bbox_inches="tight")

print("\n图表已完成修改：散点现为实心且带粗黑描边。")
