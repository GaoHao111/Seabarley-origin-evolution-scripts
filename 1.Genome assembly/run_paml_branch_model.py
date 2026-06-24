import os
import re
import subprocess
import pandas as pd
from pathlib import Path
from multiprocessing import Pool
from scipy.stats import chi2

# ================= 配置参数 =================
WORK_DIR = Path("/media/langchao/34e000e7-97b5-46f1-a86f-71b5b6c91ae4/SQF/GH/ComparativeGenomics/SubgenomeDominance/CodeML/PAML_Runs")
TREE_SRC_DIR = Path("/media/langchao/34e000e7-97b5-46f1-a86f-71b5b6c91ae4/SQF/GH/ComparativeGenomics/PhylogeneticTree/pep/OrthoFinder/Results_Oct30/GeneTrees")
THREADS = 16
OUTPUT_CSV = WORK_DIR / "BranchModel_Standard_Results.csv"

def prepare_labeled_tree(tree_src, tree_dst):
    """清理并标注树文件"""
    if not tree_src.exists(): return False
    with open(tree_src) as f:
        tree_str = f.read().strip()
    
    # 清理支持度和格式
    tree_str = re.sub(r'\)\d+(\.\d+)?(/\d+(\.\d+)?)?:', '):', tree_str)
    tree_str = tree_str.split(';')[-2] + ';' if ';' in tree_str else tree_str
    
    # 标注 Xam(#1) 和 Xag(#2)
    # 这里的正则要匹配 ID 后面的分支长度部分前插入标签
    # 寻找模式：ID:0.123 -> ID #1:0.123
    tree_str = re.sub(r'([^,()]+Xa_m[^,():\s]+)', r'\1 #1', tree_str)
    tree_str = re.sub(r'([^,()]+Xa_g[^,():\s]+)', r'\1 #2', tree_str)
    
    # 检查是否标注成功
    if "#1" not in tree_str or "#2" not in tree_str:
        return False

    with open(tree_dst, 'w') as f:
        f.write(tree_str + '\n')
    return True

def extract_omega_from_mlc(mlc_file):
    """专门从 codeml (seqtype=1) 的输出中提取 omega"""
    data = {"lnL": None, "w_bg": None, "w_Xam": None, "w_Xag": None}
    if not mlc_file.exists(): return data
    
    with open(mlc_file) as f:
        content = f.read()
        
        # 1. 提取 lnL
        lnl_match = re.search(r'lnL\(ntime:.*?\):\s+([-.\d]+)', content)
        if lnl_match: data["lnL"] = float(lnl_match.group(1))
        
        # 2. 提取 w (针对 Branch Model)
        # 在 model=2 中，PAML 会输出 "w ratios:  0.01 0.05 0.02"
        w_match = re.search(r'w ratios:\s+([\d.\s]+)', content)
        if w_match:
            ws = w_match.group(1).split()
            if len(ws) >= 3:
                data["w_bg"], data["w_Xam"], data["w_Xag"] = map(float, ws[:3])
            elif len(ws) == 1: # M0 情况
                data["w_bg"] = data["w_Xam"] = data["w_Xag"] = float(ws[0])
    return data

def run_paml_for_gene(gene):
    gene_dir = WORK_DIR / gene
    phy_file = gene_dir / f"{gene}.phy"
    tree_src = TREE_SRC_DIR / f"{gene}.treefile"
    tree_labeled = gene_dir / f"labeled.nwk"
    
    # 检查输入
    if not phy_file.exists() or not tree_src.exists(): return None

    # 1. 准备树文件
    if not prepare_labeled_tree(tree_src, tree_labeled):
        return {"Gene": gene, "Error": "Tree labeling failed"}

    # 2. 构建控制参数 (强制 seqtype = 1)
    # 注意：CodonFreq = 2 (F3x4) 是顶刊最常用的设置
    base_ctl = f"""
      seqfile = {phy_file.name}
     treefile = {tree_labeled.name}
      outfile = tmp.mlc
        noisy = 0 ; verbose = 0 ; runmode = 0
      seqtype = 1      * 1:codons; 2:AAs
    CodonFreq = 2      * 0:1/61, 1:F1X4, 2:F3X4, 3:F61
        icode = 0      * 0:universal code
    fix_kappa = 0 ; kappa = 2
    fix_omega = 0 ; omega = 0.5
    fix_alpha = 1 ; alpha = 0
    cleandata = 1      * 移除带Gap的位点
    """

    # --- 运行 M0 ---
    with open(gene_dir / "M0.ctl", 'w') as f:
        f.write(base_ctl + "\nmodel = 0\noutfile = M0.mlc")
    subprocess.run(["codeml", "M0.ctl"], cwd=gene_dir, capture_output=True)

    # --- 运行 M1 (Three-ratio) ---
    with open(gene_dir / "M1.ctl", 'w') as f:
        f.write(base_ctl + "\nmodel = 2\noutfile = M1.mlc")
    subprocess.run(["codeml", "M1.ctl"], cwd=gene_dir, capture_output=True)

    # 3. 提取数据
    m0_res = extract_omega_from_mlc(gene_dir / "M0.mlc")
    m1_res = extract_omega_from_mlc(gene_dir / "M1.mlc")
    
    if m0_res["lnL"] is None or m1_res["lnL"] is None:
        return {"Gene": gene, "Error": "PAML output empty or invalid"}

    # 计算 LRT
    lrt_stat = max(0, 2 * (m1_res["lnL"] - m0_res["lnL"]))
    p_val = 1 - chi2.cdf(lrt_stat, df=2)

    return {
        "Gene": gene,
        "lnL_M0": m0_res["lnL"],
        "lnL_M1": m1_res["lnL"],
        "P_value": p_val,
        "w_bg": m1_res["w_bg"],
        "w_Xam": m1_res["w_Xam"],
        "w_Xag": m1_res["w_Xag"],
        "Error": "None"
    }

if __name__ == "__main__":
    genes = [d.name for d in WORK_DIR.iterdir() if d.is_dir() and (d / f"{d.name}.phy").exists()]
    print(f"开始处理 {len(genes)} 个基因 (Codon Model)...")
    
    with Pool(THREADS) as pool:
        results = pool.map(run_paml_for_gene, genes)
    
    df = pd.DataFrame([r for r in results if r is not None])
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"分析完成。有效数据：{len(df[df['Error'] == 'None'])} 条。")
