import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from bioblocks.io import read_model
from bioblocks.geometry import compute_rmsd
from bioblocks.align import usalign_models
from bioblocks.transform import filter_atoms, filter_chains
from joblib import Parallel, delayed
from tqdm import tqdm

sns.set_theme(style="whitegrid")

# ========== CONFIG ==========
GROUND_TRUTH_PATH = "/home/yuyang/lb/data/IL23/5njd.pdb"   # <-- 单个 ground truth 文件（pdb 或 cif）
PREDICTED_DIR = "/home/yuyang/lb/data/IL23/cif"          # <-- 目录，里面是多个 .cif (预测结构)
OUTPUT_DIR = "/home/yuyang/lb/data/IL23"    # <-- 结果输出目录
os.makedirs(OUTPUT_DIR, exist_ok=True)

N_JOBS = 4  # 并行线程数，视你的机器调整

# ========== FUNCTIONS ==========

def filter_bb_and_chains(model, chains):
    """保留指定链与 backbone 原子 (N, CA, C, O)"""
    model = filter_chains(model, chain_filter=lambda c: c.id in chains)
    model = filter_atoms(model, atom_filter=lambda a: a.element in ["N", "CA", "C", "O"])
    return model


def compute_vh_vl_rmsd(file, truth_model):
    """
    对单个 predicted 文件：
      - 读取 predicted (file)
      - 在 A 链上把 predicted 对齐到 truth_model
      - 计算 H+L backbone RMSD（使用 bioblocks.compute_rmsd）
    返回 (file, rmsd_or_None)
    """
    pred_path = os.path.join(PREDICTED_DIR, file)
    try:
        pred = read_model(pred_path)
    except Exception as e:
        print(f"[READ ERROR] {file}: {e}")
        return file, None

    try:
        # Align predicted -> truth using antigen chain "A"
        aligned_pred = usalign_models(
            moving_model=pred,
            static_model=truth_model,
            alignment_chain_ids=["A", "A"]
        )

        # Filter H,L chains and backbone atoms
        pred_bb = filter_bb_and_chains(aligned_pred, ["H", "L"])
        truth_bb = filter_bb_and_chains(truth_model, ["H", "L"])

        rmsd_val = compute_rmsd(pred_bb, truth_bb)
        return file, rmsd_val

    except Exception as e:
        print(f"[PROCESS ERROR] {file}: {e}")
        return file, None


def main():
    # 检查 ground truth 文件存在性并读取一次
    if not os.path.exists(GROUND_TRUTH_PATH):
        raise FileNotFoundError(f"Ground truth file not found: {GROUND_TRUTH_PATH}")
    print(f"Loading ground truth model from: {GROUND_TRUTH_PATH}")
    truth_model = read_model(GROUND_TRUTH_PATH)

    # 列出预测文件（只挑 .cif，可按需扩展）
    cif_files = sorted([f for f in os.listdir(PREDICTED_DIR) if f.lower().endswith(".cif")])
    print(f"Found {len(cif_files)} predicted files in {PREDICTED_DIR}.")

    # 并行/串行处理
    results = Parallel(n_jobs=N_JOBS, verbose=1)(
        delayed(compute_vh_vl_rmsd)(file, truth_model) for file in tqdm(cif_files)
    )

    # 存 CSV（保留 None，以便 later debugging）
    df = pd.DataFrame(results, columns=["file", "RMSD"])
    # 根据 RMSD 从小到大排序（最小在前），将缺失值放到最后，重置索引
    df = df.sort_values(by="RMSD", ascending=True, na_position="last").reset_index(drop=True)
    output_csv = os.path.join(OUTPUT_DIR, "VHVL_RMSD_results.csv")
    df.to_csv(output_csv, index=False)
    print(f"Saved RMSD results to {output_csv}")

    # 只绘制有数值的分布图
    df_valid = df[df["RMSD"].notna()]
    if not df_valid.empty:
        plt.figure(figsize=(6, 4))
        sns.histplot(df_valid["RMSD"], bins=50, kde=True)
        plt.title("VH+VL RMSD distribution")
        plt.xlabel("RMSD (Å)")
        plt.ylabel("Count")
        plt.tight_layout()
        output_png = os.path.join(OUTPUT_DIR, "VHVL_RMSD_dist.png")
        plt.savefig(output_png, dpi=300)
        plt.close()
        print(f"Saved RMSD distribution plot to {output_png}")
    else:
        print("No valid RMSD values to plot (all None).")


if __name__ == "__main__":
    main()