import os
import yaml
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

# ===================== YAML 参数读取 =====================
def load_params(yaml_file="config.yaml"):
    """读取 config.yaml 参数"""
    with open(yaml_file, "r") as f:
        config = yaml.safe_load(f)
    common_params = config.get("common", {})
    rmsd_params = config.get("rmsd", {})
    params = {**common_params, **rmsd_params}
    return params

# ===================== 函数定义 =====================
def filter_bb_and_chains(model, chains):
    """保留指定链与 backbone 原子 (N, CA, C, O)"""
    model = filter_chains(model, chain_filter=lambda c: c.id in chains)
    model = filter_atoms(model, atom_filter=lambda a: a.element in ["N", "CA", "C", "O"])
    return model

def compute_vh_vl_rmsd(file, truth_model, predicted_dir):
    """对单个 predicted 文件计算 H+L backbone RMSD"""
    pred_path = os.path.join(predicted_dir, file)
    try:
        pred = read_model(pred_path)
    except Exception as e:
        print(f"[READ ERROR] {file}: {e}")
        return file, None

    try:
        aligned_pred = usalign_models(
            moving_model=pred,
            static_model=truth_model,
            alignment_chain_ids=["A", "A"]
        )

        pred_bb_H = filter_bb_and_chains(aligned_pred, ["H"])
        truth_bb_H = filter_bb_and_chains(truth_model, ["H"])

        pred_bb_L = filter_bb_and_chains(aligned_pred, ["L"])
        truth_bb_L = filter_bb_and_chains(truth_model, ["L"])

        rmsd_val_H = compute_rmsd(pred_bb_H, truth_bb_H)
        rmsd_val_L = compute_rmsd(pred_bb_L, truth_bb_L)

        rmsd_val = rmsd_val_H + rmsd_val_L
        return file, rmsd_val

    except Exception as e:
        print(f"[PROCESS ERROR] {file}: {e}")
        return file, None

# ===================== 主函数 =====================
def main(yaml_file="config.yaml"):
    params = load_params(yaml_file)

    GROUND_TRUTH_PATH = params["GROUND_TRUTH_PATH"]
    PREDICTED_DIR = params["PREDICTED_DIR"]
    OUTPUT_DIR = params["OUTPUT_DIR"]
    N_JOBS = params.get("N_JOBS", 4)
    plot_bins = params.get("plot_bins", 50)
    plot_kde = params.get("plot_kde", True)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    if not os.path.exists(GROUND_TRUTH_PATH):
        raise FileNotFoundError(f"Ground truth file not found: {GROUND_TRUTH_PATH}")
    print(f"Loading ground truth model from: {GROUND_TRUTH_PATH}")
    truth_model = read_model(GROUND_TRUTH_PATH)

    cif_files = sorted([f for f in os.listdir(PREDICTED_DIR) if f.lower().endswith(".cif")])
    print(f"Found {len(cif_files)} predicted files in {PREDICTED_DIR}.")

    results = Parallel(n_jobs=N_JOBS, verbose=1)(
        delayed(compute_vh_vl_rmsd)(file, truth_model, PREDICTED_DIR) for file in tqdm(cif_files)
    )

    df = pd.DataFrame(results, columns=["file", "RMSD"])
    df = df.sort_values(by="RMSD", ascending=True, na_position="last").reset_index(drop=True)

    output_csv = os.path.join(OUTPUT_DIR, "VHVL_RMSD_results.csv")
    df.to_csv(output_csv, index=False)
    print(f"Saved RMSD results to {output_csv}")

    df_valid = df[df["RMSD"].notna()]
    if not df_valid.empty:
        plt.figure(figsize=(6, 4))
        sns.histplot(df_valid["RMSD"], bins=plot_bins, kde=plot_kde)
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
