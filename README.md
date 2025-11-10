# 🧬 Antibody–Antigen Structural Analysis Pipeline

该仓库提供了一整套用于 **抗体–抗原复合物分析** 的自动化脚本，包括：
- PDB 与 AF3 结构文件的统一命名；
- RMSD 计算；
- PDB–MSA–AF3 残基编号映射；
- PPI（Protein–Protein Interaction）结合位点分析；
- GD 与 AF3 结果重映射及重叠区域计算。

本项目旨在支持 **抗体结构预测模型（如 AlphaFold3）与实验结构** 的系统性比较分析。

---

## 🧱 项目结构

```
.
├── config.yaml                    # 全局配置文件（统一所有路径与参数）
├── rename.ipynb                   # 重命名 PDB/CIF 文件（确保文件名一致）
├── rmsd.py                        # 计算 RMSD，结构相似性评估
├── msa_index_map.ipynb            # PDB ↔ MSA ↔ AF3 残基编号映射
├── gd_results_ppi.ipynb           # 提取 GD 预测的抗原结合位点
├── af3_results_ppi.ipynb          # 提取 AF3 预测的抗原结合位点
├── overlap_analysis.ipynb         # 计算 GD 与 AF3 的结合位点重叠比例
├── data/
│   ├── IL23/
│   │   ├── 5njd.pdb              # 实验结构（Ground Truth）
│   │   ├── cif/                  # AF3 模型预测 CIF 文件夹
│   │   ├── cif_pdb/              # 转换后的 PDB 文件夹
│   │   ├── msa_fasta/            # MSA 输入序列
│   │   ├── msa_aln/              # 对齐后的 MSA 结果
│   │   ├── ppi_csv_gd/           # GD 的 PPI 分析结果
│   │   ├── ppi_csv_af3/          # AF3 的 PPI 分析结果
│   │   ├── output/               # RMSD 输出文件
│   │   ├── pdb_to_af3_mapping.csv# PDB–AF3 残基映射表
│   │   ├── ppi_summary_gd.csv    # GD PPI 汇总
│   │   ├── ppi_summary_af3.csv   # AF3 PPI 汇总
│   │   ├── af3_gd_overlap.csv    # 简单重叠结果
│   │   └── af3_gd_overlap_mapped.csv # 基于残基映射的重叠结果
└── README.md
```

---

## ⚙️ 全局配置文件说明 (`config.yaml`)

仓库中所有脚本均通过统一的配置文件 `config.yaml` 管理路径与参数，修改后会自动生效。


---

## 🧩 核心流程

### 1️⃣ Rename 阶段  
确保 `PDB` 与 `CIF` 文件命名一致，便于后续批量分析。

```bash
jupyter nbconvert --to notebook --execute rename.ipynb
```



---

### 2️⃣ RMSD 计算  

在统一命名后，计算实验结构与预测模型的 **VH/VL RMSD**，用于衡量结构相似度。

```bash
python rmsd.py
```

输出：  
`output/VHVL_RMSD_results.csv`，示例：

| model_name | VH_RMSD | VL_RMSD | average |
|-------------|----------|----------|----------|
| MJ00D_0     | 1.82     | 2.14     | 1.98     |

> ⚠️ RMSD 阈值设定为 **60 Å**，超过该值将被视为异常模型。

---

### 3️⃣ MSA/Index Mapping  

建立 **PDB–MSA–AF3 残基编号映射表**，解决多序列比对中 `-`（gap） 引起的编号错位问题。

```bash
jupyter nbconvert --to notebook --execute msa_index_map.ipynb
```

输出：`pdb_to_af3_mapping.csv`，示例：

| pdb_residue_id | msa_index | af3_residue_id |
|----------------|------------|----------------|
| 36             | 10         | 37             |
| 37             | 11         | 38             |
| ...            | ...        | ...            |

---

### 4️⃣ PPI 分析  

从 GD 与 AF3 输出的 CSV 文件中提取抗原结合位点信息。

```bash
jupyter nbconvert --to notebook --execute gd_results_ppi.ipynb
jupyter nbconvert --to notebook --execute af3_results_ppi.ipynb
```

输出示例：

**GD结果**
```csv
model_name,antigen_binding_sites_gd
5njd,"36,37,38,39,40,107,108,109,110,111,115,124,125,126,219"
```

**AF3结果**
```csv
model_name,antigen_binding_sites_AF3
MJ00D_0,"37,38,39,40,69,78,80,81,82,83,84,106,108,109,115,217,218,219"
```

---

### 5️⃣ Overlap 分析  

利用残基映射表计算 GD 与 AF3 的结合位点重叠率：

```bash
jupyter nbconvert --to notebook --execute overlap_analysis.ipynb
```

输出：
- `af3_gd_overlap_mapped.csv`：基于 PDB–AF3 编号映射后的精确比较。

---

## 📊 输出文件汇总

| 文件 | 描述 |
|------|------|
| `VHVL_RMSD_results.csv` | RMSD 计算结果 |
| `pdb_to_af3_mapping.csv` | 残基编号映射表 |
| `ppi_summary_gd.csv` | GD PPI 汇总 |
| `ppi_summary_af3.csv` | AF3 PPI 汇总 |
| `af3_gd_overlap_mapped.csv` | 基于映射的重叠结果 |

---

## 💡 注意事项

- 所有路径均可通过 `config.yaml` 一处修改；
- 运行顺序建议：
  ```
  rename → rmsd → msa_index_map → gd_results_ppi → af3_results_ppi → overlap_analysis
  ```
- MSA 阶段若出现 gap（`-`），脚本会自动跳过并打印警告；
- 可使用多线程 (`N_JOBS`) 提升运行速度；
- 若不使用 PyMOL，请确保 `NO_PYMOL: true`。

---

## ✍️ 作者与维护

- **Author:** Yuyang  
- **Institution:** SJTU 
- **Contact:** ice292@sjtu.edu.cn 
- **Last Updated:** 2025-11-10  


