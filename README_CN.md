# SegDesign：智能蛋白质片段设计 pipeline

<div align="center">

**集序列分析、结构预测和生成建模于一体的智能蛋白质片段设计 pipeline**

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

</div>

## 📖 项目简介

SegDesign 是一个用于智能蛋白质片段设计的自动化 pipeline。它整合了多种先进的生物信息学工具和深度学习模型，可执行全面的蛋白质分析与设计：

- **序列保守性分析**：使用 HMMER 进行进化保守性分析
- **结构生成**：使用 RFdiffusion 进行靶向蛋白质骨架生成
- **序列设计**：使用 ProteinMPNN 进行氨基酸序列优化
- **序列聚类**：使用 MMSeqs2 进行序列相似性分析
- **结构验证**：使用 ESMFold 和 AlphaFold2 进行预测结构质量评估
- **二级结构分析**：使用 DSSP 进行蛋白质二级结构分析

## 🏗️ 项目架构

```
SegDesign/
├── Segdesign.py              # 主程序入口
├── Segdesign/
│   ├── hmmer/               # 序列保守性分析
│   ├── rfdiffusion/         # 结构生成
│   ├── mpnn/                # 序列设计
│   ├── esmfold/             # 结构预测
│   ├── alphafold2/          # AlphaFold2 结构预测
│   ├── dssp/                # 二级结构分析
│   └── mmseqs/              # 序列聚类
├── config/
│   ├── config.yaml          # 用户配置文件
│   └── setting.yaml         # 系统配置文件
├── environments/            # 环境安装脚本
└── example/                 # 示例输出
```

## 🚀 快速开始

### 前置条件

- **操作系统**：Linux（推荐）或 Windows+WSL2
- **Python**：3.9 或更高版本
- **Conda/Miniconda**：环境管理必需
- **GPU**：NVIDIA GPU 且支持 CUDA（强烈推荐用于 ESMFold 和 RFdiffusion）
- **内存**：至少 16GB RAM（推荐 32GB 以上）
- **存储**：至少 200GB 可用空间

### 安装步骤

#### 1. 克隆仓库

```bash
git clone https://github.com/mike114b/SegDesign.git
cd SegDesign
```

#### 2. 安装 Conda 环境

项目需要4个 conda 环境来运行不同的模块，分别是：
- **segdesign**：主环境，包含 HMMER、MMSeqs2、DSSP 等工具
- **segdesign_esmfold**：用于运行 ESMFold 模型
- **segdesign_SE3nv**：用于运行 RFdiffusion 模型 和 ProteinMPNN 模型
- **segdesign_colabfold**：用于运行 AlphaFold2 模型（基于 ColabFold）

为方便用户安装环境，我们提供了安装脚本，位于 `environments/` 目录下。在运行脚本前，请提前安装 Conda 或 Miniconda，并且确保 conda 可以正常运行。

请运行安装脚本：

```bash
# 安装主环境（HMMER、MMSeqs2、DSSP 等）
bash ./environments/segdesign_env.sh

# 安装 SE3nv 环境（包含 RFdiffusion 和 ProteinMPNN）
bash ./environments/SE3nv_env.sh

# 安装 ESMFold 环境（需要 CUDA 支持）
bash ./environments/esmfold_env.sh

# 安装 ColabFold 环境（用于 AlphaFold2，需要 CUDA 支持）
bash ./environments/colabfold_env.sh
```

dssp 工具修复（更新，修复步骤已整合到 segdesign_env.sh 中，无需单独操作）：
dssp 4.5.5 存在一些问题，需按以下步骤处理：

```bash
#使用 conda env list 查询 segdesign 环境路径
conda env list | grep segdesign

target_dir="${env_path}/share/libcifpp"
archive="${target_dir}/components.cif.gz"
gunzip -f "$archive" 

```



#### 3. 安装数据库（可选）

进行 HMMER 分析时，可能需要下载序列数据库：

```bash
# 下载 UniRef90 数据库
bash environments/download_uniref90.sh

# 下载 UniRef100 数据库
bash environments/download_uniref100.sh
```

#### 4. 配置路径

您可以编辑 `config/setting.yaml` 文件，配置以下路径：
- Anaconda 安装路径
- RFdiffusion 安装路径
- ProteinMPNN 安装路径
一般情况下，您无需修改这些路径，使用默认值即可。

## 📋 配置文件说明

### 用户配置（`config/config.yaml`）

用户配置文件控制整个工作流程的参数：

```yaml
project:
  anaconda_path:                     # Anaconda 安装路径，不写则使用 conda run 命令
  protein_file: ./example/ggtdt_af3.cif  # 输入的蛋白质结构文件（支持 .pdb 和 .cif 格式）
  output_dir: ./ggtdt_example        # 输出目录
  chain: A                           # 待分析的链
  sequence_length: 386               # 完整序列长度
  segment: 8-22                      # 设计区域（可选，如果为None则只执行profile部分分析，不跑后续设计）

profile:
  database: ./database/ggtdt_af3_A_homologous_sequences.fasta  # 序列数据库
  bitscore: 0.3                      # HMMER bit score 阈值
  n_iter: 5                          # JackHMMER 迭代次数
  cpu: 10                            # CPU 核心数
  threshold: 0.6                     # 保守性阈值

rfdiffusion:
  helix: True                        # 按 α-螺旋设计
  strand:                            # 按 β-折叠设计
  num_designs: 10                    # 生成设计的数量
  threshold: 0.5                     # 设计质量阈值

mpnn:
  num_seq_per_target: 20             # 每个设计生成的序列数
  sampling_temp: 0.2                 # MPNN 采样温度
  batch_size: 5                      # 批处理大小
  seed: 42                           # 随机种子
  top_percent: 0.5                   # 顶部百分比选择

mmseqs:  # 可选模块
  min_seq_id: 0.8                    # 最小序列同一性
  s: 7.5                             # 灵敏度参数
  c: 1.0                             # 覆盖率阈值

esmfold:
  ptm_threshold: 0.54                # PTM 分数阈值
  plddt_threshold: 70                # pLDDT 分数阈值
  #ss: helix                         # 二级结构类型
  #ss_threshold: 0.3                 # 二级结构阈值

alphafold2:
  ptm_threshold: 0.54                # PTM 分数阈值
  plddt_threshold: 70                # pLDDT 分数阈值
  #ss: helix                         # 二级结构类型
  #ss_threshold: 0.3                 # 二级结构阈值
```

## 💻 使用方法

### 基本用法

example 目录提供了一个完整的配置文件示例，待设计蛋白为ggtdt_af3，其cif文件位于 `example/ggtdt_af3.cif`，您可以参考该文件来配置您的项目。
该示例使用的蛋白质数据库为 UniRef90，为节省同源搜索时间，我们已提前将待设计蛋白的同源序列整合为 fasta 文件，位于 `database/ggtdt_af3_A_homologous_sequences.fasta`。


完整运行一次 pipeline 示例（输出目录为 `ggtdt_example/`）：

```bash
python Segdesign.py --config config/config.yaml
```

### 示例：ggtdt 蛋白质设计

`example/ggtdt_example/` 目录包含完整的输出示例，该示例展示了如何使用 CIF 格式的输入文件进行蛋白质设计：

```bash
# 运行示例工作流程
python Segdesign.py --config example/ggtdt_example/config.yaml
```

该示例包含以下特点：
- 使用 CIF 格式的蛋白质结构文件（ggtdt_af3.cif）
- 自动转换为 PDB 格式进行处理
- 完整的 HMMER、RFdiffusion、MPNN、ESMFold 和 AlphaFold2 流程
- 包含详细的日志和输出文件

### 模块单独运行

您可以单独运行各个模块：

```bash
# 仅运行序列分析（HMMER）
conda run -n segdesign \
python -u ./Segdesign/hmmer/hmmer.py \
--input_file ./example/ggtdt_af3.cif \
--select_chain A \
--output_folder ./ggtdt_example/hmmer_out \
--bitscore 0.3 \
--n_iter 5 \
--database ./database/ggtdt_af3_A_homologous_sequences.fasta \
--cpu 10 \
--minimum_sequence_coverage 50 \
--minimum_column_coverage 70 \
--identity 0.3 \
--final_report_folder ./ggtdt_example   

# 仅运行蛋白质骨架设计（RFdiffusion）
conda run -n segdesign_SE3nv \
python -u ./Segdesign/rfdiffusion/rf_diffusion.py \
--run_inference_path ./RFdiffusion/scripts/run_inference.py \
--inference.input_pdb ./ggtdt_example/ggtdt_af3.pdb \
--inference.output_prefix ./ggtdt_example/rfdiffusion_out/sample/ggtdt_af3_A \
--inference.num_designs 10 \
--contigmap.contigs '[A1-386]' \
--contigmap.inpaint_str '[A8-22]' \
--diffuser.partial_T 50 \
--contigmap.inpaint_str_helix '[A8-22]' 
            
# 仅运行序列设计（MPNN）
conda run -n segdesign_SE3nv \
python -u ./Segdesign/mpnn/mpnn.py \
--parse_multiple_chains_path ./ProteinMPNN/helper_scripts/parse_multiple_chains.py \
--assign_fixed_chains_path ./ProteinMPNN/helper_scripts/assign_fixed_chains.py \
--make_fixed_positions_dict_path ./ProteinMPNN/helper_scripts/make_fixed_positions_dict.py \
--protein_mpnn_run_path ./ProteinMPNN/protein_mpnn_run.py \
--pdb_folder ./ggtdt_example/rfdiffusion_out/filter_results \
--output_folder ./ggtdt_example/mpnn_out \
--chain_list A \
--position_list A8-22 \
--num_seq_per_target 20 \
--sampling_temp 0.2 \
--seed 42 \
--batch_size 5 

# 仅运行序列聚类分析（MMseqs2）
conda run -n segdesign \
python -u ./Segdesign/mmseqs/mmseqs.py \
--input_folder ./ggtdt_example/mpnn_out/top_filter \
--output_folder ./ggtdt_example/mmseqs_out \
--position_list A8-22 \
--threads 8 \
--min_seq_id 0.8 \
--cov_mode 0 \
--coverage 1.0 \
--sensitivity 7.5   

# 仅运行结构预测（ESMFold）
conda run -n segdesign_esmfold \
python -u /home/wangxuming/SegDesign/Segdesign/esmfold/esmfold.py \
--input_folder ./ggtdt_example/mmseqs_out/results \
--output_folder ./ggtdt_example/esmfold_out \
--mmseqs_report_path ./ggtdt_example/mmseqs_report.csv 

# 仅运行结构预测（AlphaFold2）
conda run -n segdesign_colabfold \
python -u ./Segdesign/alphafold2/af2.py 
--input_file ./ggtdt_example/esmfold_out/filter_result.fa 
--output_folder ./ggtdt_example/alphafold2_out 
--esmfold_report_path ./ggtdt_example/esmfold_report.csv 
--amber True 
--templates True 
--gpu True 
--random_seed 0 
--num_recycle 3       
```

请注意，在单独运行模块时，由于RFdiffusion、MPNN、ESMFold和AlphaFold2的运行环境不在主环境（segdesign）中，以此不会生成相应的最终报告文件（rfdiffusion_report.csv、mpnn_report.csv、esmfold_report.csv、alphafold2_report.csv）。
如需要生成这些报告文件，那么在运行完相应的模块后，还需要在主环境下运行对应的脚本（如rfdiffusion_report.py、mpnn_report.py、esmfold_report.py、alphafold2_report.py）来生成最终的报告文件。

```bash
# 生成 RFdiffusion 最终验证报告
conda run -n segdesign \
python -u ./Segdesign/rfdiffusion/rf_diffusion_report.py \
--input_pdb ./ggtdt_example/ggtdt_af3.pdb \
--rfdiffusion_prefix ./ggtdt_example/rfdiffusion_out/sample/ggtdt_af3_A \
--inpaint_str '[A8-22]' \
--threshold 0.5 \
--final_report_folder ./ggtdt_example \
--ss helix  
               
# 生成 MPNN 最终验证报告
conda run -n segdesign \
python -u ./Segdesign/mpnn/mpnn_report.py \
--seq_folder ./ggtdt_example/mpnn_out/seqs \
--output_folder ./ggtdt_example/mpnn_out \
--top_percent 0.5 \
--generate_report True \
--final_report_folder ./ggtdt_example \
--position_list A8-22 \
--protein_pdb ./ggtdt_example/ggtdt_af3.pdb  

# 生成 ESMFold 最终验证报告
conda run -n segdesign \
python -u ./Segdesign/esmfold/esmfold_report.py \
--esmfold_folder ./ggtdt_example/esmfold_out \
--original_protein_chain_path ./ggtdt_example/hmmer_out/target_chain_pdb/ggtdt_af3_A.pdb \
--seq_range_str 8-22 \
--ss helix \
--ss_threshold 0.5 \
--ptm_threshold 0.54 \
--plddt_threshold 70 \

# 生成 AlphaFold2 最终验证报告
conda run -n segdesign \
python -u ./Segdesign/alphafold2/af2_report.py \
--esmfold_report_path ./ggtdt_example/esmfold_report.csv \
--alphafold2_folder ./ggtdt_example/alphafold2_out \
--seq_range_str 8-22 \
--ss helix \
--ss_threshold 0.5 \
--ptm_threshold 0.54 \
--plddt_threshold 70                  
```

由于通过命令行的方式运行单独模块非常麻烦，因此不推荐用户使用此方法。
若想单独运行模块，建议在 config.yaml 进行设置后，使用 Segdesign.py 读取配置并运行。

```bash
# 运行 Segdesign.py 
python -u ./Segdesign/Segdesign.py --config ./config/setting.yaml
```

## 📊 输出文件结构

```
output/
├── config.yaml                    # 配置文件的副本
├── hmmer_out/                     # HMMER 分析结果
│   ├── <protein_name>_<chain>_Recommended_Design_Area.txt
│   ├── <protein_name>_<chain>_conservative_comprehensive_report.csv
│   └── jackhmmer_out/            # 原始 HMMER 比对结果
├── rfdiffusion_out/              # RFdiffusion 结果
│   ├── sample/                   # 生成的骨架结构
│   └── filter_results/           # 过滤后的结构
├── mpnn_out/                     # MPNN 序列设计结果
│   ├── seqs/                     # 设计的序列
│   └── csv_files/                # 分析 CSV 文件
├── esmfold_out/                  # ESMFold 结构预测结果
│   ├── csv_files/                # 预测质量评估文件
│   ├── dssp_files/               # 二级结构分析文件
│   ├── filter_files/             # 过滤后的序列和结构
│   └── structure_prediction_files/  # 预测的 PDB 结构文件
├── alphafold2_out/               # AlphaFold2 结构预测结果
│   ├── colabfold_batch/           # ColabFold 批处理结果
│   ├── csv_files/                # 预测质量评估文件
│   ├── dssp_files/               # 二级结构分析文件
│   ├── filter_files/             # 过滤后的序列和结构
│   └── alphafold2_report.csv     # AlphaFold2 最终验证报告
├── esmfold_report.csv            # ESMFold 最终验证报告
└── alphafold2_report.csv         # AlphaFold2 最终验证报告
```

### 输出文件列说明

| 列名 | 说明 |
|------|------|
| index | 设计编号 |
| backbone | 骨架来源结构 |
| segment | 设计区域 |
| rfdiffusion_ss8 | RFdiffusion 二级结构（8类） |
| rfdiffusion_ss3 | RFdiffusion 二级结构（3类） |
| rfdiffusion_H_prop | RFdiffusion α-螺旋比例 |
| rfdiffusion_E_prop | RFdiffusion β-折叠比例 |
| rfdiffusion_C_prop | RFdiffusion 卷曲比例 |
| backbone_pdb | 骨架PDB文件路径 |
| score | 设计分数 |
| global_score | 全局分数 |
| region | 设计区域 |
| sequence | 设计序列 |
| esmfold_ptm | ESMFold PTM 分数 |
| esmfold_plddt | ESMFold pLDDT 置信度分数 |
| esmfold_ss8 | ESMFold 二级结构（8类） |
| esmfold_ss3 | ESMFold 二级结构（3类） |
| esmfold_H_prop | ESMFold α-螺旋比例 |
| esmfold_E_prop | ESMFold β-折叠比例 |
| esmfold_C_prop | ESMFold 卷曲比例 |
| af2_ptm | AlphaFold2 PTM 分数 |
| af2_plddt | AlphaFold2 pLDDT 置信度分数 |
| af2_ss8 | AlphaFold2 二级结构（8类） |
| af2_ss3 | AlphaFold2 二级结构（3类） |
| af2_H_prop | AlphaFold2 α-螺旋比例 |
| af2_E_prop | AlphaFold2 β-折叠比例 |
| af2_C_prop | AlphaFold2 卷曲比例 |
| ss_filter | 二级结构过滤结果 |
| whether_pass | 质量控制通过状态 |

## 🔧 模块详细说明

### 1. HMMER 模块
- 使用 JackHMMER 进行序列保守性分析
- 识别保守区域以智能选择设计区域
- 生成综合保守性报告

### 2. RFdiffusion 模块
- 为设计区域生成新的蛋白质骨架
- 支持二级结构约束（螺旋/折叠）
- 生成多个设计候选

### 3. ProteinMPNN 模块
- 为生成的骨架设计氨基酸序列
- 优化序列的稳定性和表达性
- 支持固定骨架位置

### 4. MMSeqs2 模块
- 进行序列聚类分析
- 过滤冗余序列
- 生成聚类报告

### 5. ESMFold 模块
- 使用深度学习预测验证设计结构
- 评估 pLDDT 和 PTM 分数
- 过滤低质量设计
- 使用 DSSP 进行二级结构分析

### 6. AlphaFold2 模块（基于 ColabFold）
- 使用 AlphaFold2 模型进行高精度结构预测
- 支持使用 PDB 模板提高预测准确性
- 可选 AMBER 力场优化
- 支持 GPU 加速
- 使用 DSSP 进行二级结构分析

## ⚠️ 常见问题处理

### GPU 内存不足
```bash
# 减小批量大小或设计数量
# 设置 GPU 内存限制环境变量
export CUDA_VISIBLE_DEVICES=0
```

### Conda 环境激活问题
```bash
# 确保 CONDA_PATH 设置正确
export CONDA_PATH="/path/to/anaconda3"
source $CONDA_PATH/etc/profile.d/conda.sh
```

### 数据库错误
- 验证 `config/setting.yaml` 中的数据库路径
- 确保数据库格式正确
- 检查文件权限

## 📧 联系方式

如有疑问或建议，请提交 issue 或联系作者。

---

<div align="center">

**祝您蛋白质设计愉快！🔬🧬**

</div>
