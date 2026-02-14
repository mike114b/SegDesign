# SegDesign: Intelligent Protein Segment Design Pipeline

<div align="center">

**An integrated pipeline for intelligent protein segment design combining sequence analysis, structure prediction, and generative modeling**

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[English](README_EN.md) | [中文](README_CN.md)

</div>

## 📖 Overview

SegDesign is an automated pipeline for intelligent protein segment design. It integrates multiple state-of-the-art bioinformatics tools and deep learning models to perform comprehensive protein analysis and design:

- **Sequence Conservation Analysis**: Using HMMER for evolutionary conservation analysis
- **Structure Generation**: Using RFdiffusion for targeted protein backbone generation
- **Sequence Design**: Using ProteinMPNN for amino acid sequence optimization
- **Structure Validation**: Using ESMFold and AlphaFold2 for predicted structure quality assessment
- **Secondary Structure Analysis**: Using DSSP for protein secondary structure analysis
- **Sequence Clustering**: Using MMSeqs2 for sequence similarity analysis

## 🏗️ Architecture

```
SegDesign/
├── Segdesign.py              # Main entry point
├── Segdesign/
│   ├── hmmer/               # Sequence conservation analysis
│   ├── rfdiffusion/         # Structure generation
│   ├── mpnn/                # Sequence design
│   ├── esmfold/             # Structure prediction
│   ├── alphafold2/          # AlphaFold2 structure prediction
│   ├── dssp/                # Secondary structure analysis
│   └── mmseqs/              # Sequence clustering
├── config/
│   ├── config.yaml          # User configuration
│   └── setting.yaml         # System settings
├── environments/            # Environment installation scripts
└── example/                 # Example outputs
```

## 🚀 Quick Start

### Prerequisites

- **Operating System**: Linux (recommended) or Windows with WSL2
- **Python**: 3.9 or higher
- **Conda/Miniconda**: Required for environment management
- **GPU**: NVIDIA GPU with CUDA support (strongly recommended for ESMFold and RFdiffusion)
- **Memory**: At least 16GB RAM (32GB+ recommended)
- **Storage**: At least 200GB free space

### Installation

#### 1. Clone the Repository

```bash
git clone https://github.com/mike114b/SegDesign2.git
cd SegDesign2
```

#### 2. Install Conda Environments

The project requires 4 conda environments to run different modules:
- **segdesign**: Main environment containing HMMER, MMSeqs2, DSSP, etc.
- **segdesign_esmfold**: For running ESMFold model
- **segdesign_SE3nv**: For running RFdiffusion model and ProteinMPNN model
- **segdesign_colabfold**: For running AlphaFold2 model (based on ColabFold)

We provide installation scripts in the `environments/` directory for user convenience. Before running the scripts, ensure that Conda or Miniconda is installed. You can use CONDA_PATH to specify the Conda installation path. If not specified, ensure conda runs properly and the script will use the conda run command to install environments.

Please run the installation scripts:

```bash
# You can set CONDA_PATH to specify Conda installation path, but this is optional
CONDA_PATH="/path/to/your/anaconda3"

# Install main environment (HMMER, MMSeqs2, DSSP, etc.)
bash ./environments/segdesign_env.sh

# Install SE3nv environment (containing RFdiffusion and ProteinMPNN)
bash ./environments/SE3nv_env.sh

# Install ESMFold environment (requires CUDA support)
bash ./environments/esmfold_env.sh

# Install ColabFold environment (for AlphaFold2, requires CUDA support)
bash ./environments/colabfold_env.sh
```

DSSP Tool Fix:
DSSP 4.5.5 has some issues, follow these steps to fix:

```bash
# Use conda env list to query segdesign environment path
conda env list | grep segdesign

target_dir="${env_path}/share/libcifpp"
archive="${target_dir}/components.cif.gz"
gunzip -f "$archive" 

```



#### 3. Install Databases (Optional)

For HMMER analysis, you may need to download sequence databases:

```bash
# Download UniRef90 database
bash environments/download_uniref90.sh

# Download UniRef100 database
bash environments/download_uniref100.sh
```

#### 4. Configure Paths

You can edit `config/setting.yaml` file to configure the following paths:
- Anaconda installation path
- RFdiffusion installation path
- ProteinMPNN installation path
In general, you don't need to modify these paths, using default values is fine.

## 📋 Configuration

### User Configuration (`config/config.yaml`)

The user configuration file controls the workflow parameters:

```yaml
project:
  anaconda_path:                     # Anaconda installation path, leave empty to use conda run command
  protein_file: ./example/ggtdt_af3.cif  # Input protein structure file (supports .pdb and .cif formats)
  output_dir: ./ggtdt_example        # Output directory
  chain: A                           # Chain to analyze
  sequence_length: 386               # Full sequence length
  segment: 8-22                      # Design region (optional, if None then only profile analysis is performed)

profile:
  database: ./database/ggtdt_af3_A_homologous_sequences.fasta  # Sequence database
  bitscore: 0.3                      # HMMER bit score threshold
  n_iter: 5                          # JackHMMER iterations
  cpu: 10                            # Number of CPU cores
  threshold: 0.6                     # Conservation threshold

rfdiffusion:
  helix: True                        # Design as alpha-helix
  strand:                            # Design as beta-strand
  num_designs: 10                    # Number of designs to generate
  threshold: 0.5                     # Design quality threshold

mpnn:
  num_seq_per_target: 20             # Sequences per design
  sampling_temp: 0.2                 # MPNN sampling temperature
  batch_size: 5                      # Batch size
  seed: 42                           # Random seed
  top_percent: 0.5                   # Top percentage selection

mmseqs:  # Optional module
  min_seq_id: 0.8                    # Minimum sequence identity
  s: 7.5                             # Sensitivity parameter
  c: 1.0                             # Coverage threshold

esmfold:
  ptm_threshold: 0.54                # PTM score threshold
  plddt_threshold: 70                # pLDDT score threshold
  #ss: helix                         # Secondary structure type
  #ss_threshold: 0.3                 # Secondary structure threshold

alphafold2:
  ptm_threshold: 0.54                # PTM score threshold
  plddt_threshold: 70                # pLDDT score threshold
  #ss: helix                         # Secondary structure type
  #ss_threshold: 0.3                 # Secondary structure threshold
```

## 💻 Usage

### Basic Usage

The example directory provides a complete configuration file example that you can reference to configure your project.
This example uses the protein database uniprot_sprot.fasta, which you can download using the download_uniprot_sprot.sh script.

```bash
bash example/download_uniprot_sprot.sh
```

Run the complete pipeline:

```bash
python Segdesign.py --config config/config.yaml
```

### Individual Module Execution

Individual modules can be run separately:

```bash
# Run sequence analysis only (HMMER)
conda run -n segdesign python ./Segdesign/hmmer/hmmer.py \
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

# Run protein backbone design only (RFdiffusion)
conda run -n segdesign_SE3nv python ./Segdesign/rfdiffusion/rf_diffusion.py \
--run_inference_path ./RFdiffusion/scripts/run_inference.py \
--inference.input_pdb ./ggtdt_example/ggtdt_af3.pdb \
--inference.output_prefix ./ggtdt_example/rfdiffusion_out/sample/ggtdt_af3_A \
--inference.num_designs 10 \
--contigmap.contigs '[A1-386]' \
--contigmap.inpaint_str '[A8-22]' \
--diffuser.partial_T 50 \
--contigmap.inpaint_str_helix '[A8-22]' 
            
# Run sequence design only (MPNN)
conda run -n segdesign_SE3nv python ./Segdesign/mpnn/mpnn.py \
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

# Run sequence clustering analysis only (MMSeqs2)
conda run -n segdesign python ./Segdesign/mmseqs/mmseqs.py \
--input_folder ./ggtdt_example/mpnn_out/top_filter \
--output_folder ./ggtdt_example/mmseqs_out \
--position_list A8-22 \
--threads 8 \
--min_seq_id 0.8 \
--cov_mode 0 \
--coverage 1.0 \
--sensitivity 7.5   

# Run structure prediction only (ESMFold)
conda run -n segdesign_esmfold python ./Segdesign/esmfold/esmfold.py \
--input_folder ./ggtdt_example/mmseqs_out/results \
--output_folder ./ggtdt_example/esmfold_out \
--mmseqs_report_path ./ggtdt_example/mmseqs_report.csv 

# Run structure prediction only (AlphaFold2)
conda run -n segdesign_colabfold python ./Segdesign/alphafold2/af2.py \
--input_file ./ggtdt_example/esmfold_out/filter_result.fa \
--output_folder ./ggtdt_example/alphafold2_out \
--esmfold_report_path ./ggtdt_example/esmfold_report.csv \
--amber True \
--templates True \
--gpu True \
--random_seed 0 \
--num_recycle 3       
```

**Note**: When running modules individually, since RFdiffusion, MPNN, ESMFold and AlphaFold2 run in environments other than the main environment (segdesign), the corresponding final report files (rfdiffusion_report.csv, mpnn_report.csv, esmfold_report.csv, alphafold2_report.csv) will not be generated. If you need to generate these report files, you need to run the corresponding scripts (such as rfdiffusion_report.py, mpnn_report.py, esmfold_report.py, alphafold2_report.py) in the main environment after running the corresponding modules.

```bash
# Generate RFdiffusion final validation report
conda run -n segdesign python ./Segdesign/rfdiffusion/rf_diffusion_report.py \
--input_pdb ./ggtdt_example/ggtdt_af3.pdb \
--rfdiffusion_prefix ./ggtdt_example/rfdiffusion_out/sample/ggtdt_af3_A \
--inpaint_str '[A8-22]' \
--threshold 0.5 \
--final_report_folder ./ggtdt_example \
--ss helix  
               
# Generate MPNN final validation report
conda run -n segdesign python ./Segdesign/mpnn/mpnn_report.py \
--seq_folder ./ggtdt_example/mpnn_out/seqs \
--output_folder ./ggtdt_example/mpnn_out \
--top_percent 0.5 \
--generate_report True \
--final_report_folder ./ggtdt_example \
--position_list A8-22 \
--protein_pdb ./ggtdt_example/ggtdt_af3.pdb  

# Generate ESMFold final validation report
conda run -n segdesign python ./Segdesign/esmfold/esmfold_report.py \
--esmfold_folder ./ggtdt_example/esmfold_out \
--original_protein_chain_path ./ggtdt_example/hmmer_out/target_chain_pdb/ggtdt_af3_A.pdb \
--seq_range_str 8-22 \
--ss helix \
--ss_threshold 0.5 \
--ptm_threshold 0.54 \
--plddt_threshold 70 \

# Generate AlphaFold2 final validation report
conda run -n segdesign python ./Segdesign/alphafold2/af2_report.py \
--esmfold_report_path ./ggtdt_example/esmfold_report.csv \
--alphafold2_folder ./ggtdt_example/alphafold2_out \
--seq_range_str 8-22 \
--ss helix \
--ss_threshold 0.5 \
--ptm_threshold 0.54 \
--plddt_threshold 70                  
```

Since running individual modules via command line is very cumbersome, this method is not recommended. If you want to run individual modules, it is recommended to set them in config.yaml and use Segdesign.py to read the configuration and run.

```bash
# Run Segdesign.py
python Segdesign.py --config ./config/config.yaml
```

## 📊 Output Structure

```
output/
├── config.yaml                    # Copy of configuration
├── hmmer_out/                     # HMMER analysis results
│   ├── <protein_name>_<chain>_Recommended_Design_Area.txt
│   ├── <protein_name>_<chain>_conservative_comprehensive_report.csv
│   └── jackhmmer_out/            # Raw HMMER alignments
├── rfdiffusion_out/              # RFdiffusion results
│   ├── sample/                   # Generated backbones
│   └── filter_results/           # Filtered structures
├── mpnn_out/                     # MPNN sequence designs
│   ├── seqs/                     # Designed sequences
│   └── csv_files/                # Analysis CSVs
├── esmfold_out/                  # ESMFold structure prediction results
│   ├── csv_files/                # Prediction quality assessment files
│   ├── dssp_files/               # Secondary structure analysis files
│   ├── filter_files/             # Filtered sequences and structures
│   └── structure_prediction_files/  # Predicted PDB structure files
├── alphafold2_out/               # AlphaFold2 structure prediction results
│   ├── colabfold_batch/           # ColabFold batch processing results
│   ├── csv_files/                # Prediction quality assessment files
│   ├── dssp_files/               # Secondary structure analysis files
│   ├── filter_files/             # Filtered sequences and structures
│   └── alphafold2_report.csv     # AlphaFold2 final validation report
├── esmfold_report.csv            # ESMFold final validation report
└── alphafold2_report.csv         # AlphaFold2 final validation report
```

### Output Columns Description

| Column | Description |
|--------|-------------|
| index | Design identifier |
| backbone | Source backbone structure |
| segment | Designed region |
| rfdiffusion_ss8 | RFdiffusion secondary structure (8 classes) |
| rfdiffusion_ss3 | RFdiffusion secondary structure (3 classes) |
| rfdiffusion_H_prop | RFdiffusion α-helix proportion |
| rfdiffusion_E_prop | RFdiffusion β-sheet proportion |
| rfdiffusion_C_prop | RFdiffusion coil proportion |
| backbone_pdb | Backbone PDB file path |
| score | Design score |
| global_score | Global score |
| region | Design region |
| sequence | Designed sequence |
| esmfold_ptm | ESMFold PTM score |
| esmfold_plddt | ESMFold pLDDT confidence score |
| esmfold_ss8 | ESMFold secondary structure (8 classes) |
| esmfold_ss3 | ESMFold secondary structure (3 classes) |
| esmfold_H_prop | ESMFold α-helix proportion |
| esmfold_E_prop | ESMFold β-sheet proportion |
| esmfold_C_prop | ESMFold coil proportion |
| af2_ptm | AlphaFold2 PTM score |
| af2_plddt | AlphaFold2 pLDDT confidence score |
| af2_ss8 | AlphaFold2 secondary structure (8 classes) |
| af2_ss3 | AlphaFold2 secondary structure (3 classes) |
| af2_H_prop | AlphaFold2 α-helix proportion |
| af2_E_prop | AlphaFold2 β-sheet proportion |
| af2_C_prop | AlphaFold2 coil proportion |
| ss_filter | Secondary structure filter result |
| whether_pass | Quality control pass status |

## 🔧 Module Details

### 1. HMMER Module
- Performs sequence conservation analysis using JackHMMER
- Identifies conserved regions for intelligent design area selection
- Generates comprehensive conservation reports

### 2. RFdiffusion Module
- Generates novel protein backbones for the design region
- Supports secondary structure constraints (helix/strand)
- Produces multiple design candidates

### 3. ProteinMPNN Module
- Designs amino acid sequences for generated backbones
- Optimizes sequences for stability and expression
- Supports fixed backbone positions

### 4. MMSeqs2 Module (Optional)
- Performs sequence clustering analysis
- Filters redundant sequences
- Generates cluster reports

### 5. ESMFold Module
- Validates designed structures using deep learning prediction
- Assesses pLDDT and PTM scores
- Filters low-quality designs
- Performs secondary structure analysis using DSSP

### 6. AlphaFold2 Module (Based on ColabFold)
- Performs high-accuracy structure prediction using AlphaFold2 model
- Supports using PDB templates to improve prediction accuracy
- Optional AMBER force field optimization
- Supports GPU acceleration
- Performs secondary structure analysis using DSSP

## ⚠️ Troubleshooting

### GPU Memory Issues
```bash
# Reduce batch size or number of designs
# Set environment variable for GPU memory limit
export CUDA_VISIBLE_DEVICES=0
```

### Conda Environment Activation
```bash
# Ensure CONDA_PATH is set correctly
export CONDA_PATH="/path/to/anaconda3"
source $CONDA_PATH/etc/profile.d/conda.sh
```

### Database Errors
- Verify database paths in `config/setting.yaml`
- Ensure databases are properly formatted
- Check file permissions

## 📧 Contact

For questions or suggestions, please open an issue or contact the author.

---

<div align="center">

**Happy Protein Designing! 🔬🧬**

</div>
