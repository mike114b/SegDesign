#!/bin/bash

set -e  # 出错立即退出，避免静默失败
# ESMFold CUDA 12.4 分步安装脚本

echo "=== ESMFold 安装脚本 for CUDA 12.4 ==="

ENV_NAME="segdesign_esmfold"
# 1. 创建基础环境
echo "创建基础环境..."
echo "conda create -n '$ENV_NAME' python=3.9 -y"
conda create -n "$ENV_NAME" python=3.9 -y

# ===================== 获取 Anaconda 安装路径 =====================
# 尝试直接执行 conda info --base 获取根路径
CONDA_PATH=$(conda info --base 2>/dev/null)

# 异常处理1：conda 命令未找到（未初始化 conda）
if [ -z "$CONDA_PATH" ]; then
    echo "提示：未检测到 conda 命令，尝试初始化 conda..."
    # 尝试常见的 Anaconda 安装路径（适配大多数用户）
    common_CONDA_PATHs=(
        "$HOME/anaconda3"
        "$HOME/miniconda3"
        "$HOME/conda"
        "/opt/anaconda3"
        "/usr/local/anaconda3"
    )
    # 遍历常见路径，检查是否存在 conda 可执行文件
    for path in "${common_CONDA_PATHs[@]}"; do
        if [ -f "$path/bin/conda" ]; then
            # 初始化 conda 并重新获取路径
            source "$path/bin/activate" >/dev/null 2>&1
            CONDA_PATH=$(conda info --base 2>/dev/null)
            break
        fi
    done
fi

# 异常处理2：仍未获取到路径（Anaconda 未安装）
if [ -z "$CONDA_PATH" ] || [ ! -d "$CONDA_PATH" ]; then
    echo "错误：未找到 Anaconda/Miniconda 安装路径！"
    exit 1
fi


# ===================== 验证并输出结果 =====================
echo "✅ Anaconda 安装路径已获取："
echo "CONDA_PATH = $CONDA_PATH"
#--------------------------------------------------

if [ -n "$CONDA_PATH" ]; then
    echo "检测到 CONDA_PATH，使用环境激活方式安装"
    #写入你anaconda的安装路径
    #CONDA_PATH="/opt/software/anaconda3"
    # 加载conda环境
    if [ -f "$CONDA_PATH/etc/profile.d/conda.sh" ]; then
        source "$CONDA_PATH/etc/profile.d/conda.sh"
    elif [ -f "$CONDA_PATH/bin/activate" ]; then
        source "$CONDA_PATH/bin/activate"
    else
        echo "找不到conda激活脚本" >&2
        exit 1
    fi

    echo "进入虚拟环境..."
    echo "conda activate esmfold"
    conda activate esmfold

    #2.安装前置
    echo "安装前置"
    echo "pip install setuptools==59.5.0 --no-cache-dir"
    pip install setuptools==59.5.0 --no-cache-dir
    echo "conda install conda-forge::openmm=7.5.1 -y"
    conda install conda-forge::openmm=7.5.1 -y
    echo "conda install -c conda-forge pdbfixer einops fairscale omegaconf=2.1.1 hydra-core=1.1.0 pytest -y "
    conda install -c conda-forge pdbfixer einops fairscale omegaconf=2.1.1 hydra-core=1.1.0 pytest -y

    # 3. 安装生物信息学工具
    echo "安装生物信息学工具..."
    echo "conda install -c conda-forge -c bioconda hhsuite=3.3.0 -y"
    conda install -c conda-forge -c bioconda hhsuite=3.3.0 -y
    echo "conda install bioconda::hmmer==3.3.2 -y"
    conda install bioconda::hmmer==3.3.2 -y
    echo "conda install bioconda::kalign2==2.04 -y"
    conda install bioconda::kalign2==2.04 -y

    # 4. 安装 PyTorch 1.12 (CUDA 11.8版本，但兼容 CUDA 12.4)
    echo "安装 PyTorch 1.12..."
    echo "pip install torch==2.4.0 torchvision==0.19.0 torchaudio==2.4.0 --index-url https://download.pytorch.org/whl/cu124"
    pip install torch==2.4.0 torchvision==0.19.0 torchaudio==2.4.0 --index-url https://download.pytorch.org/whl/cu124

    echo "降级pip为24.0版本"
    echo "install pip==24.0"
    pip install pip==24.0

    # 5. 安装 pip 包
    echo "安装 pip 包..."
    echo "pip install \
    biopython==1.79 \
    dm-tree==0.1.6 \
    ml-collections==0.1.0 \
    PyYAML==5.4.1 \
    requests==2.26.0 \
    scipy  \
    tqdm==4.62.2 \
    typing-extensions==4.8.0 \
    pytorch_lightning==1.6.5 \
    wandb==0.12.21 \
    pandas==1.4.4 \
    modelcif \
    "
    pip install \
    biopython==1.79 \
    dm-tree==0.1.6 \
    ml-collections==0.1.0 \
    PyYAML==5.4.1 \
    requests==2.26.0 \
    scipy  \
    tqdm==4.62.2 \
    typing-extensions==4.8.0 \
    pytorch_lightning==1.6.5 \
    wandb==0.12.21 \
    pandas==1.4.4 \
    modelcif \


    #7.安装esmfold
    echo "安装esmfold"
    echo "pip install 'fair-esm[esmfold]'"
    pip install "fair-esm[esmfold]"

    # 8. 安装 NVIDIA dllogger
    echo "安装 NVIDIA dllogger..."
    echo "pip install git+https://github.com/NVIDIA/dllogger.git"
    pip install git+https://github.com/NVIDIA/dllogger.git

    #9.安装openfold
    echo "安装openfold"
    echo "cd openfold"
    cd openfold
    echo "pip install --no-build-isolation ."
    pip install . --no-build-isolation
    echo "pip uninstall deepspeed -y"
    pip uninstall deepspeed -y

    pip install numpy==1.22.4
else
    echo "未检测到 CONDA_PATH，使用 conda run 方式安装"

    #2.安装前置
    echo "安装前置"
    echo "conda run -n '$ENV_NAME' pip install setuptools==59.5.0 --no-cache-dir"
    conda run -n "$ENV_NAME" pip install setuptools==59.5.0 --no-cache-dir
    echo "conda run -n '$ENV_NAME' conda install conda-forge::openmm=7.5.1 -y"
    conda run -n "$ENV_NAME" conda install conda-forge::openmm=7.5.1 -y
    echo "conda run -n '$ENV_NAME' conda install -c conda-forge pdbfixer einops fairscale omegaconf=2.1.1 hydra-core=1.1.0 pytest -y "
    conda run -n "$ENV_NAME" conda install -c conda-forge pdbfixer einops fairscale omegaconf=2.1.1 hydra-core=1.1.0 pytest -y

    # 3. 安装生物信息学工具
    echo "安装生物信息学工具..."
    echo "conda run -n '$ENV_NAME' conda install -c conda-forge -c bioconda hhsuite=3.3.0 -y"
    conda run -n "$ENV_NAME" conda install -c conda-forge -c bioconda hhsuite=3.3.0 -y
    echo "conda run -n '$ENV_NAME' conda install bioconda::hmmer==3.3.2 -y"
    conda run -n "$ENV_NAME" conda install bioconda::hmmer==3.3.2 -y
    echo "conda run -n '$ENV_NAME' conda install bioconda::kalign2==2.04 -y"
    conda run -n "$ENV_NAME" conda install bioconda::kalign2==2.04 -y

    # 4. 安装 PyTorch 1.12 (CUDA 11.8版本，但兼容 CUDA 12.4)
    echo "安装 PyTorch 1.12..."
    echo "conda run -n '$ENV_NAME' pip install torch==2.4.0 torchvision==0.19.0 torchaudio==2.4.0 --index-url https://download.pytorch.org/whl/cu124"
    conda run -n "$ENV_NAME" pip install torch==2.4.0 torchvision==0.19.0 torchaudio==2.4.0 --index-url https://download.pytorch.org/whl/cu124

    echo "降级pip为24.0版本"
    echo "conda run -n '$ENV_NAME' pip install pip==24.0"
    conda run -n "$ENV_NAME" pip install pip==24.0

    # 5. 安装 pip 包
    echo "安装 pip 包..."
    echo "conda run -n '$ENV_NAME' pip install \
    biopython==1.79 \
    dm-tree==0.1.6 \
    ml-collections==0.1.0 \
    PyYAML==5.4.1 \
    requests==2.26.0 \
    scipy  \
    tqdm==4.62.2 \
    typing-extensions==4.8.0 \
    pytorch_lightning==1.6.5 \
    wandb==0.12.21 \
    pandas==1.4.4 \
    modelcif \
    "
    conda run -n "$ENV_NAME" pip install \
    biopython==1.79 \
    dm-tree==0.1.6 \
    ml-collections==0.1.0 \
    PyYAML==5.4.1 \
    requests==2.26.0 \
    scipy  \
    tqdm==4.62.2 \
    typing-extensions==4.8.0 \
    pytorch_lightning==1.6.5 \
    wandb==0.12.21 \
    pandas==1.4.4 \
    modelcif \


    #7.安装esmfold
    echo "安装esmfold"
    echo "conda run -n '$ENV_NAME' pip install 'fair-esm[esmfold]'"
    conda run -n "$ENV_NAME" pip install "fair-esm[esmfold]"

    # 8. 安装 NVIDIA dllogger
    echo "安装 NVIDIA dllogger..."
    echo "conda run -n '$ENV_NAME' pip install git+https://github.com/NVIDIA/dllogger.git"
    conda run -n "$ENV_NAME" pip install git+https://github.com/NVIDIA/dllogger.git

    #9.安装openfold
    echo "安装openfold"
    echo "cd openfold"
    cd openfold
    echo "conda run -n '$ENV_NAME' pip install --no-build-isolation ."
    conda run -n "$ENV_NAME" pip install . --no-build-isolation
    echo "conda run -n '$ENV_NAME' pip uninstall deepspeed -y"
    conda run -n "$ENV_NAME" pip uninstall deepspeed -y

    conda run -n "$ENV_NAME" pip install numpy==1.22.4

fi