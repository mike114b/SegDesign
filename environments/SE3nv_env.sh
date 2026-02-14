#!/bin/bash

set -e  # 出错立即退出，避免静默失败

echo "=== 安装ProteinMPNN库 ==="
ProteinMPNN_DIR="./ProteinMPNN"
if [ ! -d "${ProteinMPNN_DIR}" ]; then
    # 不存在则创建
   echo "=== 克隆ProteinMPNN库到本地 ==="
   echo "git clone https://github.com/dauparas/ProteinMPNN.git"
   git clone https://github.com/dauparas/ProteinMPNN.git

else
    echo "ProteinMPNN库已存在，跳过"
fi

echo "=== 安装RFdiffusion库 ==="
RFdiffusion_DIR="./RFdiffusion"
if [ ! -d "${RFdiffusion_DIR}" ]; then
    # 不存在则创建
   echo "克隆RFdiffusion库到本地"
   echo "git clone https://github.com/RosettaCommons/RFdiffusion.git"
   git clone https://github.com/RosettaCommons/RFdiffusion.git
else
    echo "RFdiffusion库已存在，跳过"
fi

echo "将模型权重下载到 RFDiffusion/models 目录"
echo "cd RFdiffusion"
cd RFdiffusion
echo "mkdir models && cd models"
if [ ! -d "models" ]; then
  mkdir models
fi

cd models
echo "下载模型权重"
FILES=("Base_ckpt.pt" "Complex_base_ckpt.pt" "Complex_Fold_base_ckpt.pt"
"InpaintSeq_ckpt.pt" "InpaintSeq_Fold_ckpt.pt" "ActiveSite_ckpt.pt"
"Base_epoch8_ckpt.pt" "Complex_beta_ckpt.pt" "RF_structure_prediction_weights.pt")

ADDRESS=(
"http://files.ipd.uw.edu/pub/RFdiffusion/6f5902ac237024bdd0c176cb93063dc4/Base_ckpt.pt"
"http://files.ipd.uw.edu/pub/RFdiffusion/e29311f6f1bf1af907f9ef9f44b8328b/Complex_base_ckpt.pt"
"http://files.ipd.uw.edu/pub/RFdiffusion/60f09a193fb5e5ccdc4980417708dbab/Complex_Fold_base_ckpt.pt"
"http://files.ipd.uw.edu/pub/RFdiffusion/74f51cfb8b440f50d70878e05361d8f0/InpaintSeq_ckpt.pt"
"http://files.ipd.uw.edu/pub/RFdiffusion/76d00716416567174cdb7ca96e208296/InpaintSeq_Fold_ckpt.pt"
"http://files.ipd.uw.edu/pub/RFdiffusion/5532d2e1f3a4738decd58b19d633b3c3/ActiveSite_ckpt.pt"
"http://files.ipd.uw.edu/pub/RFdiffusion/12fc204edeae5b57713c5ad7dcb97d39/Base_epoch8_ckpt.pt"
"http://files.ipd.uw.edu/pub/RFdiffusion/f572d396fae9206628714fb2ce00f72e/Complex_beta_ckpt.pt"
"http://files.ipd.uw.edu/pub/RFdiffusion/1befcb9b28e2f778f53d47f18b7597fa/RF_structure_prediction_weights.pt"
)

DOWNLOAD_DIR="."

for ((i=0; i<${#FILES[@]}; i++)); do
    # 获取当前索引对应的文件名和下载地址
    FILE_NAME="${FILES[$i]}"
    FILE_PATH="${DOWNLOAD_DIR}/${FILE_NAME}"
    DOWNLOAD_URL="${ADDRESS[$i]}"
    # 拼接文件的完整路径
    echo "检查${FILE_NAME}是否存在"
    # 3. 检查文件是否存在
    #echo "${FILE_PATH}"
    #echo "[ -f "${FILE_PATH}" ]"
    if [ -f "${FILE_PATH}" ]; then
        echo -e "[SKIP] 文件 ${FILE_NAME} 已存在于 ${DOWNLOAD_DIR}，跳过下载\n"
        continue  # 跳过当前文件，继续下一个
    else
      echo -e "[DOWNLOAD] 开始下载 ${FILE_NAME} ..."
      # wget 参数说明：
      # -O：指定保存的文件路径
      # --continue：断点续传（网络中断后重新下载不会从头开始）
      # -q：静默模式（可选，去掉则显示下载进度）
      # -t 3：重试3次（避免网络波动导致下载失败）
      echo " wget "${DOWNLOAD_URL}""
      wget "${DOWNLOAD_URL}"
    fi
    # 4. 文件不存在，执行下载


done


echo "返回上一级目录"
echo "cd .."
cd ..




ENV_NAME="segdesign_SE3nv"
echo "=== RFdiffusion模块环境的安装 ==="
echo "创建环境..."
echo "conda create -n '$ENV_NAME' python=3.9 -y"
conda create -n "$ENV_NAME" python=3.9 -y


#写入你anaconda的安装路径
#CONDA_PATH="/opt/software/anaconda3"

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
    echo "conda activate '$ENV_NAME'"
    conda activate "$ENV_NAME"

    echo "安装cudatoolkit 11.8"
    echo "conda install -y cudatoolkit=11.8 -c conda-forge -y"
    conda install -y cudatoolkit=11.8 -c conda-forge -y

    echo "安装pytorch 2.1.0"
    echo "pip install torch==2.1.0 torchvision==0.16.0 torchaudio==2.1.0 --index-url https://download.pytorch.org/whl/cu118"
    pip install torch==2.1.0 torchvision==0.16.0 torchaudio==2.1.0 --index-url https://download.pytorch.org/whl/cu118

    echo "pip install  dgl -f https://data.dgl.ai/wheels/torch-2.1/cu118/repo.html"
    pip install  dgl -f https://data.dgl.ai/wheels/torch-2.1/cu118/repo.html

    echo "pip install hydra-core pyrsistent"
    pip install hydra-core pyrsistent

    echo "安装 SE3-Transformer"
    echo "cd env/SE3Transformer"
    cd env/SE3Transformer
    echo "安装环境依赖项..."
    echo "pip install \
    e3nn==0.3.3 \
    wandb==0.12.0 \
    pynvml==11.0.0 \
    decorator==5.1.0 \
    --no-cache-dir "
    pip install \
    e3nn==0.3.3 \
    wandb==0.12.0 \
    pynvml==11.0.0 \
    decorator==5.1.0 \
    --no-cache-dir
    echo "pip install git+https://github.com/NVIDIA/dllogger#egg=dllogger --no-cache-dir "
    pip install git+https://github.com/NVIDIA/dllogger#egg=dllogger --no-cache-dir

    echo "python setup.py install"
    python setup.py install
    echo "切换到仓库的根目录"
    echo "cd ../.. "
    cd ../..
    echo "从代码库的根目录安装 rfdiffusion 模块"
    echo "pip install -e ."
    pip install -e .

    echo "安装pandas库"
    echo "conda install pandas numpy=1.26.4 -y"
    conda install pandas numpy=1.26.4 -y
else
    echo "未检测到 CONDA_PATH，使用 conda run 方式安装"
    echo "安装cudatoolkit 11.8"
    echo "conda run -n '$ENV_NAME' conda install -y cudatoolkit=11.8 -c conda-forge -y"
    conda run -n "$ENV_NAME" conda install -y cudatoolkit=11.8 -c conda-forge -y

    echo "安装pytorch 2.1.0"
    echo "conda run -n '$ENV_NAME' pip install torch==2.1.0 torchvision==0.16.0 torchaudio==2.1.0 --index-url https://download.pytorch.org/whl/cu118"
    conda run -n "$ENV_NAME" pip install torch==2.1.0 torchvision==0.16.0 torchaudio==2.1.0 --index-url https://download.pytorch.org/whl/cu118

    echo "conda run -n '$ENV_NAME' pip install  dgl -f https://data.dgl.ai/wheels/torch-2.1/cu118/repo.html"
    conda run -n "$ENV_NAME" pip install  dgl -f https://data.dgl.ai/wheels/torch-2.1/cu118/repo.html

    echo "conda run -n '$ENV_NAME' pip install hydra-core pyrsistent"
    conda run -n "$ENV_NAME" pip install hydra-core pyrsistent

    echo "安装 SE3-Transformer"
    echo "cd env/SE3Transformer"
    cd env/SE3Transformer
    echo "安装环境依赖项..."
    echo "conda run -n '$ENV_NAME' pip install \
    e3nn==0.3.3 \
    wandb==0.12.0 \
    pynvml==11.0.0 \
    decorator==5.1.0 \
    --no-cache-dir "
    conda run -n "$ENV_NAME" pip install \
    e3nn==0.3.3 \
    wandb==0.12.0 \
    pynvml==11.0.0 \
    decorator==5.1.0 \
    --no-cache-dir
    echo "conda run -n '$ENV_NAME' pip install git+https://github.com/NVIDIA/dllogger#egg=dllogger --no-cache-dir "
    conda run -n "$ENV_NAME" pip install git+https://github.com/NVIDIA/dllogger#egg=dllogger --no-cache-dir

    echo "conda run -n '$ENV_NAME' python setup.py install"
    conda run -n "$ENV_NAME" python setup.py install
    echo "切换到仓库的根目录"
    echo "cd ../.. "
    cd ../..
    echo "从代码库的根目录安装 rfdiffusion 模块"
    echo "conda run -n '$ENV_NAME' pip install -e ."
    conda run -n "$ENV_NAME" pip install -e .

    echo "安装pandas库"
    echo "conda run -n '$ENV_NAME' conda install pandas numpy=1.26.4 -y"
    conda run -n "$ENV_NAME" conda install pandas numpy=1.26.4 -y
fi