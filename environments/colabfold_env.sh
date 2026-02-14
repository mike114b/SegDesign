#!/bin/bash -e

type wget 2>/dev/null || { echo "wget is not installed. Please install it using apt or yum." ; exit 1 ; }

CURRENTPATH=`pwd`
COLABFOLDDIR="${CURRENTPATH}/localcolabfold"
mkdir -p "${COLABFOLDDIR}"
cd "${COLABFOLDDIR}"

ENV_NAME="segdesign_colabfold"

# 1. 创建基础环境
echo "创建基础环境..."
echo "conda create -n '$ENV_NAME' python=3.9 -y"
conda create -n "$ENV_NAME" -c conda-forge -c bioconda \
python=3.10 \
openmm==8.2.0 pdbfixer kalign2=2.04 hhsuite=3.3.0 mmseqs2 -y

ENVPATH=$(conda run -n "$ENV_NAME" python -c "import sys; print(sys.prefix)")

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
    echo "conda activate '$ENV_NAME'"
    conda activate "$ENV_NAME"

    echo "安装 ColabFold 但不包含 JAX（后面单独安装）"
    echo 'pip install --no-warn-conflicts "colabfold[alphafold-minus-jax] @ git+https://github.com/sokrypton/ColabFold"'
    pip install --no-warn-conflicts "colabfold[alphafold-minus-jax] @ git+https://github.com/sokrypton/ColabFold"

    echo "安装完整的 AlphaFold 组件"
    echo 'pip install "colabfold[alphafold]"'
    pip install "colabfold[alphafold]"

    echo "安装 ColabFold 的核心计算引擎 jax"
    echo 'pip install --upgrade "jax[cuda12]==0.5.3"'
    pip install --upgrade "jax[cuda12]==0.5.3"

    echo '安装 TensorFlow（某些组件依赖）'
    echo 'pip install --upgrade tensorflow'
    pip install --upgrade tensorflow

    echo '安装 TensorFlow 警告抑制库'
    echo 'pip install silence_tensorflow'
    pip install silence_tensorflow

    #wget -qnc -O "$COLABFOLDDIR/update_linux.sh" \
    #https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/update_linux.sh

    pushd "${ENVPATH}/lib/python3.10/site-packages/colabfold"
    # Use 'Agg' for non-GUI backend
    sed -i -e "s#from matplotlib import pyplot as plt#import matplotlib\nmatplotlib.use('Agg')\nimport matplotlib.pyplot as plt#g" plot.py
    # modify the default params directory
    sed -i -e "s#appdirs.user_cache_dir(__package__ or \"colabfold\")#\"${COLABFOLDDIR}/colabfold\"#g" download.py
    # suppress warnings related to tensorflow
    sed -i -e "s#from io import StringIO#from io import StringIO\nfrom silence_tensorflow import silence_tensorflow\nsilence_tensorflow()#g" batch.py
    # remove cache directory
    rm -rf __pycache__
    popd

    "$ENVPATH/bin/python3" -m colabfold.download
    echo "Download of alphafold2 weights finished."
    echo "-----------------------------------------"
    echo "Installation of ColabFold finished."


















