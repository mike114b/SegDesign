import shutil
import subprocess
import os
import logging
from typing import Dict, Optional, List
import shlex
import argparse
from pathlib import Path
import yaml
import sys
import threading
from Bio.PDB import MMCIFParser, PDBIO
from numpy.ma.core import identity


# 配置项（可根据实际情况修改）
CONFIG = {
    "MODULES":{
        'hmmer': {"path":'./Segdesign/hmmer/hmmer.py'},
        'rfdiffusion': {"path":'./Segdesign/rfdiffusion/rf_diffusion.py'},
        'rfdiffusion_report': {"path":'./Segdesign/rfdiffusion/rf_diffusion_report.py'},
        'mpnn': {"path":'./Segdesign/mpnn/mpnn.py'},
        'mmseqs': {"path":'./Segdesign/mmseqs/mmseqs.py'},
        'mpnn_report': {"path":'./Segdesign/mpnn/mpnn_report.py'},
        'esmfold': {"path":'./Segdesign/esmfold/esmfold.py'},
        'esmfold_report': {"path":'./Segdesign/esmfold/esmfold_report.py'},
        'alphafold2': {"path":'./Segdesign/alphafold2/af2.py'},
        'alphafold2_report': {"path":'./Segdesign/alphafold2/af2_report.py'},
        'dssp': {"path":'./dssp/dssp.py'},
        'cluster_analysis':{"path":'./Segdesign/mpnn/cluster_analysis.py'},
    },
    "CONFIG_PATH": {
        "MAIN": "./config/config.yaml",
        "SETTING": "./config/setting.yaml"
    }
}


def setup_logger(log_path, console_output=True):
    # 1. 处理日志目录（创建不存在的目录）
    log_dir = os.path.dirname(log_path)
    if log_dir and not os.path.exists(log_dir):
        os.makedirs(log_dir, exist_ok=True)

    # 2. 创建/获取命名日志器（__name__）
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)  # 设置日志器级别（必须≥输出日志的级别）
    logger.propagate = False  # 禁用传播（避免重复输出）

    # 3. 清空已有Handler（避免多次调用重复输出）
    if logger.handlers:
        logger.handlers.clear()

    # 4. 给命名日志器绑定文件Handler（核心：让日志有输出目标）
    file_handler = logging.FileHandler(
        log_path,
        mode='w',  # 覆盖模式，清空旧内容
        encoding='utf-8'  # 避免中文乱码
    )
    # 配置日志格式
    formatter = logging.Formatter(
        "%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # 可选：添加控制台Handler（日志同时输出到控制台+文件）
    if console_output:
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)


    return logger



def validate_environment(env_name: str) -> bool:
    """验证Conda环境是否存在"""
    conda_info_cmd = [
        f"{CONFIG['MINICONDA_PATH']}/bin/conda",
        "info",
        "--envs"
    ]

    try:
        result = subprocess.run(
            conda_info_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
            timeout=30
        )
        # 检查环境是否在输出中（支持完整名称匹配）
        return any(f"*{env_name}" in line or f"  {env_name} " in line for line in result.stdout.splitlines())
    except subprocess.TimeoutExpired:
        logger.warning(f"验证环境 {env_name} 超时")
        logger.warning(f"Environment {env_name} validation timeout")
        return False
    except subprocess.CalledProcessError as e:
        logger.error(f"验证环境失败: {e.stderr}")
        logger.error(f"Environment validation failed: {e.stderr}")
        return False


def validate_module(module_name: str) -> str:
    """验证模块是否存在并返回完整路径"""
    if module_name not in CONFIG['MODULES']:
        raise ModuleRunnerError(f"模块 {module_name} 未在配置中定义，可用模块: {list(CONFIG['MODULES'].keys())}")
        raise ModuleRunnerError(f"Module {module_name} is not defined in configuration, available modules: {list(CONFIG['MODULES'].keys())}")

    module_path = os.path.abspath(CONFIG['MODULES'][module_name]['path'])
    if not os.path.exists(module_path):
        raise ModuleRunnerError(f"模块文件不存在: {module_path}")
        raise ModuleRunnerError(f"Module file does not exist: {module_path}")

    if not os.access(module_path, os.R_OK):
        raise ModuleRunnerError(f"模块文件无读取权限: {module_path}")
        raise ModuleRunnerError(f"Module file has no read permission: {module_path}")

    return module_path


def build_command(module_name: str, module_path: str, anaconda_path, env_name: str, custom_args: List[str], module_log_config='') -> str:
    """构建安全的执行命令"""


    # 合并默认参数和自定义参数（自定义参数优先级更高）
    #default_args = MODULE_CONFIG[module_name]["default_args"]
    #final_args = default_args + custom_args

    # 安全转义所有参数，防止命令注入
    escaped_args = [shlex.quote(arg) for arg in custom_args]
    args_str = " ".join(escaped_args)

    # 构建命令（使用set -e确保任一命令失败即退出）
    if anaconda_path is not None:
        anaconda_path = os.path.expanduser(anaconda_path)
        command = f"""
            #!/bin/bash
            set -euo pipefail
            PS1="${{PS1:-}}"
            # Load conda environment
            if [ -f "{shlex.quote(anaconda_path)}/etc/profile.d/conda.sh" ]; then
                source "{shlex.quote(anaconda_path)}/etc/profile.d/conda.sh"
            elif [ -f "{shlex.quote(anaconda_path)}/bin/activate" ]; then
                source "{shlex.quote(anaconda_path)}/bin/activate"
            else
                echo "Conda activation script not found" >&2
                exit 1
            fi

            # Activate environment and run module
            conda activate {shlex.quote(env_name)}
            python -u {shlex.quote(module_path)} {args_str} {module_log_config}
            """
    else:
        command = f"""
            # Activate environment and run module
            conda run -n {shlex.quote(env_name)} python -u {shlex.quote(module_path)} {args_str} {module_log_config}
            """

    return command

def run_command(command):
    # 创建子进程，捕获标准输出和错误
    logger.info('*'*10)
    logger.info(f"Now starting to execute the command:\n{command}")
    logger.info('*'*10)
    process = subprocess.Popen(
            command,
            shell=True,
            executable="/bin/bash",
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding='utf-8',  # 显式指定编码，解决子进程输出中文乱码
            errors='ignore'  # 忽略无法解码的字符（避免崩溃
        )

    # 实时打印输出的函数
    def print_output():
        for line in iter(process.stdout.readline, ''):
            # 移除行尾换行符后打印
            #print(line, end='')
            logger.info(line)
            sys.stdout.flush()  # 确保立即显示
        process.stdout.close()

    # 启动输出打印线程
    output_thread = threading.Thread(target=print_output)
    output_thread.daemon = True  # 主程序退出时自动结束线程
    output_thread.start()
    # 等待进程结束
    process.wait()

    # 检查退出状态
    if process.returncode == 0:
        # 场景1：退出码0，完全正常执行，无额外操作，正常返回
        logger.info("\n=== 命令执行成功 ===")
        logger.info("\n=== Command executed successfully ===")
    elif process.returncode == 100:
        # 场景2：退出码100，约定正常终止（无有效结果），不抛异常，提示信息
        #print("\n=== 命令执行完成，正常终止===")
        logger.info(f"退出码：{process.returncode}")
        logger.info(f"Exit code: {process.returncode}")
        sys.exit(0)
    elif process.returncode == 101:
        logger.info(f"退出码：{process.returncode}")
        logger.info(f"Exit code: {process.returncode}")
        sys.exit(0)
    elif process.returncode == 102:
        logger.info(f"退出码：{process.returncode}")
        logger.info(f"Exit code: {process.returncode}")
        sys.exit(0)
    else:
        # 场景3：其他非0/非100退出码，真正的执行失败，抛出原有异常
        raise RuntimeError(f"Command execution failed，exit code: {process.returncode}")
    return


def run_module(
        module_name: str,
        anaconda_path,
        params,
        module_log_config='',
        retry_count: int = 0
) :
    """
    在指定Conda环境中运行模块（支持重试）

    Args:
        module_name: 模块名称
        args: 模块的命令行参数
        retry_count: 当前重试次数

    Returns:
        退出代码（0表示成功）

    Raises:
        ModuleRunnerError: 模块验证或运行失败时抛出
    """
    # 验证模块
    try:
        module_path = validate_module(module_name)
    except ModuleRunnerError as e:
        logger.error(f"模块验证失败: {e}")
        logger.error(f"Module validation failed: {e}")
        raise

    # 获取环境名称
    env_name = params['env_name']
    logger.info(f"🚀 启动模块: {module_name} (环境: {env_name}, 路径: {module_path})")
    logger.info(f"🚀 Starting module: {module_name} (Environment: {env_name}, Path: {module_path})")

    args = [elem for k, v in params['args'].items() for elem in (f'--{k}', str(v))]
    # 构建命令
    command = build_command(
        module_name=module_name,
        module_path=module_path,
        anaconda_path=anaconda_path,
        env_name=env_name,
        custom_args=list(args),
        module_log_config=module_log_config,
    )

    run_command(command)
    return



def run_module_old(
        module_name: str,
        anaconda_path,
        params,
        retry_count: int = 0
) -> int:
    """
    在指定Conda环境中运行模块（支持重试）

    Args:
        module_name: 模块名称
        args: 模块的命令行参数
        retry_count: 当前重试次数

    Returns:
        退出代码（0表示成功）

    Raises:
        ModuleRunnerError: 模块验证或运行失败时抛出
    """
    # 验证模块
    try:
        module_path = validate_module(module_name)
    except ModuleRunnerError as e:
        logger.error(f"模块验证失败: {e}")
        logger.error(f"Module validation failed: {e}")
        raise

    # 获取环境名称
    env_name = params['env_name']
    logger.info(f"🚀 启动模块: {module_name} (环境: {env_name}, 路径: {module_path})")
    logger.info(f"🚀 Starting module: {module_name} (Environment: {env_name}, Path: {module_path})")

    args = [elem for k, v in params['args'].items() for elem in (f'--{k}', str(v))]
    # 构建命令
    command = build_command(
        module_name=module_name,
        module_path=module_path,
        anaconda_path=os.path.expanduser(anaconda_path),
        env_name=env_name,
        custom_args=list(args)
    )

    try:
        # 执行命令
        result = subprocess.run(
            command,
            shell=True,
            executable="/bin/bash",
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
            timeout=CONFIG["COMMAND_TIMEOUT"]
        )

        # 记录输出
        logger.info(f"=== 模块 {module_name} 输出 ===")
        logger.info(f"=== Module {module_name} output ===")
        # 正常停止逻辑：检测约定退出码，不触发异常
        if result.stdout:
            logger.info(result.stdout)
        if result.stderr:
            logger.error(f"模块 {module_name} 错误输出: {result.stderr}")
            logger.error(f"Module {module_name} error output: {result.stderr}")

        logger.info(f"模块 {module_name} 退出代码: {result.returncode}")
        logger.info(f"Module {module_name} exit code: {result.returncode}")

        # 重试逻辑
        #if result.returncode != 0 and retry_count < CONFIG["MAX_RETRIES"]:
            #retry_count += 1
            #logger.warning(f"模块 {module_name} 运行失败，将进行第 {retry_count}/{CONFIG['MAX_RETRIES']} 次重试...")
            #return run_module(module_name, *args, retry_count=retry_count)

        return result.returncode

    except subprocess.TimeoutExpired:
        error_msg = f"模块 {module_name} 运行超时（{CONFIG['COMMAND_TIMEOUT']}秒）"
        logger.error(error_msg)
        logger.error(f"Module {module_name} execution timeout ({CONFIG['COMMAND_TIMEOUT']} seconds)")
        raise ModuleRunnerError(error_msg) from None
    except subprocess.CalledProcessError as e:
        error_msg = f"模块 {module_name} 运行失败: {e.stderr}"
        logger.error(error_msg)
        logger.error(f"Module {module_name} execution failed: {e.stderr}")
        raise ModuleRunnerError(error_msg) from e
    except Exception as e:
        error_msg = f"模块 {module_name} 运行异常: {str(e)}"
        logger.error(error_msg, exc_info=True)
        logger.error(f"Module {module_name} runtime exception: {str(e)}", exc_info=True)
        raise ModuleRunnerError(error_msg) from e


def read_yaml_file(yaml_path: str) -> dict:
    """
    读取YAML文件并返回字典格式数据

    Args:
        yaml_path: YAML文件的路径（相对路径或绝对路径）

    Returns:
        解析后的字典数据

    Raises:
        FileNotFoundError: 文件不存在
        yaml.YAMLError: YAML格式错误
        PermissionError: 无文件读取权限
    """
    # 转换为Path对象，方便路径处理
    file_path = Path(yaml_path)

    # 检查文件是否存在
    if not file_path.exists():
        raise FileNotFoundError(f"错误：文件不存在 → {yaml_path}")

    # 检查是否是文件（不是目录）
    if not file_path.is_file():
        raise IsADirectoryError(f"错误：{yaml_path} 是目录，不是文件")

    # 读取并解析YAML文件
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            # yaml.safe_load() 避免执行恶意代码，更安全
            data = yaml.safe_load(f)
        return data or {}
    except PermissionError:
        raise PermissionError(f"错误：无权限读取文件 → {yaml_path}")
    except yaml.YAMLError as e:
        raise yaml.YAMLError(f"错误：YAML格式无效 → {e}")
    except Exception as e:
        raise Exception(f"未知错误：{e}")

def merge_configs(config_path: str, setting_path: str) -> dict:
    """
    合并用户配置和系统配置
    
    Args:
        config_path: 用户配置文件路径
        setting_path: 系统配置文件路径
        
    Returns:
        合并后的配置字典
    """
    # 读取配置文件
    user_config = read_yaml_file(config_path)
    setting_config = read_yaml_file(setting_path)
    
    # 合并配置
    merged = {}
    
    # 转换为模块配置
    #merged["modules"] = convert_to_module_config(user_config, setting_config)
    global_parameters = {}
    modules = {}
    project = user_config.get("project", {})
    profile = user_config.get("profile")

    rfdiffusion = user_config.get("rfdiffusion")
    mpnn = user_config.get("mpnn")
    mmseqs = user_config.get("mmseqs")
    esmfold = user_config.get("esmfold")
    alphafold2 = user_config.get("alphafold2")
    output_dir = project.get("output_dir", "./output")

    hmmer_setting = setting_config.get("hmmer", {})  # 无"hmmer"则返回{}
    hmmer_args = hmmer_setting.get("args", {})  # 无"args"则返回{}
    hmmer_user = profile or {}
    hmmer_args.update(hmmer_user)

    rfdiffusion_setting = setting_config.get("rfdiffusion", {})
    rfdiffusion_args = rfdiffusion_setting.get("args", {})
    rfdiffusion_user = rfdiffusion or {}
    rfdiffusion_args.update(rfdiffusion_user)

    mpnn_setting = setting_config.get("mpnn", {})
    mpnn_args = mpnn_setting.get("args", {})
    mpnn_user = mpnn or {}
    mpnn_args.update(mpnn_user)

    mmseqs_setting = setting_config.get("mmseqs", {})
    mmseqs_args = mmseqs_setting.get("args", {})
    mmseqs_user = mmseqs or {}
    mmseqs_args.update(mmseqs_user)


    esmfold_setting = setting_config.get("esmfold", {})
    esmfold_args = esmfold_setting.get("args", {})
    esmfold_user = esmfold or {}
    esmfold_args.update(esmfold_user)

    alphafold2_setting = setting_config.get("alphafold2", {})
    alphafold2_args = alphafold2_setting.get("args", {})
    alphafold2_user = alphafold2 or {}
    alphafold2_args.update(alphafold2_user)



    main_env = setting_config["environments"]["main_env"]

    protein_file = project.get("protein_file", '')
    proteinfile = os.path.basename(protein_file)
    protein_name = os.path.splitext(proteinfile)[0]
    os.makedirs(output_dir, exist_ok=True)
    if not os.path.isfile(protein_file):
        raise ValueError(f"蛋白质文件不存在，请检查路径： {protein_file}")
        raise ValueError(f"Protein file does not exist, please check the path: {protein_file}")
    if protein_file.endswith('.pdb'):
        input_pdb = f'{output_dir}/{protein_name}.pdb'
        print(f'蛋白质文件是pdb文件，将该文件复制到工作目录中：{input_pdb}')
        print(f'Protein file is in PDB format, copying to working directory: {input_pdb}')
        shutil.copy2(input_file, input_pdb)
    elif protein_file.endswith('.cif'):
        input_pdb = f'{output_dir}/{protein_name}.pdb'
        print(f'蛋白质文件是cif文件，将该文件转换为pdb文件，保存路径为：{input_pdb}\n')
        print(f'Protein file is in CIF format, converting to PDB format, save path: {input_pdb}\n')
        cif_to_pdb_biopython(
            cif_file_path=protein_file,
            pdb_file_path=input_pdb,
        )
    else:
        raise ValueError(f'蛋白质文件类型错误，目前仅支持.pdb和.cif文件，请检查输入文件：{input_file}')
        raise ValueError(f'Incorrect protein file type, currently only .pdb and .cif files are supported, please check input file: {input_file}')



    # 全局参数配置 (profile)
    if project.get("anaconda_path") is not None:
        global_parameters['anaconda_path'] = project['anaconda_path']
    global_parameters['work_dir'] = output_dir
    merged['global parameters'] = global_parameters

    chain = project.get("chain", "A")

    log_config = {}

    # 1. 获取当前脚本的路径（可能是相对路径）
    script_path = __file__
    # 2. 转换为绝对路径（避免相对路径的歧义）
    absolute_script_path = os.path.abspath(script_path)
    # 3. 提取文件所在的目录路径
    script_dir_path = os.path.dirname(absolute_script_path)
    print(f'Segdesign.py所在的目录：{script_dir_path}\n')
    print(f'Directory where Segdesign.py is located: {script_dir_path}\n')

    # hmmer 配置 (profile)
    if profile is not None:
        hmmer_env = setting_config["environments"].get("hmmer", main_env)
        if hmmer_env is None:
            hmmer_env = main_env
        hmmer_output_folder = os.path.join(output_dir, hmmer_args.get("output_folder", "hmmer_out"))
        hmmer_bitscore = hmmer_args.get("bitscore", 0.3)
        hmmer_n_iter = hmmer_args.get("n_iter", 5)
        hmmer_database = hmmer_args.get("database", "")
        hmmer_cpu = hmmer_args.get("cpu", 10)
        hmmer_minimum_sequence_coverage = hmmer_args.get("minimum_sequence_coverage", 50)
        hmmer_minimum_column_coverage = hmmer_args.get("minimum_column_coverage", 70)
        identity = hmmer_args.get("identity", 0.3)
        modules["hmmer"] = {
            "env_name": hmmer_env,
            "args": {
                "input_file": input_pdb,
                "select_chain": chain,
                "output_folder": hmmer_output_folder,
                "bitscore": hmmer_bitscore,
                "n_iter": hmmer_n_iter,
                "database": hmmer_database,
                "cpu": hmmer_cpu,
                "minimum_sequence_coverage": hmmer_minimum_sequence_coverage,
                "minimum_column_coverage": hmmer_minimum_column_coverage,
                "identity": identity,
                "final_report_folder": output_dir,  # 新增：最终报告输出到总工作目录
            }
        }

        log_config['hmmer'] = " 2>&1 | tee " + os.path.join(hmmer_output_folder,'hmmer_out.log')


    if project.get("segment") is not None:

        # rfdiffusion 配置
        if rfdiffusion is not None:

            run_inference_path = rfdiffusion_args["run_inference_path"]
            print(f'正在检测 run_inference.py 的路径是否正确... ')
            print(f'Checking if run_inference.py path is correct... ')
            if not os.path.isfile(run_inference_path):
                if not os.path.isabs(run_inference_path):
                    run_inference_path = os.path.join(script_dir_path, run_inference_path)
                    if os.path.isfile(run_inference_path):
                        print(f'检测完毕，路径正确\n')
                        print(f'Check completed, path is correct\n')
                    else:
                        raise ValueError(f'run_inference_path路径错误，请检查：{run_inference_path}')
                        raise ValueError(f'run_inference_path path is incorrect, please check: {run_inference_path}')
                else:
                    raise ValueError(f'run_inference_path路径错误，请检查：{run_inference_path}')
                    raise ValueError(f'run_inference_path path is incorrect, please check: {run_inference_path}')
            else:
               print(f'检测完毕，路径正确\n')
               print(f'Check completed, path is correct\n')

            rfdiffusion_output_folder = os.path.join(output_dir, rfdiffusion_args.get("output_folder", "rfdiffusion_out"))
            output_prefix = os.path.join(rfdiffusion_output_folder, f"sample/{protein_name}_{chain}")
            num_designs = rfdiffusion_args.get("num_designs", 10)
            contigs = f"[{chain}1-{project.get('sequence_length', '')}]"
            inpaint_str = f"[{chain}{project.get('segment', '')}]"
            partial_T = rfdiffusion_args["diffuser.partial_T"]
            rfdiffusion_env = setting_config["environments"]["rfdiffusion"]

            modules["rfdiffusion"] = {
                "env_name": rfdiffusion_env,
                "args": {
                    "run_inference_path": run_inference_path,
                    "inference.input_pdb": input_pdb,
                    "inference.output_prefix": output_prefix,
                    "inference.num_designs": num_designs,
                    "contigmap.contigs": contigs,
                    "contigmap.inpaint_str": inpaint_str,
                    "diffuser.partial_T": partial_T,
                }
            }
            if rfdiffusion_args.get("contigmap.inpaint_seq") is not None:
                modules["rfdiffusion"]["args"]["contigmap.inpaint_seq"] = rfdiffusion_args.get("contigmap.inpaint_seq")

            # RFdiffusion_report 配置
            rfdiffusion_report_env = setting_config["environments"].get("rfdiffusion_report", main_env)
            if rfdiffusion_report_env is None:
                rfdiffusion_report_env = main_env
            rfdiffusion_threshold = rfdiffusion_args.get("threshold", 0.6)
            modules["rfdiffusion_report"] = {
                "env_name": rfdiffusion_report_env,
                "args": {
                    "input_pdb": input_pdb,
                    "rfdiffusion_prefix": output_prefix,
                    "inpaint_str": inpaint_str,
                    "threshold": rfdiffusion_threshold,
                    "final_report_folder": output_dir,  # 新增：最终报告输出到总工作目录
                }

            }

            # 添加结构约束
            select_helix = rfdiffusion_args.get("helix")
            select_strand = rfdiffusion_args.get("strand")
            if select_helix and select_strand is not True:
                modules["rfdiffusion"]["args"]["contigmap.inpaint_str_helix"] = \
                    f"[{chain}{project.get('segment', '')}]"
                modules["rfdiffusion_report"]["args"]['ss'] = f"helix"
            elif select_strand and select_helix is not True:
                modules["rfdiffusion"]["args"]["contigmap.inpaint_str_strand"] = \
                    f"[{chain}{project.get('segment', '')}]"
                modules["rfdiffusion_report"]["args"]['ss'] = "strand"
            else:
                raise ModuleRunnerError(
                    f"Abnormal setting of secondary structure in the design area of module rfdiffusion")

            log_config['rfdiffusion'] = " 2>&1 | tee "  + os.path.join(rfdiffusion_output_folder, 'rfdiffusion_out.log')
            log_config['rfdiffusion_report'] = " 2>&1 | tee -a " + os.path.join(rfdiffusion_output_folder, 'rfdiffusion_out.log')



        # mpnn 配置
        if mpnn is not None:
            mpnn_env = setting_config["environments"]["mpnn"]

            parse_multiple_chains_path = mpnn_args["parse_multiple_chains_path"]
            assign_fixed_chains_path = mpnn_args["assign_fixed_chains_path"]
            make_fixed_positions_dict_path = mpnn_args["make_fixed_positions_dict_path"]
            protein_mpnn_run_path = mpnn_args["protein_mpnn_run_path"]


            print(f'正在检测 parse_multiple_chains_path 的路径是否正确... ')
            print(f'Checking if parse_multiple_chains_path path is correct... ')
            if not os.path.isfile(parse_multiple_chains_path):
                if not os.path.isabs(parse_multiple_chains_path):
                    parse_multiple_chains_path = os.path.join(script_dir_path, parse_multiple_chains_path)
                    if os.path.isfile(parse_multiple_chains_path):
                        print(f'检测完毕，路径正确\n')
                        print(f'Check completed, path is correct\n')
                    else:
                        raise ValueError(f'parse_multiple_chains_path路径错误，请检查：{parse_multiple_chains_path}')
                        raise ValueError(f'parse_multiple_chains_path path is incorrect, please check: {parse_multiple_chains_path}')
                else:
                    raise ValueError(f'parse_multiple_chains_path路径错误，请检查：{parse_multiple_chains_path}')
                    raise ValueError(f'parse_multiple_chains_path path is incorrect, please check: {parse_multiple_chains_path}')
            else:
                print(f'检测完毕，路径正确\n')
                print(f'Check completed, path is correct\n')

            print(f'正在检测 assign_fixed_chains_path 的路径是否正确... ')
            print(f'Checking if assign_fixed_chains_path path is correct... ')
            if not os.path.isfile(assign_fixed_chains_path):
                if not os.path.isabs(assign_fixed_chains_path):
                    assign_fixed_chains_path = os.path.join(script_dir_path, assign_fixed_chains_path)
                    if os.path.isfile(assign_fixed_chains_path):
                        print(f'检测完毕，路径正确\n')
                        print(f'Check completed, path is correct\n')
                    else:
                        raise ValueError(f'assign_fixed_chains_path路径错误，请检查：{assign_fixed_chains_path}')
                        raise ValueError(f'assign_fixed_chains_path path is incorrect, please check: {assign_fixed_chains_path}')
                else:
                    raise ValueError(f'assign_fixed_chains_path路径错误，请检查：{assign_fixed_chains_path}')
                    raise ValueError(f'assign_fixed_chains_path path is incorrect, please check: {assign_fixed_chains_path}')
            else:
                print(f'检测完毕，路径正确\n')
                print(f'Check completed, path is correct\n')

            print(f'正在检测 make_fixed_positions_dict_path 的路径是否正确... ')
            print(f'Checking if make_fixed_positions_dict_path path is correct... ')
            if not os.path.isfile(make_fixed_positions_dict_path):
                if not os.path.isabs(make_fixed_positions_dict_path):
                    make_fixed_positions_dict_path = os.path.join(script_dir_path, make_fixed_positions_dict_path)
                    if os.path.isfile(make_fixed_positions_dict_path):
                        print(f'检测完毕，路径正确\n')
                        print(f'Check completed, path is correct\n')
                    else:
                        raise ValueError(f'make_fixed_positions_dict_path路径错误，请检查：{make_fixed_positions_dict_path}')
                        raise ValueError(f'make_fixed_positions_dict_path path is incorrect, please check: {make_fixed_positions_dict_path}')
                else:
                    raise ValueError(f'make_fixed_positions_dict_path路径错误，请检查：{make_fixed_positions_dict_path}')
                    raise ValueError(f'make_fixed_positions_dict_path path is incorrect, please check: {make_fixed_positions_dict_path}')
            else:
                print(f'检测完毕，路径正确\n')
                print(f'Check completed, path is correct\n')

            print(f'正在检测 protein_mpnn_run_path 的路径是否正确... ')
            print(f'Checking if protein_mpnn_run_path path is correct... ')
            if not os.path.isfile(protein_mpnn_run_path):
                if not os.path.isabs(protein_mpnn_run_path):
                    protein_mpnn_run_path = os.path.join(script_dir_path, protein_mpnn_run_path)
                    if os.path.isfile(protein_mpnn_run_path):
                        print(f'检测完毕，路径正确\n')
                        print(f'Check completed, path is correct\n')
                    else:
                        raise ValueError(f'protein_mpnn_run_path路径错误，请检查：{protein_mpnn_run_path}')
                        raise ValueError(f'protein_mpnn_run_path path is incorrect, please check: {protein_mpnn_run_path}')
                else:
                    raise ValueError(f'protein_mpnn_run_path路径错误，请检查：{protein_mpnn_run_path}')
                    raise ValueError(f'protein_mpnn_run_path path is incorrect, please check: {protein_mpnn_run_path}')
            else:
                print(f'检测完毕，路径正确\n')
                print(f'Check completed, path is correct\n')


            if mpnn_args.get("pdb_folder") is not None:
                pdb_foler = mpnn_args.get("pdb_folder")
            else:
                pdb_foler = os.path.join(output_dir, f"rfdiffusion_out/filter_results")
            mpnn_output_folder = os.path.join(output_dir, mpnn_args.get("output_folder","mpnn_out"))
            chain_list = chain
            position_list =  f"{chain}{project.get('segment', '')}"
            num_seq_per_target = mpnn_args.get("num_seq_per_target", 20)
            sampling_temp = mpnn_args.get("sampling_temp", 0.3)
            seed = mpnn_args.get("seed", 42)
            batch_size = mpnn_args.get("batch_size", 1)

            modules["mpnn"] = {
                "env_name": mpnn_env,
                "args": {
                    "parse_multiple_chains_path": parse_multiple_chains_path,
                    "assign_fixed_chains_path": assign_fixed_chains_path,
                    "make_fixed_positions_dict_path": make_fixed_positions_dict_path,
                    "protein_mpnn_run_path": protein_mpnn_run_path,
                    "pdb_folder": pdb_foler,
                    "output_folder": mpnn_output_folder,
                    "chain_list": chain_list,
                    "position_list": position_list,
                    "num_seq_per_target": num_seq_per_target,
                    "sampling_temp": sampling_temp,
                    "seed": seed,
                    'batch_size': batch_size,
                    #"top_percent": int(proteinmpnn.get("threshold", 0.9))
                }
            }

            # mpnn_report 配置
            mpnn_report_env = setting_config["environments"].get("mpnn_report", main_env)
            if mpnn_report_env is None:
                mpnn_report_env = main_env
            seq_folder = os.path.join(mpnn_output_folder, "seqs")
            mpnn_report_output_folder = mpnn_output_folder
            top_percent = mpnn_args.get("top_percent", 0.5)
            rfdiffusion_report_path = mpnn_args.get("rfdiffusion_report_path")

            modules["mpnn_report"] = {
                "env_name": mpnn_report_env,
                "args": {
                    "seq_folder": seq_folder,
                    "output_folder": mpnn_report_output_folder,
                    "top_percent": top_percent,
                    "generate_report": True,  # 添加生成报告标志
                    "final_report_folder": output_dir,  # 新增：最终报告输出到总工作目录
                    "rfdiffusion_report_path": rfdiffusion_report_path,
                    "position_list": position_list,
                    'protein_pdb': input_pdb

                }
            }
            log_config['mpnn'] = " 2>&1 | tee " + os.path.join(mpnn_report_output_folder, 'mpnn_out.log')
            log_config['mpnn_report'] = " 2>&1 | tee -a " + os.path.join(mpnn_report_output_folder, 'mpnn_out.log')



        # 聚类分析配置
        if mmseqs is not None:
            mmseqs_env = setting_config["environments"].get("mmseqs", main_env)
            if mmseqs_env is None:
                mmseqs_env = main_env
            if mmseqs_args.get("input_folder") is not None:
                mmseqs_input_folder = mmseqs_args.get("input_folder")
            else:
                if mpnn is not None:
                    mmseqs_input_folder = os.path.join(mpnn_report_output_folder, 'top_filter')
                else:
                    mmseqs_input_folder =os.path.join(output_dir, 'mpnn_out', 'top_filter')

            mmseqs_output_folder = os.path.join(output_dir, mmseqs_args.get("output_folder", "mmseqs_out"))
            threads = mmseqs_args.get("threads", 8)
            min_seq_id = mmseqs_args.get("min_seq_id")
            cov_mode = mmseqs_args.get("cov_mode", 0)
            coverage = mmseqs_args.get("c", mmseqs_args.get("coverage", 0.8))
            mmseqs_path = mmseqs_args.get("mmseqs_path")
            sensitivity = mmseqs_args.get("s", mmseqs_args.get("sensitivity", 4.0))
            position_list = f"{chain}{project.get('segment', '')}"

            modules["mmseqs"] = {
                "env_name": mmseqs_env,
                "args": {
                    'input_folder': mmseqs_input_folder,
                    "output_folder": mmseqs_output_folder,
                    "position_list": position_list,
                    "threads": threads,
                    "min_seq_id": min_seq_id,
                    "cov_mode": cov_mode,
                    "coverage": coverage,
                    "mmseqs_path": mmseqs_path,
                    "sensitivity": sensitivity,
                }
            }

            log_config['mmseqs'] = " 2>&1 | tee " + os.path.join(mmseqs_output_folder, 'mmseqs_out.log')

            '''
            modules["mpnn_report"] = {
                "env_name": mpnn_report_env,
                "args": {
                    "seq_folder": seq_folder,
                    "output_folder": mpnn_report_output_folder,
                    "top_percent": top_percent,
                    "position_list": position_list,
                    "threads": threads,
                    "min_seq_id": min_seq_id,
                    "cov_mode": cov_mode,
                    "coverage": coverage,
                    "mmseqs_path": mmseqs_path,
                    "sensitivity": sensitivity,

                }
            }
            '''

        # esmfold 配置
        if esmfold is not None:
            esmfold_env = setting_config["environments"]["esmfold"]
            if esmfold_args.get("input_folder") is not None:
                esmfold_input_folder = esmfold_args.get("input_folder")
            else:
                esmfold_input_folder = os.path.join(output_dir, f"mmseqs_out/results")
            esmfold_output_folder = os.path.join(output_dir, esmfold_args.get("output_folder","esmfold_out"))
            mmseqs_report_path = os.path.join(output_dir, "mmseqs_report.csv")


            modules["esmfold"] = {
                "env_name": esmfold_env,
                "args": {
                    "input_folder": esmfold_input_folder,
                    "output_folder": esmfold_output_folder,
                    "mmseqs_report_path": mmseqs_report_path,
                }
            }

            # esmfold_report 配置
            esmfold_report_env = setting_config["environments"].get("esmfold_report", main_env)
            if esmfold_report_env is None:
                esmfold_report_env = main_env
            fasta_folder = esmfold_input_folder
            esmfold_folder = esmfold_output_folder
            plddt_threshold = esmfold_args.get("plddt_threshold")
            ptm_threshold = esmfold_args.get("ptm_threshold")
            if esmfold_args.get("original_protein_chain_path") is not None:
                original_protein_chain_path = esmfold_args.get("original_protein_chain_path")
            else:
                chain_folder = os.path.join(output_dir, f"hmmer_out/target_chain_pdb")
                filenames = f"{protein_name}_{chain}.pdb"
                original_protein_chain_path = os.path.join(chain_folder, filenames)

            if esmfold_args.get("seq_range_str") is not None:
                seq_range_str = esmfold_args.get("seq_range_str")
            else:
                seq_range_str = project.get("segment")

            esmfold_ss = esmfold_args.get("ss")
            if esmfold_ss is not None:
                ss_threshold = esmfold_args.get("ss_threshold")
                if ss_threshold is not None:
                    pass
                else:
                    ss_threshold = rfdiffusion_args.get("threshold", 0.6)
                modules["esmfold_report"] = {
                    "env_name": esmfold_report_env,
                    "args": {
                        "esmfold_folder": esmfold_folder,
                        "original_protein_chain_path": original_protein_chain_path,
                        "seq_range_str": seq_range_str,
                        'ss': esmfold_ss,
                        'ss_threshold': ss_threshold,
                    }
                }
            else:
                if rfdiffusion is not None:
                    esmfold_ss = modules["rfdiffusion_report"]["args"]['ss']
                    ss_threshold = rfdiffusion_threshold
                    modules["esmfold_report"] = {
                        "env_name": esmfold_report_env,
                        "args": {
                            "esmfold_folder": esmfold_folder,
                            "original_protein_chain_path": original_protein_chain_path,
                            "seq_range_str": seq_range_str,
                            'ss': esmfold_ss,
                            'ss_threshold': ss_threshold,
                        }
                    }
                else:
                    modules["esmfold_report"] = {
                        "env_name": esmfold_report_env,
                        "args": {
                            "esmfold_folder": esmfold_folder,
                            "original_protein_chain_path": original_protein_chain_path,
                            "seq_range_str": seq_range_str,
                        }
                    }
                if ptm_threshold is not None:
                    modules['esmfold_report']['args']['ptm_threshold'] = ptm_threshold
                if plddt_threshold is not None:
                    modules['esmfold_report']['args']['plddt_threshold'] = plddt_threshold

            esmfold_ss_filter = esmfold_args.get("ss_filter", True)
            if esmfold_ss_filter is None:
                esmfold_ss_filter = True
            modules["esmfold_report"]["args"]["ss_filter"] = esmfold_ss_filter

            log_config['esmfold'] = " 2>&1 | tee " + os.path.join(esmfold_output_folder, 'esmfold_out.log')
            log_config['esmfold_report'] = " 2>&1 | tee -a " + os.path.join(esmfold_output_folder, 'esmfold_out.log')
        
        # alphafold2 配置
        if alphafold2 is not None:
            alphafold2_env = setting_config["environments"]["alphafold2"]
            if alphafold2_args.get("input_file") is not None:
                alphafold2_input_file = alphafold2_args.get("input_file")
            else:
                if esmfold is not None:
                    alphafold2_input_file = os.path.join(esmfold_output_folder, "filter_result.fa")
                else:
                    alphafold2_input_file = os.path.join(output_dir, "esmfold_out", "filter_result.fa")
            alphafold2_output_folder = os.path.join(output_dir, alphafold2_args.get("output_folder", "alphafold2_out"))
            esmfold_report_path = os.path.join(output_dir, "esmfold_report.csv")
            num_recycle = alphafold2_args.get("num_recycle", None)
            amber = alphafold2_args.get("amber", True)
            templates = alphafold2_args.get("templates", True)
            gpu = alphafold2_args.get("gpu", False)
            random_seed = alphafold2_args.get("random_seed", 0)

            modules["alphafold2"] = {
                "env_name": alphafold2_env,
                "args": {
                    "input_file": alphafold2_input_file,
                    "output_folder": alphafold2_output_folder,
                    "esmfold_report_path": esmfold_report_path,
                    "amber": amber,
                    "templates": templates,
                    'gpu': gpu,
                    'random_seed': random_seed,
                }
            }
            if num_recycle is not None:
                modules["alphafold2"]["args"]["num_recycle"] = num_recycle



            # alphafold2_report 配置
            alphafold2_report_env = setting_config["environments"].get("alphafold2_report", main_env)
            if alphafold2_report_env is None:
                alphafold2_report_env = main_env
            #fasta_folder = esmfold_input_folder
            alphafold2_folder = alphafold2_output_folder
            af2_plddt_threshold = alphafold2_args.get("plddt_threshold")
            af2_ptm_threshold = alphafold2_args.get("ptm_threshold")
            esmfold_report_path = os.path.join(output_dir, "esmfold_report.csv")
            

            if alphafold2_args.get("seq_range_str") is not None:
                seq_range_str = alphafold2_args.get("seq_range_str")
            else:
                seq_range_str = project.get("segment")

            af2_ss = alphafold2_args.get("ss")
            if af2_ss is not None:
                af2_ss_threshold = alphafold2_args.get("ss_threshold")
                if af2_ss_threshold is not None:
                    pass
                else:
                    af2_ss_threshold = rfdiffusion_args.get("threshold", 0.6)
                modules["alphafold2_report"] = {
                    "env_name": alphafold2_report_env,
                    "args": {
                        "esmfold_report_path":esmfold_report_path,
                        'alphafold2_folder': alphafold2_folder,
                        "seq_range_str": seq_range_str,
                        'ss': af2_ss,
                        'ss_threshold': af2_ss_threshold,
                    }
                }
            else:
                if rfdiffusion is not None:
                    af2_ss = modules["rfdiffusion_report"]["args"]['ss']
                    af2_ss_threshold = rfdiffusion_threshold
                    modules["alphafold2_report"] = {
                        "env_name": alphafold2_report_env,
                        "args": {
                            "esmfold_report_path":esmfold_report_path,
                            'alphafold2_folder': alphafold2_folder,
                            "seq_range_str": seq_range_str,
                            'ss': af2_ss,
                            'ss_threshold': af2_ss_threshold,
                        }
                    }
                else:
                    modules["alphafold2_report"] = {
                        "env_name": alphafold2_report_env,
                        "args": {
                            "esmfold_report_path":esmfold_report_path,
                            'alphafold2_folder': alphafold2_folder,
                            "seq_range_str": seq_range_str,
                        }
                    }

            if af2_ptm_threshold is not None:
                modules['alphafold2_report']['args']['ptm_threshold'] = af2_ptm_threshold
            if af2_plddt_threshold is not None:
                modules['alphafold2_report']['args']['plddt_threshold'] = af2_plddt_threshold
            af2_ss_filter = alphafold2_args.get("ss_filter", True)
            if af2_ss_filter is None:
                af2_ss_filter = True
            modules['alphafold2_report']['args']['ss_filter'] = af2_ss_filter

            log_config['alphafold2'] = " 2>&1 | tee " + os.path.join(alphafold2_output_folder, 'alphafold2_out.log')
            log_config['alphafold2_report'] = " 2>&1 | tee -a " + os.path.join(alphafold2_output_folder, 'alphafold2_out.log')

    # 聚类分析配置
    """
        if project.get("segment") is not None and mmseqs is not None:
        # 动态计算Top百分比文件夹路径
        top_percent_value = mpnn_args.get("top_percent", 0.5)
        top_percent_str = f"{top_percent_value*100:.1f}%"
        
        # 获取mpnn_output_folder，如果不存在则使用默认值
        mpnn_output_folder = os.path.join(output_dir, mpnn_args.get("output_folder", "mpnn_out"))
        top_sequences_folder = os.path.join(mpnn_output_folder, f"top_{top_percent_str}")
        
        # 解析区域位置
        position_range = project.get("segment", "")
        start_pos = int(position_range.split('-')[0]) if '-' in position_range else 1
        end_pos = int(position_range.split('-')[1]) if '-' in position_range else 100
        
        modules["cluster_analysis"] = {
            "env_name": setting_config["environments"].get("cluster_analysis", setting_config["environments"]["main_env"]),
            "args": {
                "input_folder": top_sequences_folder,
                "output_folder": os.path.join(output_dir, "cluster_analysis_out"),
                "start": start_pos,
                "end": end_pos,
                "min_seq_id": mmseqs_args.get("min_seq_id", 0.8),
                "cov_mode": mmseqs_args.get("cov_mode", 0),
                "coverage": mmseqs_args.get("coverage", 0.8),
                "mmseqs_path": mmseqs_args.get("mmseqs_path", "mmseqs"),
                "threads": mmseqs_args.get("threads", 8)
            }
        }

    """

    merged["modules"] = modules
    merged["log_config"] = log_config
    return merged


def global_work_dir_handling(yaml_data):
    """处理工作目录"""
    work_dir = os.path.expanduser(yaml_data.get('global parameters', {}).get("work_dir", "./output"))
    if not os.path.exists(work_dir):
        os.makedirs(work_dir, exist_ok=True)
    return work_dir


def cif_to_pdb_biopython(cif_file_path, pdb_file_path):
    """
    使用Biopython将CIF文件转换为PDB文件
    :param cif_file_path: 输入CIF文件路径（绝对/相对路径）
    :param pdb_file_path: 输出PDB文件路径（绝对/相对路径）
    """
    try:
        # 1. 初始化CIF解析器（QUIET=True关闭无关日志输出）
        cif_parser = MMCIFParser(QUIET=True)

        # 2. 解析CIF文件，获取结构对象（第一个参数为结构名称，可自定义）
        structure = cif_parser.get_structure("target_structure", cif_file_path)

        # 3. 初始化PDB写入器
        pdb_writer = PDBIO()

        # 4. 设置要写入的结构对象
        pdb_writer.set_structure(structure)

        # 5. 写入PDB文件（可选：select参数筛选原子，默认写入全部原子）
        pdb_writer.save(pdb_file_path)

        print(f"转换成功！PDB文件已保存至：{pdb_file_path}")
        print(f"Conversion successful! PDB file saved to: {pdb_file_path}")

    except FileNotFoundError:
        raise ValueError(f"错误：找不到CIF文件 '{cif_file_path}'")
        raise ValueError(f"Error: CIF file not found '{cif_file_path}'")
    except Exception as e:
        print(f"转换失败：{str(e)}")
        print(f"Conversion failed: {str(e)}")
    return




if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="SegDesign: 蛋白质设计工具",
        epilog="示例：python Segdesign.py --config ./config/config.yaml --setting ./config/setting.yaml"
    )

    # 添加参数
    parser.add_argument(
        "--config",
        type=str,
        default=CONFIG["CONFIG_PATH"]["MAIN"],
        help="用户配置文件路径（相对路径或绝对路径）"
    )
    parser.add_argument(
        "--setting",
        type=str,
        default=CONFIG["CONFIG_PATH"]["SETTING"],
        help="系统配置文件路径（相对路径或绝对路径）"
    )
    
    args = parser.parse_args()

    # 合并配置
    merged_config = merge_configs(args.config, args.setting)
    print("✅ 配置文件读取成功！")
    print("✅ Configuration files read successfully!")
    print("📊 解析后的数据：")
    print("📊 Parsed data:")
    print(yaml.dump(merged_config, allow_unicode=True, sort_keys=False))

    # 处理工作目录
    output_dir = global_work_dir_handling(merged_config)

    log_path = os.path.join(output_dir, "Segdesign.log")
    logger = setup_logger(log_path)

    logger.info(f"工作目录: {output_dir}")
    logger.info(f"Working directory: {output_dir}")

    # 将config.yaml复制到工作目录下

    shutil.copy(args.config, f"{output_dir}/config.yaml")

    # 获取anaconda路径
    anaconda_path = merged_config["global parameters"].get("anaconda_path")

    # 运行模块
    for module_name, params in merged_config["modules"].items():
        module_log_config = merged_config["log_config"][module_name]
        if module_name in CONFIG['MODULES']:
            try:
                logger.info(f"正在运行模块: {module_name}")
                logger.info(f"Running module: {module_name}")
                run_module(
                    module_name=module_name,
                    anaconda_path=anaconda_path,
                    params=params,
                    module_log_config=module_log_config
                )
                logger.info(f"✅ 模块 {module_name} 运行成功")
                logger.info(f"✅ Module {module_name} executed successfully")
            except ModuleRunnerError as e:
                logger.critical(f"❌ 模块 {module_name} 运行失败: {e}")
                logger.critical(f"❌ Module {module_name} execution failed: {e}")
                exit(1)
            except KeyboardInterrupt:
                logger.info("程序被用户中断")
                logger.info("Program interrupted by user")
                exit(0)
            except Exception as e:
                logger.critical(f"❌ 模块 {module_name} 未预期的错误: {str(e)}", exc_info=True)
                logger.critical(f"❌ Module {module_name} unexpected error: {str(e)}", exc_info=True)
                exit(1)
    logger.info("🎉 所有模块运行完成！")
    logger.info("🎉 All modules completed!")


