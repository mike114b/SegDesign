import os
import argparse
import shutil
import pandas as pd
import numpy as np
from typing import Dict

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re
from colabfold.batch import run
import sys
from contextlib import contextmanager


def parse_args():
    parser = argparse.ArgumentParser(description='Protein 3D Structure Prediction(alphafold2)',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--input_file', type=str,
                        help='The folder for storing sequence files or the path where the sequence file is located')
    parser.add_argument('--output_folder', type=str,
                        help='Folder for storing output files')
    parser.add_argument("--esmfold_report_path", type=str, default=None,
                        help="The path to esmfold_report.csv. If not entered, the default path will be used: {work_dir}/esmfold_report.csv")
    parser.add_argument("--num_recycle", type=int, default=None,
                        help="Number of model iteration optimizations (default is None, using the model's default value).")
    parser.add_argument("--amber", type=bool, default=True,
                        help="Whether to perform AMBER force field relaxation on the predicted structure, default is True. After enabling, only optimize the best model")
    parser.add_argument("--templates", type=bool, default=True,
                        help="Whether to use PDB templates (default is True). Setting it to True will search for homologous structure templates to improve prediction accuracy.")
    parser.add_argument("--gpu", type=bool, default=False,
                        help="Whether to use GPU to accelerate Relax (default is False)")
    parser.add_argument("--random_seed", type=int, default=0,
                        help="Random seed for reproducibility (default is 0)")

    
    return parser.parse_args()


#读取esmfold_report.csv文件，生成序列字典
def read_esmfold_report(esmfold_report_path) -> Dict[str, Dict[str, str]]:
    
    # 1. 输入参数类型校验
    if not isinstance(esmfold_report_path, str):
        raise TypeError(f"文件路径必须是字符串，当前类型：{type(esmfold_report_path)}")
    if not esmfold_report_path.strip():
        raise ValueError("文件路径不能为空字符串")

    # 2. 读取CSV文件，捕获解析相关异常
    try:
        df = pd.read_csv(esmfold_report_path)
    except FileNotFoundError:
        raise FileNotFoundError(f"未找到mmseqs报告文件：{esmfold_report_path}")
    except pd.errors.EmptyDataError:
        raise pd.errors.EmptyDataError(f"mmseqs报告文件为空：{esmfold_report_path}")
    except pd.errors.ParserError:
        raise pd.errors.ParserError(f"CSV格式错误，无法解析：{esmfold_report_path}")

    # 4. 筛选whether_pass为True的行（兼容布尔/字符串/数值类型，提升鲁棒性）
    filtered_df = df.copy()
    # 把所有'-'替换为False，再统一转布尔筛选
    filtered_df['whether_pass'] = filtered_df['whether_pass'].replace('-', False)
    filtered_df['whether_pass'] = filtered_df['whether_pass'].replace('False', False)
    filtered_df['whether_pass'] = filtered_df['whether_pass'].replace('True', True)
    filtered_df = filtered_df[filtered_df['whether_pass'].astype(bool)]
    #print('filtered_df:',filtered_df['whether_pass'])
    # print('filtered_df:', filtered_df)
    # 无符合条件数据时直接返回空字典
    if filtered_df.empty:
        return {}

    # 5. 构建结果字典（保持backbone的原始出现顺序，drop_duplicates比unique更直观）
    result_dict = {}
    # 按原始顺序获取唯一的backbone值
    unique_backbones = filtered_df["backbone"].drop_duplicates().tolist()
    for backbone in unique_backbones:
        # 筛选当前backbone的所有数据
        backbone_subdf = filtered_df[filtered_df["backbone"] == backbone]
        # 构建内层字典：index -> sequence
        inner_dict = {row["index"]: row["sequence"] for _, row in backbone_subdf.iterrows()}
        result_dict[backbone] = inner_dict
    #print('序列字典：',result_dict)
    return result_dict

#自然顺序排列
def natural_sort_key(filename):
    """生成自然排序的key：将文件名拆分为字符串和数字部分，数字转整数"""
    parts = re.split(r'(\d+)', os.path.splitext(filename)[0])
    key = []
    for part in parts:
        if part.isdigit():
            key.append(int(part))
        else:
            key.append(part)
    return key

#从fasta文件中读取数据，生成序列字典
def read_fasta_file(input_file: str) -> Dict[str, Dict[str, str]]:
    seqs_dict = {}
    filename = os.path.basename(input_file)
    if filename.endswith(".fasta") or filename.endswith(".fa"):
        for record in SeqIO.parse(input_file, "fasta"):
            seqs_dict[record.id] = record.seq

        #print('seqs_dict: ', seqs_dict)


        pattern = re.compile(r'^(.+)_(\d+)_(.+)$')

        # 2. 初始化分组字典：临时存储各组数据
        groups = {}
        protein_chains = {}  # 存储不符合格式的条目

        # 3. 遍历输入字典，按格式分组
        for name, seq in seqs_dict.items():
            match = pattern.match(name)
            if match:
                # 匹配成功：提取"任意字符_任意数字"作为组名（第一部分+第二部分）
                part1, part2, _ = match.groups()
                group_name = f"{part1}_{part2}"
                # 将当前名称和序列加入对应组
                if group_name not in groups:
                    groups[group_name] = {}
                groups[group_name][name] = seq
            else:
                # 匹配失败：加入protein_chains
                protein_chains[name] = seq

        # 4. 构建最终结果字典：所有组按自然顺序排序，组内元素也按自然顺序排序
        final_result = {}

        # 4.1 处理匹配成功的组：组名自然排序，组内名称自然排序
        sorted_group_names = sorted(groups.keys(), key=natural_sort_key)  # 组名自然排序
        for group in sorted_group_names:
            # 组内名称自然排序，构建有序子字典
            sorted_names = sorted(groups[group].keys(), key=natural_sort_key)
            final_result[group] = {name: groups[group][name] for name in sorted_names}

        # 4.2 处理protein_chains：名称自然排序
        if protein_chains:
            sorted_protein_names = nsorted(protein_chains.keys())
            final_result["protein_chains"] = {name: protein_chains[name] for name in sorted_protein_names}

    else:
        final_result = {}
        files = os.listdir(input_file)
        files = sorted(files, key=natural_sort_key)
        files = [file for file in files if file.endswith(".fasta") or file.endswith(".fa")]
        for file in files:
            backbone = os.path.splitext(file)[0]
            file_path = os.path.join(input_file, file)
            fasta_dict = {}
            for record in SeqIO.parse(file_path, "fasta"):
                fasta_dict[record.id] = record.seq

            fasta_ids = sorted(fasta_dict.keys(), key=natural_sort_key)
            seqs_dict = {name: fasta_dict[name] for name in fasta_ids}
            final_result[backbone] = seqs_dict

    print('final_result: ', final_result)
    return final_result


def run_alphafold2(seqs_dict,output_folder, templates, gpu, amber, random_seed, num_recycle=None):
    def run_function_silently(func, *args, **kwargs):
        """
        静默执行函数，抑制所有stdout/stderr输出
        参数:
            func: 要执行的函数
            *args: 函数位置参数
            **kwargs: 函数关键字参数
        返回:
            函数的返回值（如果有）
        """
        # 方案1：将输出捕获到内存（可后续查看，适合调试）
        # stdout_buf = io.StringIO()
        # stderr_buf = io.StringIO()

        # 方案2：直接丢弃输出（推荐，完全静默）
        stdout_buf = open(os.devnull, 'w')
        stderr_buf = open(os.devnull, 'w')

        try:
            # 重定向stdout和stderr
            with redirect_stdout(stdout_buf), redirect_stderr(stderr_buf):
                # 执行目标函数
                result = func(*args, **kwargs)
            return result
        finally:
            # 关闭文件/缓冲区（避免资源泄漏）
            stdout_buf.close()
            stderr_buf.close()

    for backbone, backbone_value in seqs_dict.items():
        backbone_folder = os.path.join(output_folder, 'colabfold_batch', backbone)
        if not os.path.exists(backbone_folder):
            os.makedirs(backbone_folder, exist_ok=True)
        for name, seq in backbone_value.items():
            seq_folder = os.path.join(backbone_folder, name)
            if not os.path.exists(seq_folder):
                os.makedirs(seq_folder, exist_ok=True)

            print(f'现在开始预测 {name} 的三维结构...')
            print(f'Starting 3D structure prediction for {name}...')
            if amber:
                if num_recycle is not None:
                    run(
                        queries=[(name, seq, None, None)],  # 输入序列列表
                        result_dir=seq_folder,  # 输出目录
                        num_recycle=num_recycle,
                        num_models=5,
                        is_complex=False,
                        random_seed=random_seed,
                        num_relax=1,
                        use_templates=templates,  # 是否使用模板
                        use_gpu_relax=gpu,  # GPU 松弛
                    )
                else:
                    run(
                        queries=[(name, seq, None, None)],  # 输入序列列表
                        result_dir=seq_folder,  # 输出目录
                        num_models=5,
                        is_complex=False,
                        random_seed=random_seed,
                        num_relax=1,
                        use_templates=templates,  # 是否使用模板
                        use_gpu_relax=gpu,  # GPU 松弛
                    )

            else:
                if num_recycle is not None:
                    run(
                        queries=[(name, seq, None, None)],  # 输入序列列表
                        result_dir=seq_folder,  # 输出目录
                        num_recycle=num_recycle,
                        num_models=5,
                        is_complex=False,
                        random_seed=random_seed,
                        use_templates=templates,  # 是否使用模板
                        use_gpu_relax=gpu,  # GPU 松弛
                    )

                else:
                    run(
                        queries=[(name, seq, None, None)],  # 输入序列列表
                        result_dir=seq_folder,  # 输出目录
                        num_models=5,
                        is_complex=False,
                        random_seed=random_seed,
                        use_templates=templates,  # 是否使用模板
                        use_gpu_relax=gpu,  # GPU 松弛
                    )

            print(f'序列 {name} 预测完毕')
            print(f'Sequence {name} prediction completed')

    return

def extract_best_model(colabfold_folder, out_folder, result_dict=None):
    if result_dict is None:
        backbone_folders = os.listdir(colabfold_folder)
        backbone_folders = sorted(backbone_folders, key=natural_sort_key)
    else:
        backbone_folders = result_dict.keys()

    for backbone in backbone_folders:
        if result_dict is None:
            seq_folders = os.listdir(os.path.join(colabfold_folder, backbone))
            seq_folders = sorted(seq_folders, key=natural_sort_key)
        else:
            seq_folders = result_dict[backbone].keys()

        for seq_folder in seq_folders:
            folder_path = os.path.join(colabfold_folder, backbone, seq_folder)

            pattern0 = re.compile(r'.*_relaxed_rank_001_.*\.pdb$')
            pattern1 = re.compile(r'.*_scores_rank_001_.*\.json$')

            # 存储匹配的文件路径
            pdb_file = ''
            json_file = ''

            # 遍历目标目录下的所有条目
            for entry in os.listdir(folder_path):
                # 拼接完整路径
                entry_path = os.path.join(folder_path, entry)
                # 1. 仅处理文件（排除目录）；2. 文件名匹配正则模式
                if os.path.isfile(entry_path) and pattern0.match(entry):
                    pdb_file=entry_path
                elif os.path.isfile(entry_path) and pattern1.match(entry):
                    json_file=entry_path

            print(f"找到评分最高的pdb文件：")
            print(f"Found the highest scoring pdb file:")
            print(f"{pdb_file}")

            new_file_folder = os.path.join(out_folder, backbone)
            new_file_path = os.path.join(out_folder, backbone, f'{seq_folder}.pdb')
            new_json_file = os.path.join(out_folder, backbone, f'{seq_folder}.json')
            os.makedirs(new_file_folder, exist_ok=True)
            shutil.copy(pdb_file, new_file_path)
            print(f'将评分最高的pdb文件拷贝到路径：\n{new_file_path}')
            print(f'Copying the highest scoring pdb file to path:\n{new_file_path}')
            shutil.copy(json_file, new_json_file)
    return


def main():
    args = parse_args()
    input_file = os.path.expanduser(args.input_file)
    output_folder = os.path.expanduser(args.output_folder)
    esmfold_report_path = os.path.expanduser(args.esmfold_report_path)
    random_seed = args.random_seed
    num_recycle = args.num_recycle
    templates = args.templates
    gpu = args.gpu
    amber = args.amber

    work_dir = output_folder.rsplit('/', 1)[0]
    if not os.path.isfile(esmfold_report_path):
        esmfold_report_path_default = os.path.join(work_dir, 'esmfold_report.csv')
        if not os.path.isfile(esmfold_report_path_default):
            seqs_dict = read_fasta_file(input_file)
        else:
            esmfold_report_path = esmfold_report_path_default
            seqs_dict = read_esmfold_report(esmfold_report_path)
    else:
        seqs_dict = read_esmfold_report(esmfold_report_path)

    run_alphafold2(
        seqs_dict=seqs_dict,
        output_folder=output_folder,
        templates=templates,
        gpu=gpu,
        amber=amber,
        random_seed=random_seed,
        num_recycle=num_recycle,
    )

    extract_best_model(
        colabfold_folder=os.path.join(output_folder, 'colabfold_batch'),
        out_folder=os.path.join(output_folder, 'structure_prediction_files'),
        result_dict=seqs_dict,
    )

if __name__ == '__main__':
    main()












                








    






