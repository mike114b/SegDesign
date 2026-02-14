#!/usr/bin/env python3
"""
MMseqs2 聚类分析工具
功能：从指定文件夹读取FASTA文件，对序列进行聚类分析，生成mmseqs_report.csv报告
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import re


def arg_parser():
    parser = argparse.ArgumentParser(
        description="Perform MMseqs2 clustering on FASTA files and generate mmseqs_report.csv"
    )
    parser.add_argument("--input_folder", required=True, type=Path,
                        help="Input folder containing FASTA files")
    parser.add_argument("--output_folder", required=True, type=Path,
                        help="Output folder")
    parser.add_argument("-t", "--threads", type=int, default=8,
                        help="Number of threads for MMseqs2 (default: 8)")
    parser.add_argument("--min_seq_id", type=float, default=0.8,
                        help="Minimum sequence identity (default: 0.8)")
    parser.add_argument("--cov_mode", type=int, default=0,
                        help="Coverage mode (0=bidirectional, 1=query, default: 0)")
    parser.add_argument("-c", "--coverage", type=float, default=0.8,
                        help="Coverage threshold (default: 0.8)")
    parser.add_argument("-s", "--sensitivity", type=float, default=4.0,
                        help="Sensitivity: 1.0 faster; 4.0 fast; 7.5 sensitive [4.000]")
    parser.add_argument("--mmseqs_path", type=str, default="mmseqs",
                        help="Path to mmseqs command (default: mmseqs)")
    parser.add_argument('--position_list', type=str, default=None,
                        help='Specified design area, such as A1-5 or 1-5')
    return parser.parse_args()


def extract_subregions(
        input_file: Path,
        output_fasta: Path,
        start_pos: int,
        end_pos: int,
) -> Dict[str, str]:
    """
    从 FASTA 文件中提取特定区域，并记录 ID 映射关系

    参数:
        input_file: 输入 FASTA 文件
        output_fasta: 输出的子区域 FASTA 文件
        start_pos: 起始位置 (1-based, 包含)
        end_pos: 结束位置 (1-based, 包含)

    返回:
        sub_to_orig: 子序列ID -> 原始序列ID 的字典
    """
    sub_to_orig = {}
    sub_records = []
    ndx = 0
    
    print(f"检测到FASTA文件: {input_file}")
    print(f"Detected FASTA file: {input_file}")
    for record in SeqIO.parse(input_file, "fasta"):
        orig_id = record.id
        #print('orig_id:', orig_id)
        # 创建子序列ID
        ndx += 1
        sub_id = f"{ndx}"
        sub_to_orig[sub_id] = orig_id

        # 提取子序列 (转换为0-based索引)
        start_idx = max(0, start_pos - 1)
        end_idx = min(len(record.seq), end_pos)

        if start_idx >= end_idx:
            print(f"警告: 序列 {orig_id} 长度 {len(record.seq)} 小于指定区域，跳过", file=sys.stderr)
            print(f"Warning: Sequence {orig_id} length {len(record.seq)} is less than specified region, skipping", file=sys.stderr)
            continue

        sub_seq = Seq(str(record.seq[start_idx:end_idx]))

        # 创建新记录
        sub_record = SeqRecord(
            seq=sub_seq,
            id=sub_id,
            description=f""
        )
        sub_records.append(sub_record)

    # 写入子序列FASTA
    with open(output_fasta, 'w') as f:
        SeqIO.write(sub_records, f, 'fasta')

    print(f"提取完成: {len(sub_records)} 条序列 -> {output_fasta}")
    print(f"Extraction completed: {len(sub_records)} sequences -> {output_fasta}")
    return sub_to_orig


def run_mmseqs_cluster(
        input_fasta: Path,
        output_prefix: Path,
        threads: int = 8,
        min_seq_id: float = 0.5,
        cov_mode: int = 0,
        coverage: float = 0.8,
        mmseqs_path: str = "mmseqs"
) -> Path:
    """
    运行 MMseqs2 聚类

    参数:
        input_fasta: 输入FASTA文件
        output_prefix: 输出文件前缀
        threads: 线程数
        min_seq_id: 最小序列相似度
        cov_mode: 覆盖度模式
        coverage: 覆盖度阈值
        mmseqs_path: mmseqs 命令路径

    返回:
        cluster_rep: 代表序列FASTA文件路径
    """
    # 确保输出前缀是绝对路径
    output_prefix_abs = output_prefix.resolve()
    input_fasta_abs = input_fasta.resolve()
    
    # 确保输出目录存在
    output_prefix_abs.parent.mkdir(parents=True, exist_ok=True)
    
    # 使用骨架文件夹作为工作目录
    import shutil
    original_cwd = os.getcwd()
    
    try:
        # 切换到输出目录
        os.chdir(output_prefix_abs.parent)
        
        # 执行mmseqs
        cmd = [
            mmseqs_path, "easy-cluster",
            str(input_fasta_abs),  # 使用输入文件的绝对路径
            output_prefix_abs.name,     # 使用输出前缀的名称
            "tmp_mmseqs",               # mmseqs临时文件夹
            "--threads", f"{threads}",
            "--min-seq-id", f"{min_seq_id}",
            "--cov-mode", f"{cov_mode}",
            "-c", f"{coverage}",
        ]

        print(f"运行 MMseqs2: {' '.join(cmd)}")
        print(f"Running MMseqs2: {' '.join(cmd)}")
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        
        # 检查结果文件
        cluster_rep = output_prefix_abs.parent / f"{output_prefix_abs.name}_rep_seq.fasta"
        if not cluster_rep.exists():
            print(f"错误: 代表序列文件 {cluster_rep} 未生成", file=sys.stderr)
            print(f"Error: Representative sequence file {cluster_rep} not generated", file=sys.stderr)
            sys.exit(1)
        
        return cluster_rep
        
    except Exception as e:
        print(f"MMseqs2 执行失败: {e}", file=sys.stderr)
        print(f"MMseqs2 execution failed: {e}", file=sys.stderr)
        sys.exit(1)
    finally:
        # 恢复原始工作目录
        os.chdir(original_cwd)


def output_representative_sequences(
        orig_fasta: Path,
        cluster_rep: Path,
        sub_to_orig: Dict[str, str],
        output_fasta: Path
) -> List[str]:
    """
    输出代表序列的原始完整序列 FASTA

    参数:
        orig_fasta: 原始完整序列 FASTA
        cluster_rep: 聚类代表序列FASTA
        sub_to_orig: 子序列ID -> 原始序列ID 的字典
        output_fasta: 输出的代表序列 FASTA 文件

    返回:
        representative_ids: 代表序列的ID列表
    """
    print('sub_to_orig', sub_to_orig)
    result_records = []
    representative_ids = []
    # 加载原始序列到字典
    orig_records = {record.id: record for record in SeqIO.parse(orig_fasta, "fasta")}
    rep_id_l = [record.id for record in SeqIO.parse(cluster_rep, "fasta")]
    #print('orig_records', orig_records)
    #print('rep_id_l', rep_id_l)
    
    with open(output_fasta, 'w') as f:
        f.truncate(0)
        for rep_id in rep_id_l:
            result_id = sub_to_orig[rep_id]
            result_record = orig_records[result_id]
            result_records.append(result_record)
            representative_ids.append(result_id)
            
            f.write(f'>{result_id}\n')
            f.write(f'{result_record.seq}\n')

    print(f"代表序列输出完成: {len(representative_ids)} 条序列 -> {output_fasta}")
    print(f"Representative sequences output completed: {len(representative_ids)} sequences -> {output_fasta}")
    return representative_ids


def generate_mmseqs_report(
        input_folder: Path,
        output_folder: Path,
        representative_ids_dict: Dict[str, List[str]]
) -> Path:
    """
    生成mmseqs_report.csv报告文件

    参数:
        input_folder: 输入FASTA文件所在文件夹
        output_folder: 输出报告文件夹
        representative_ids_dict: 每个FASTA文件对应的代表序列ID列表

    返回:
        report_path: 报告文件路径
    """
    report_data = []
    
    # 尝试读取mpnn_report.csv
    mpnn_report_path = os.path.join(output_folder, 'mpnn_report.csv')
    print('mpnn_report_path:', mpnn_report_path)
    mpnn_data = None
    
    if os.path.exists(mpnn_report_path):
        print(f"检测到mpnn_report.csv，将继承其中的序列信息")
        print(f"Detected mpnn_report.csv, will inherit sequence information from it")
        try:
            mpnn_data = pd.read_csv(mpnn_report_path)
            print(f"mpnn_report.csv包含 {len(mpnn_data)} 条记录")
            print(f"mpnn_report.csv contains {len(mpnn_data)} records")
        except Exception as e:
            print(f"读取mpnn_report.csv时出错: {e}")
            print(f"Error reading mpnn_report.csv: {e}")
            mpnn_data = None
    else:
        print("未检测到mpnn_report.csv，将使用FASTA文件中的序列信息")
        print("mpnn_report.csv not detected, will use sequence information from FASTA files")
    
    if mpnn_data is not None:
        # 使用mpnn_report.csv中的数据
        for i, row in mpnn_data.iterrows():
            seq_id = row['index']
            
            # 只选择第一行（原始蛋白）或whether_pass为True的行
            if i == 0 or (row.get('whether_pass', False) is True or row.get('whether_pass', 'False') == 'True'):
                # 确定序列所属的backbone
                backbone = row.get('backbone', '')
                if not backbone:
                    # 如果没有backbone信息，尝试从index中提取
                    if '_' in seq_id and not seq_id.startswith('original_protein'):
                        backbone = seq_id.split('_mpnn_')[0]
                
                # 获取当前backbone的代表序列ID列表
                representative_ids = representative_ids_dict.get(backbone, [])
                
                # 确定whether_pass值
                if i == 0 or row.get('whether_pass', '') == '-':
                    # 原始蛋白的whether_pass值设为'-'
                    whether_pass = "-"
                else:
                    # 其他序列的whether_pass值根据是否为代表序列确定
                    whether_pass = "True" if seq_id in representative_ids else "False"
                
                # 创建报告条目，继承mpnn_report.csv中的所有信息（除了whether_pass）
                report_entry = {}
                
                # 复制mpnn_report.csv中的所有列（除了whether_pass）
                for col in mpnn_data.columns:
                    if col != 'whether_pass':
                        report_entry[col] = row.get(col, '')
                
                # 在最后一列添加whether_pass
                report_entry['whether_pass'] = whether_pass
                
                report_data.append(report_entry)
    else:
        # 使用FASTA文件中的数据
        for fa_file in os.listdir(input_folder):
            if fa_file.endswith('.fa') or fa_file.endswith('.fasta'):
                fa_path = os.path.join(input_folder, fa_file)
                base_name = os.path.splitext(fa_file)[0]
                
                # 加载所有序列
                all_sequences = {record.id: str(record.seq) for record in SeqIO.parse(fa_path, "fasta")}
                
                # 获取当前文件的代表序列ID列表
                representative_ids = representative_ids_dict.get(base_name, [])
                
                # 为每个序列生成报告条目
                for seq_id, sequence in all_sequences.items():
                    # 假设原始蛋白是第一个序列，或者ID中不包含'_mpnn_'或'_rfdiffusion_'
                    if '_mpnn_' not in seq_id and '_rfdiffusion_' not in seq_id:
                        whether_pass = "-"
                    else:
                        whether_pass = "True" if seq_id in representative_ids else "False"
                    
                    report_entry = {
                        'index': seq_id,
                        'backbone': base_name,
                        'sequence': sequence,
                        'whether_pass': whether_pass
                    }
                    report_data.append(report_entry)
    
    # 生成报告文件
    report_path = os.path.join(output_folder, 'mmseqs_report.csv')
    df = pd.DataFrame(report_data)
    df.to_csv(report_path, index=False)
    
    print(f"MMseqs报告已生成: {report_path}")
    print(f"共包含 {len(report_data)} 条记录")
    print(f"MMseqs report generated: {report_path}")
    print(f"Total records: {len(report_data)}")
    
    return report_path


def process_fasta_files(
        input_folder: Path,
        output_folder: Path,
        start: int,
        end: int,
        threads: int = 8,
        min_seq_id: float = 0.8,
        cov_mode: int = 0,
        coverage: float = 0.8,
        mmseqs_path: str = "mmseqs"
) -> Tuple[List[Path], Dict[str, List[str]]]:
    """
    处理输入文件夹中的所有FASTA文件，进行聚类分析

    参数:
        input_folder: 输入FASTA文件所在文件夹
        output_folder: 输出文件夹
        start: 起始位置
        end: 结束位置
        threads: 线程数
        min_seq_id: 最小序列相似度
        cov_mode: 覆盖度模式
        coverage: 覆盖度阈值
        mmseqs_path: mmseqs命令路径

    返回:
        generated_files: 生成的文件路径列表
        representative_ids_dict: 每个FASTA文件对应的代表序列ID列表
    """
    # 创建结果文件夹
    results_folder = os.path.join(output_folder, 'results')
    os.makedirs(results_folder, exist_ok=True)
    
    # 创建cluster_data文件夹
    cluster_data_folder = os.path.join(output_folder, 'cluster_data')
    os.makedirs(cluster_data_folder, exist_ok=True)
    
    generated_files = []
    representative_ids_dict = {}
    
    # 遍历所有FASTA文件
    for fa_file in os.listdir(input_folder):
        if fa_file.endswith('.fa') or fa_file.endswith('.fasta'):
            fa_path = os.path.join(input_folder, fa_file)
            base_name = os.path.splitext(fa_file)[0]
            
            print(f"\n处理文件: {fa_file}")
            print(f"\nProcessing file: {fa_file}")
            
            # 创建骨架特定的子文件夹
            skeleton_folder = os.path.join(cluster_data_folder, base_name)
            os.makedirs(skeleton_folder, exist_ok=True)
            
            # 提取子区域
            subregion_fasta = os.path.join(skeleton_folder, "subregion_sequences.fasta")
            sub_to_orig = extract_subregions(
                input_file=Path(fa_path),
                output_fasta=Path(subregion_fasta),
                start_pos=start,
                end_pos=end
            )
            
            if not sub_to_orig:
                print(f"文件 {fa_file} 中没有有效序列，跳过")
                print(f"No valid sequences in file {fa_file}, skipping")
                continue
            
            # 运行聚类
            cluster_prefix = os.path.join(skeleton_folder, "cluster_output")
            cluster_rep = run_mmseqs_cluster(
                input_fasta=Path(subregion_fasta),
                output_prefix=Path(cluster_prefix),
                threads=threads,
                min_seq_id=min_seq_id,
                cov_mode=cov_mode,
                coverage=coverage,
                mmseqs_path=mmseqs_path
            )
            
            # 输出代表序列
            output_fasta = os.path.join(results_folder, f"{base_name}.fa")
            representative_ids = output_representative_sequences(
                orig_fasta=Path(fa_path),
                cluster_rep=Path(cluster_rep),
                sub_to_orig=sub_to_orig,
                output_fasta=Path(output_fasta)
            )
            
            # 记录生成的文件和代表序列ID
            generated_files.append(output_fasta)
            representative_ids_dict[base_name] = representative_ids
    
    return generated_files, representative_ids_dict

def get_start_end(input_str):
    """
    提取输入中的开始数字和结束数字
    :param input_str: 输入字符串（格式：A1-3 / A6 / 1 2 3 等）
    :return: (chain_id, start_num, end_num) 链，开始数字和结束数字
    """
    # 情况1：空格分隔的连续数字（如1 2 3）
    if " " in input_str:
        num_list = [int(num) for num in input_str.split()]  # 按空格分割（支持多空格）
        return 'A', num_list[0], num_list[-1]

    # 情况2：连字符分隔格式（如A1-3、1-5）
    elif "-" in input_str:
        match = re.match(r"^([A-Za-z])*(\d+)-(\d+)$", input_str)
        if match:
            if not match.group(1):
                chain_id = 'A'
            else:
                chain_id = match.group(1)
            start = int(match.group(2))
            end = int(match.group(3))
            return chain_id, start, end

    # 情况3：字母+数字（如A6）或纯数字（如6）
    else:
        match = re.match(r"^([A-Za-z])*(\d+)$", input_str)
        if match:
            if not match.group(1):
                chain_id = 'A'
            else:
                chain_id = match.group(1)
            num = int(match.group(2))
            return chain_id, num, num

    # 无效输入返回None（可选）
    return None, None, None


def main():
    args = arg_parser()
    input_folder = os.path.expanduser(args.input_folder)
    output_folder = os.path.expanduser(args.output_folder)
    threads = args.threads
    min_seq_id = args.min_seq_id
    cov_mode = args.cov_mode
    coverage = args.coverage
    mmseqs_path = args.mmseqs_path

    position_list = args.position_list

    _, start, end = get_start_end(position_list)

    work_dir = output_folder.rsplit("/",1)[0]

    
    print("=== MMseqs2 聚类分析工具 ===")
    print("=== MMseqs2 Clustering Analysis Tool ===")
    print(f"输入文件夹: {input_folder}")
    print(f"Input folder: {input_folder}")
    print(f"输出文件夹: {output_folder}")
    print(f"Output folder: {output_folder}")
    print(f"分析区域: {start}-{end}")
    print(f"Analysis region: {start}-{end}")
    print(f"线程数: {threads}")
    print(f"Number of threads: {threads}")
    print(f"最小序列相似度: {min_seq_id}")
    print(f"Minimum sequence identity: {min_seq_id}")
    print(f"覆盖度模式: {cov_mode}")
    print(f"Coverage mode: {cov_mode}")
    print(f"覆盖度阈值: {coverage}")
    print(f"Coverage threshold: {coverage}")
    print(f"MMseqs路径: {mmseqs_path}")
    print(f"MMseqs path: {mmseqs_path}")
    
    # 处理所有FASTA文件
    generated_files, representative_ids_dict = process_fasta_files(
        input_folder=Path(input_folder),
        output_folder=Path(output_folder),
        start=start,
        end=end,
        threads=threads,
        min_seq_id=min_seq_id,
        cov_mode=cov_mode,
        coverage=coverage,
        mmseqs_path=mmseqs_path
    )
    
    if generated_files:
        print(f"\n成功处理 {len(generated_files)} 个FASTA文件")
        print(f"\nSuccessfully processed {len(generated_files)} FASTA files")
        
        # 生成mmseqs_report.csv报告
        report_path = generate_mmseqs_report(
            input_folder=Path(input_folder),
            output_folder=Path(work_dir),
            representative_ids_dict=representative_ids_dict
        )
        
        print(f"\n[SUCCESS] MMseqs聚类分析完成！")
        print(f"\n[SUCCESS] MMseqs clustering analysis completed!")
        print(f"[OUTPUT] 主要输出文件:")
        print(f"[OUTPUT] Main output files:")
        print(f"   - 聚类代表序列: {output_folder}/results/")
        print(f"   - Cluster representative sequences: {output_folder}/results/")
        print(f"   - 聚类分析数据: {output_folder}/cluster_data/")
        print(f"   - Cluster analysis data: {output_folder}/cluster_data/")
        print(f"   - MMseqs报告: {report_path}")
        print(f"   - MMseqs report: {report_path}")
    else:
        print("[ERROR] 没有成功处理任何文件")
        print("[ERROR] No files were successfully processed")
        sys.exit(1)


if __name__ == "__main__":
    main()
