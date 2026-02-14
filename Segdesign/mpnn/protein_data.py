#!/usr/bin/env python3
"""
蛋白质数据处理工具
功能：从PDB文件中提取蛋白质序列和二级结构信息，生成FASTA文件和DSSP文件，并返回结构化数据
"""

import os
import re
import sys
from pathlib import Path
import os
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))  # rfdiffusion/目录
PARENT_DIR = os.path.dirname(ROOT_DIR)  # 上级目录（包含dssp/和rfdiffusion/）
sys.path.append(PARENT_DIR)  # 把上级目录加入搜索路径
# 导入dssp模块
from dssp.dssp import run_dssp
from dssp.dsspcsv import parse_dssp


def extract_sequence_from_pdb(pdb_path, chain_id='A'):
    """
    从PDB文件中提取指定链的蛋白质序列

    参数:
        pdb_path: PDB文件路径
        chain_id: 链ID，默认为'A'

    返回:
        完整的蛋白质序列
    """
    sequence = []
    processed_residues = set()

    # 标准和特殊氨基酸的三字母到单字母映射
    aa_dict = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
        # 特殊氨基酸处理
        'SEC': 'U', 'PYL': 'O', 'MSE': 'M', 'SEP': 'S', 'TPO': 'T',
        'PTR': 'Y', 'XLE': 'L', 'GLX': 'Z', 'ASX': 'B'
    }

    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                # 提取链ID和残基信息
                chain = line[21].strip()
                if chain != chain_id:
                    continue

                # 提取残基名称、序号和插入码
                residue_name = line[17:20].strip()
                residue_num = line[22:26].strip()
                icode = line[26].strip()

                # 创建唯一标识符 (链ID + 残基序号 + 插入码)
                residue_id = (chain, residue_num, icode)

                # 跳过已处理的残基
                if residue_id in processed_residues:
                    continue

                # 确定氨基酸单字母代码
                if residue_name in aa_dict:
                    aa = aa_dict[residue_name]
                else:
                    # 非标准氨基酸用X表示
                    aa = 'X'

                # 添加到序列中
                sequence.append(aa)

                # 标记此残基已处理
                processed_residues.add(residue_id)

    return ''.join(sequence)


def get_segment_sequence(sequence, segment_range):
    """
    从完整序列中提取指定片段的序列

    参数:
        sequence: 完整的蛋白质序列
        segment_range: 片段范围，如'1-5'

    返回:
        指定片段的蛋白质序列
    """
    try:
        start, end = map(int, segment_range.split('-'))
        # 转换为0-based索引
        start_idx = max(0, start - 1)
        end_idx = min(len(sequence), end)
        return sequence[start_idx:end_idx]
    except:
        return ''


def calculate_ss_properties(ss3):
    """
    计算二级结构中H、E、C的占比

    参数:
        ss3: ss3类型的二级结构序列

    返回:
        H_prop: H的占比
        E_prop: E的占比
        C_prop: C的占比
    """
    if not ss3:
        return 0.0, 0.0, 0.0

    total = len(ss3)
    h_count = ss3.count('H')
    e_count = ss3.count('E')
    c_count = ss3.count('C')

    h_prop = round(h_count / total, 4)
    e_prop = round(e_count / total, 4)
    c_prop = round(c_count / total, 4)

    return h_prop, e_prop, c_prop


def generate_fasta(pdb_path, chain_id='A', output_folder=None):
    """
    从PDB文件生成FASTA序列文件

    参数:
        pdb_path: PDB文件路径
        chain_id: 链ID，默认为'A'
        output_folder: 输出文件夹，默认为PDB文件所在目录

    返回:
        fasta文件路径
    """
    if output_folder is None:
        output_folder = os.path.dirname(pdb_path)

    # 提取蛋白质名称（不含后缀）
    protein_name = os.path.splitext(os.path.basename(pdb_path))[0]
    fasta_filename = f"{protein_name}_chain{chain_id}.fa"
    fasta_path = os.path.join(output_folder, fasta_filename)

    # 提取序列
    sequence = extract_sequence_from_pdb(pdb_path, chain_id)

    # 写入FASTA文件
    os.makedirs(output_folder, exist_ok=True)
    with open(fasta_path, 'w') as f:
        f.write(f">{protein_name}_chain{chain_id}\n")
        f.write(f"{sequence}\n")

    return fasta_path


def process_protein_data(pdb_path, chain_id='A', segment_range=None, output_path=None):
    """
    处理蛋白质数据，生成FASTA文件和DSSP文件，并返回结构化数据

    参数:
        pdb_path: 蛋白质PDB文件路径
        chain_id: 链ID，默认为'A'
        segment_range: 片段范围，如'1-5'
        output_path: 输出路径，默认为None（使用PDB文件所在目录的output子目录）

    返回:
        包含蛋白质数据的字典
    """
    # 验证PDB文件是否存在
    if not os.path.exists(pdb_path):
        print(f"错误: PDB文件不存在: {pdb_path}")
        print(f"Error: PDB file does not exist: {pdb_path}")
        return None

    # 提取蛋白质名称（不含后缀）
    protein_name = os.path.splitext(os.path.basename(pdb_path))[0]

    # 确定输出文件夹
    if output_path:
        output_folder = output_path
    else:
        output_folder = os.path.join(os.path.dirname(pdb_path), 'output')
    os.makedirs(output_folder, exist_ok=True)

    # 提取完整序列
    full_sequence = extract_sequence_from_pdb(pdb_path, chain_id)

    # 提取指定片段的序列
    segment_seq = get_segment_sequence(full_sequence, segment_range)

    # 生成FASTA文件
    fasta_path = generate_fasta(pdb_path, chain_id, output_folder)

    # 生成DSSP文件
    dssp_filename = f"{protein_name}_chain{chain_id}.dssp"
    dssp_path = os.path.join(output_folder, dssp_filename)
    run_dssp(pdb_path, dssp_path)

    # 解析DSSP文件获取二级结构信息
    ss8_full = ''
    ss3_full = ''
    
    if os.path.exists(dssp_path):
        with open(dssp_path, 'r') as f:
            dssp_content = f.read()
        
        dssp_data = parse_dssp(dssp_content)
        
        # 提取指定链的二级结构
        for item in dssp_data:
            res_num, chain, _, ss8, ss3, _ = item
            if chain == chain_id:
                ss8_full += ss8
                ss3_full += ss3
    
    # 提取指定片段的二级结构
    try:
        start, end = map(int, segment_range.split('-'))
        start_idx = max(0, start - 1)
        end_idx = min(len(ss8_full), end)
        ss8_segment = ss8_full[start_idx:end_idx]
        ss3_segment = ss3_full[start_idx:end_idx]
    except:
        ss8_segment = ''
        ss3_segment = ''

    # 计算二级结构中H、E、C的占比
    H_prop, E_prop, C_prop = calculate_ss_properties(ss3_segment)

    # 构建返回字典
    protein_data = {
        'index': f'{protein_name}_{chain_id}',
        'backbone': protein_name,
        'segment': segment_range,
        'ss8': ss8_segment,
        'ss3': ss3_segment,
        'H_prop': H_prop,
        'E_prop': E_prop,
        'C_prop': C_prop,
        'backbone_pdb': pdb_path,
        'score': '-',
        'global_score': '-',
        'region': segment_seq,
        'sequence': full_sequence,
        'whether_pass': '-'  
    }

    return protein_data


def main():
    """
    主函数，用于测试
    """
    if len(sys.argv) < 2:
        print("用法: python protein_data.py <pdb_path> [chain_id] [segment_range] [output_path]")
        print("Usage: python protein_data.py <pdb_path> [chain_id] [segment_range] [output_path]")
        print("示例: python protein_data.py example.pdb A 1-50 ./output")
        print("Example: python protein_data.py example.pdb A 1-50 ./output")
        sys.exit(1)

    pdb_path = sys.argv[1]
    chain_id = sys.argv[2] if len(sys.argv) > 2 else 'A'
    segment_range = sys.argv[3] if len(sys.argv) > 3 else '1-10'
    output_path = sys.argv[4] if len(sys.argv) > 4 else None

    print(f"处理PDB文件: {pdb_path}")
    print(f"Processing PDB file: {pdb_path}")
    print(f"链ID: {chain_id}")
    print(f"Chain ID: {chain_id}")
    print(f"片段范围: {segment_range}")
    print(f"Segment range: {segment_range}")
    if output_path:
        print(f"输出路径: {output_path}")
        print(f"Output path: {output_path}")

    # 处理蛋白质数据
    protein_data = process_protein_data(pdb_path, chain_id, segment_range, output_path)

    if protein_data:
        print("\n成功提取蛋白质数据:")
        print("\nSuccessfully extracted protein data:")
        for key, value in protein_data.items():
            print(f"{key}: {value}")
        print("\n处理完成！")
        print("\nProcessing completed!")
    else:
        print("\n处理失败！")
        print("\nProcessing failed!")


if __name__ == "__main__":
    main()
