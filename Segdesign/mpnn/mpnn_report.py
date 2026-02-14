import threading
import sys
import subprocess
import shlex
import argparse
import re
from pathlib import Path
import os
import pandas as pd
import shutil
import math
from Bio import SeqIO
import csv
from protein_data import process_protein_data




def parse_args():
    parser = argparse.ArgumentParser(description='Protein sequence prediction and report generation')
    parser.add_argument('--protein_pdb', type=str, default=None,
                        help='Path to the original protein PDB file')
    parser.add_argument("--seq_folder", type=str,
                        help="Folder containing MPNN generated fasta files")
    parser.add_argument("--output_folder", type=str,
                        help="Folder for storing output files")
    parser.add_argument("--final_report_folder", type=str, default=None,
                        help="Folder for storing final mpnn_report.csv (default: same as output_folder)")
    parser.add_argument('--top_percent', type=float, default=0.2,
                        help='Filter sequences with the lowest global_score by percentage')
    parser.add_argument("--generate_report", type=bool, default=True,
                        help="Generate comprehensive MPNN report")
    parser.add_argument("--rfdiffusion_report_path", type=str, default=None,
                        help="The path to rfdiffusion_report.csv. If not entered, the default path will be used: {{work_dir}}/rfdiffusion_report.csv")
    parser.add_argument('--position_list', type=str, default=None,
                        help='Specified design area, such as A1-5')
    
    return parser.parse_args()


def extract_sequences_from_fasta(file_path):
    """
    从FASTA文件中提取序列数据，区分初始序列和生成序列
    """
    sequences = []
    try:
        with open(file_path, 'r') as f:
            for record in SeqIO.parse(f, "fasta"):
                # 解析头部信息
                header = record.description
                
                # 提取属性
                attributes = {}
                header_parts = header.split(', ')
                for part in header_parts:
                    if '=' in part:
                        key, value = part.split('=', 1)
                        attributes[key.strip()] = value.strip()
                
                sequence_data = {
                    'header': header,
                    'attributes': attributes,
                    'sequence': str(record.seq),
                    'id': record.id
                }
                sequences.append(sequence_data)
    except Exception as e:
        print(f"读取文件 {file_path} 时出错：{e}")
        print(f"Error reading file {file_path}: {e}")
        return []
    
    return sequences


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


def load_backbone_data_from_rfdiffusion(working_dir, rfdiffusion_report_path = None):
    """
    从rfdiffusion_report.csv中加载骨架数据

    
    """
    backbone_data = {}
    try:
        # 构建rfdiffusion_report.csv的完整路径
        if rfdiffusion_report_path == 'None' or rfdiffusion_report_path is None:
            rf_report_path = os.path.join(working_dir, 'rfdiffusion_report.csv')
        else:
            rf_report_path = rfdiffusion_report_path
        #if not os.path.exists(rf_report_path):
            # 尝试其他可能的路径
           # rf_report_path = os.path.join(working_dir, 'rfdiffusion_out', 'rfdiffusion_report.csv')
        
        if os.path.exists(rf_report_path):
            df_rf = pd.read_csv(rf_report_path)
            for _, row in df_rf.iterrows():
                backbone_index = row['index']
                backbone_data[backbone_index] = {
                    'ss8': row.get('design_ss8', ''),
                    'ss3': row.get('design_ss3', ''),
                    'H_prop': row.get('H_prop', 0.0),
                    'E_prop': row.get('E_prop', 0.0),
                    'C_prop': row.get('C_prop', 0.0),
                    'backbone': row.get('backbone', ''),
                    'success_backbone': row.get('success_backbone', ''),
                    'Success': row.get('Success', '')
                }
            print(f"已加载 {len(backbone_data)} 个骨架的数据")
            print(f"Loaded data for {len(backbone_data)} backbones")
        else:
            print(f"警告：找不到rfdiffusion_report.csv文件: {rf_report_path}")
            print(f"Warning: rfdiffusion_report.csv file not found: {rf_report_path}")
            
    except Exception as e:
        print(f"读取rfdiffusion_report.csv时出错：{e}")
        print(f"Error reading rfdiffusion_report.csv: {e}")
    
    return backbone_data




def generate_csv_for_fasta(seq_file_path, output_folder, fa_filename, position_list, working_dir, rfdiffusion_report_path = None):
    """
    为单个FASTA文件生成CSV文件，包含完整的骨架信息和MPNN数据
    """
    print(f"处理文件：{fa_filename}")
    print(f"Processing file: {fa_filename}")
    
    # 提取所有序列
    sequences = extract_sequences_from_fasta(seq_file_path)
    
    if not sequences:
        print(f"文件 {fa_filename} 中没有找到有效序列")
        print(f"No valid sequences found in file {fa_filename}")
        return None
    
    # 第一个序列是初始序列，后续是生成序列
    generated_sequences = sequences[1:] if len(sequences) > 1 else []
    
    if not generated_sequences:
        print(f"文件 {fa_filename} 中没有找到生成序列")
        print(f"No generated sequences found in file {fa_filename}")
        return None
    
    # 从rfdiffusion_report.csv加载骨架数据
    backbone_data = load_backbone_data_from_rfdiffusion(working_dir, rfdiffusion_report_path)
    #print(f"backbone data: {backbone_data}")
    
    # 提取骨架ID（从文件名如"Dusp4_A_2"）
    backbone_id = fa_filename.replace('.fa', '')
    
    # 获取设计区域位置
    _, design_start, design_end = get_start_end(position_list)
    
    # 获取对应的骨架数据
    backbone_info = backbone_data.get(backbone_id, {
        'segment':'',
        'ss8': '',
        'ss3': '',
        'H_prop': 0.0,
        'E_prop': 0.0,
        'C_prop': 0.0,
        'backbone': '',
        'success_backbone': '',
        'Success': ''
    })
    print(f"backbone info: {backbone_info}")
    
    # 准备CSV数据
    csv_data = []
    
    for idx, seq_data in enumerate(generated_sequences):
        # 提取MPNN属性
        score = float(seq_data['attributes'].get('score', '0.0'))
        global_score = float(seq_data['attributes'].get('global_score', '0.0'))
        
        # 计算设计区域序列（从设计区域位置提取）
        full_sequence = seq_data['sequence']
        if len(full_sequence) >= design_end:
            design_region = full_sequence[design_start-1:design_end]  # Python索引从0开始
        else:
            design_region = full_sequence
        
        csv_row = {
            'index': f"{backbone_id}_mpnn_{idx}",
            'backbone': backbone_id,
            'segment': backbone_info.get('segment', f'{design_start}-{design_end}'),
            'ss8': backbone_info['ss8'],
            'ss3': backbone_info['ss3'],
            'H_prop': backbone_info['H_prop'],
            'E_prop': backbone_info['E_prop'],
            'C_prop': backbone_info['C_prop'],
            'backbone_pdb': backbone_info['success_backbone'] if backbone_info['success_backbone'] != '-' else backbone_info['backbone'],
            'score': score,
            'global_score': global_score,
            'region': design_region,
            'sequence': full_sequence
        }
        csv_data.append(csv_row)
    
    # 生成CSV文件名
    csv_filename = f"mpnn_{backbone_id}.csv"
    csv_path = os.path.join(output_folder, csv_filename)
    
    # 保存CSV文件
    df = pd.DataFrame(csv_data)
    df.to_csv(csv_path, index=False)
    
    print(f"已生成CSV文件：{csv_filename}，包含 {len(csv_data)} 个序列")
    print(f"CSV file generated: {csv_filename}, containing {len(csv_data)} sequences")
    return csv_path, csv_data


def filter_top_sequences(csv_data, top_percent):
    """
    根据global_score筛选最低的top_percent百分比序列（保持原始index顺序）
    """
    if not csv_data:
        return []
    
    # 按global_score排序（升序，数值越低越好）以找出需要保留的序列
    sorted_data = sorted(csv_data, key=lambda x: x['global_score'])
    
    # 计算需要保留的序列数量
    total_sequences = len(sorted_data)
    n = max(1, math.ceil(total_sequences * top_percent))
    
    # 获取需要保留的序列的index（按原始顺序）
    top_indices = {seq['index'] for seq in sorted_data[:n]}
    
    # 按原始顺序返回筛选后的序列
    filtered_sequences = [seq for seq in csv_data if seq['index'] in top_indices]
    
    return filtered_sequences


def process_all_fasta_files(seq_folder, output_folder, top_percent, position_list, rfdiffusion_report_path = None):
    """
    处理所有FASTA文件并生成相应的CSV文件
    """
    print(f"开始处理FASTA文件...")
    print(f"Starting to process FASTA files...")
    print(f"输入文件夹：{seq_folder}")
    print(f"Input folder: {seq_folder}")
    print(f"输出文件夹：{output_folder}")
    print(f"Output folder: {output_folder}")
    
    # 获取工作目录（假设seq_folder在工作目录下的mpnn_out/seqs）
    working_dir = output_folder.rsplit('/', 1)[0]
    print(f"工作目录：{working_dir}")
    print(f"Working directory: {working_dir}")
    
    # 创建输出文件夹
    os.makedirs(output_folder, exist_ok=True)
    
    # 获取所有FASTA文件
    fa_files = sorted([f for f in os.listdir(seq_folder) if f.endswith('.fa')], 
                     key=natural_sort_key)
    #print('fa_files:', fa_files)
    if not fa_files:
        print(f"在文件夹 {seq_folder} 中没有找到FASTA文件")
        print(f"No FASTA files found in folder {seq_folder}")
        return [], []
    
    print(f"找到 {len(fa_files)} 个FASTA文件")
    print(f"Found {len(fa_files)} FASTA files")
    
    # 创建seqs_csv文件夹
    seqs_csv_folder = os.path.join(output_folder, 'seqs_csv')
    os.makedirs(seqs_csv_folder, exist_ok=True)
    
    # 处理每个FASTA文件
    all_csv_data = []
    generated_files = []
    
    # 创建top_percent文件夹
    top_percent_str = f"{top_percent*100:.1f}%"
    top_folder = os.path.join(output_folder, f'top_{top_percent_str}')
    os.makedirs(top_folder, exist_ok=True)
    
    # 创建top_filter文件夹用于保存fasta文件
    top_filter_folder = os.path.join(output_folder, 'top_filter')
    os.makedirs(top_filter_folder, exist_ok=True)
    
    # 对每个FASTA文件独立进行筛选
    top_generated_files = []
    top_fasta_files = []
    
    for fa_file in fa_files:
        fa_file_path = os.path.join(seq_folder, fa_file)
        
        result = generate_csv_for_fasta(fa_file_path, seqs_csv_folder, fa_file, position_list, working_dir, rfdiffusion_report_path)
        if result:
            csv_path, csv_data = result
            generated_files.append(csv_path)
            all_csv_data.extend(csv_data)
            
            # 对当前文件的序列进行独立筛选
            top_sequences_current = filter_top_sequences(csv_data, top_percent)
            
            if top_sequences_current:
                base_name = os.path.splitext(fa_file)[0]
                top_csv_filename = f"top_mpnn_{base_name}.csv"
                top_csv_path = os.path.join(top_folder, top_csv_filename)
                
                df = pd.DataFrame(top_sequences_current)
                df.to_csv(top_csv_path, index=False)
                top_generated_files.append(top_csv_path)
                
                print(f"已生成Top序列CSV文件：{top_csv_filename}，包含 {len(top_sequences_current)} 个序列")
                print(f"Top sequence CSV file generated: {top_csv_filename}, containing {len(top_sequences_current)} sequences")
                
                # 将top序列保存为fasta文件
                base_name = os.path.splitext(fa_file)[0]
                fasta_filename = f"{base_name}.fa"
                fasta_path = os.path.join(top_filter_folder, fasta_filename)
                
                with open(fasta_path, 'w') as f:
                    for seq_data in top_sequences_current:
                        # 创建fasta头部
                        header = f">{seq_data['index']} score={seq_data['score']} global_score={seq_data['global_score']}"
                        # 写入头部和序列
                        f.write(f"{header}\n{seq_data['sequence']}\n")
                
                print(f"已生成Top序列FASTA文件：{fasta_filename}，包含 {len(top_sequences_current)} 个序列")
                print(f"Top sequence FASTA file generated: {fasta_filename}, containing {len(top_sequences_current)} sequences")
                top_fasta_files.append(fasta_path)
            else:
                print(f"文件 {fa_file} 中没有符合筛选条件的序列")
                print(f"No sequences in file {fa_file} meet the filtering criteria")
    
    print(f"Top序列CSV文件已保存到文件夹：{top_folder}")
    print(f"Top sequence CSV files saved to folder: {top_folder}")
    print(f"Top序列FASTA文件已保存到文件夹：{top_filter_folder}")
    print(f"Top sequence FASTA files saved to folder: {top_filter_folder}")
    
    return generated_files, top_generated_files, top_fasta_files


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


def generate_final_mpnn_report(output_folder, top_percent, protein_pdb, position_list, rfdiffusion_report_path=None, final_report_folder=None):
    """
    生成最终的mpnn_report.csv文件，包含所有序列
    
    参数:
        output_folder: 输出文件夹
        top_percent: top筛选百分比
        rfdiffusion_report_path: rfdiffusion报告路径（默认为None）
        final_report_folder: 最终报告输出文件夹（默认为output_folder）
    """
    print("生成最终的MPNN报告（包含所有序列）...")
    print("Generating final MPNN report (including all sequences)...")
    segment = ""

    # 确定最终报告输出路径
    if final_report_folder is None:
        final_report_folder = output_folder.rsplit('/', 1)[0]
    
    # 只读取所有原始序列CSV文件
    seqs_csv_folder = os.path.join(output_folder, 'seqs_csv')
    top_percent_str = f"{top_percent*100:.1f}%"
    top_folder = os.path.join(output_folder, f'top_{top_percent_str}')
    
    # 获取所有序列的index集合（用于标记是否为Top序列）
    top_sequence_indices = set()
    if os.path.exists(top_folder):
        top_csv_files = sorted([f for f in os.listdir(top_folder) if f.endswith('.csv')],
                               key=natural_sort_key)
        for csv_file in top_csv_files:
            csv_path = os.path.join(top_folder, csv_file)
            df_top = pd.read_csv(csv_path)
            top_sequence_indices.update(df_top['index'].tolist())
    
    report_data = []
    
    # 读取rfdiffusion_report.csv的第一行并添加到report_data开头
    try:
        if rfdiffusion_report_path is None or rfdiffusion_report_path == 'None':
            rfdiffusion_report_path = os.path.join(final_report_folder, 'rfdiffusion_report.csv')
        if os.path.exists(rfdiffusion_report_path):
            df_rf = pd.read_csv(rfdiffusion_report_path)
            if not df_rf.empty:
                # 获取第一行数据
                rf_first_row = df_rf.iloc[0].to_dict()

                chain_id, start_num, end_num = get_start_end(position_list)
                # 创建原始蛋白信息行
                original_protein_row = process_protein_data(
                    pdb_path=protein_pdb,
                    chain_id=chain_id,
                    segment_range=f'{start_num}-{end_num}',
                    output_path=os.path.join(output_folder, 'original_protein_data')
                )
                
                # 添加到report_data开头
                report_data.append(original_protein_row)
                print("已添加原始蛋白信息到报告开头")
                print("Original protein information added to the beginning of the report")
    except Exception as e:
        print(f"读取rfdiffusion_report.csv时出错：{e}")
        print(f"Error reading rfdiffusion_report.csv: {e}")
    
    # 处理所有原始序列CSV文件
    if os.path.exists(seqs_csv_folder):
        csv_files = sorted([f for f in os.listdir(seqs_csv_folder) if f.endswith('.csv')],
                           key=natural_sort_key)
        
        for csv_file in csv_files:
            csv_path = os.path.join(seqs_csv_folder, csv_file)
            df = pd.read_csv(csv_path)
            
            for _, row in df.iterrows():
                # 检查该序列是否在Top筛选中
                is_top_sequence = row['index'] in top_sequence_indices
                
                report_entry = {
                    'index': row['index'],
                    'backbone': row.get('backbone', ''),
                    'segment': row.get('segment', ''),
                    'ss8': row.get('ss8', ''),
                    'ss3': row.get('ss3', ''),
                    'H_prop': row.get('H_prop', ''),
                    'E_prop': row.get('E_prop', ''),
                    'C_prop': row.get('C_prop', ''),
                    'backbone_pdb': row.get('backbone_pdb', ''),
                    'score': row['score'],
                    'global_score': row['global_score'],
                    'region': row.get('region', ''),
                    'sequence': row['sequence'],
                    'whether_pass': 'True' if is_top_sequence else 'False'
                }
                report_data.append(report_entry)
    
    # 生成最终报告
    if report_data:
        final_report_path = os.path.join(final_report_folder, 'mpnn_report.csv')
        df_final = pd.DataFrame(report_data)
        df_final.to_csv(final_report_path, index=False)
        
        print(f"最终MPNN报告已生成：{final_report_path}")
        print(f"包含 {len(report_data)} 条记录")
        print(f"Final MPNN report generated: {final_report_path}")
        print(f"Contains {len(report_data)} records")
        

        return final_report_path
    else:
        print("没有数据生成最终报告")
        print("No data to generate final report")
        return None






if __name__ == "__main__":
    args = parse_args()
    seq_folder = os.path.expanduser(args.seq_folder)
    output_folder = os.path.expanduser(args.output_folder)
    top_percent = args.top_percent
    rfdiffusion_report_path = args.rfdiffusion_report_path
    position_list = args.position_list
    protein_pdb = args.protein_pdb
    
    print("=== MPNN序列处理和报告生成 ===")
    print("=== MPNN Sequence Processing and Report Generation ===")
    print(f"输入序列文件夹: {seq_folder}")
    print(f"Input sequence folder: {seq_folder}")
    print(f"输出文件夹: {output_folder}")
    print(f"Output folder: {output_folder}")
    print(f"Top筛选百分比: {top_percent*100:.1f}%")
    print(f"Top filtering percentage: {top_percent*100:.1f}%")
    
    # 处理所有FASTA文件并生成CSV
    all_csv_files, top_csv_files, top_fasta_files = process_all_fasta_files(seq_folder, output_folder, top_percent, position_list, rfdiffusion_report_path)
    
    if all_csv_files:
        print(f"成功处理 {len(all_csv_files)} 个CSV文件")
        print(f"Successfully processed {len(all_csv_files)} CSV files")
        if top_csv_files:
            print(f"成功生成 {len(top_csv_files)} 个Top序列CSV文件")
            print(f"Successfully generated {len(top_csv_files)} top sequence CSV files")
        if top_fasta_files:
            print(f"成功生成 {len(top_fasta_files)} 个Top序列FASTA文件")
            print(f"Successfully generated {len(top_fasta_files)} top sequence FASTA files")



        # 生成最终报告
        if args.generate_report:
            final_report_folder = args.final_report_folder
            final_report_path = generate_final_mpnn_report(
                output_folder,
                top_percent,
                protein_pdb,
                position_list,
                rfdiffusion_report_path,
                final_report_folder
            )
            if final_report_path:
                print(f"\n[SUCCESS] 完整MPNN报告生成完成！")
                print(f"\n[SUCCESS] Complete MPNN report generation completed!")
                print(f"[OUTPUT] 主要输出文件:")
                print(f"[OUTPUT] Main output files:")
                print(f"   - 原始序列CSV: {output_folder}/seqs_csv/")
                print(f"   - Original sequence CSV: {output_folder}/seqs_csv/")
                print(f"   - Top序列CSV: {output_folder}/top_{top_percent * 100:.1f}%/")
                print(f"   - Top sequence CSV: {output_folder}/top_{top_percent * 100:.1f}%/")
                print(f"   - Top序列FASTA: {output_folder}/top_filter/")
                print(f"   - Top sequence FASTA: {output_folder}/top_filter/")
                print(f"   - 最终报告: {final_report_path}")
                print(f"   - Final report: {final_report_path}")
    else:
        print("[ERROR] 没有成功处理任何文件")
        print("[ERROR] No files were successfully processed")
        sys.exit(1)