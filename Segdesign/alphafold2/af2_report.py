import os
import re
import pandas as pd
import shutil
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys
import argparse
root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# 将根目录添加到 Python 的模块搜索路径中
sys.path.append(root_dir)
from dssp.dssp import run_dssp
from dssp.dsspcsv import dssp_to_csv
from hmmer.pdb_to_fasta import pdb_to_fasta
import json
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(description='Protein 3D Structure Prediction(alphafold2)', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--esmfold_report_path', type=str, default=None,
                        help="Folder for storing ESMFold report")
    parser.add_argument('--alphafold2_folder','-i', type=str,
                        help="The folder where alphafold2's output data is stored")
    parser.add_argument('--plddt_threshold', type=float, default=None,
                        help='pLDDT selection threshold')
    parser.add_argument('--ptm_threshold', type=float, default=None,
                        help='ptm selection threshold')
    parser.add_argument('--seq_range_str', type=str,
                        help='Enter the area to be modified, in the format: start position-end position, such as 1-10')
    parser.add_argument('--ss', type=str, default=None,
                        help='Secondary structure type, such as helix and strand')
    parser.add_argument('--ss_threshold', type=float, default=None,
                        help='Threshold for DSSP analysis')
    parser.add_argument('--ss_filter', type=bool, default=True,
                        help='Whether to enable filtering based on secondary structure, default is True')
    return parser.parse_args()


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

#为预测的pdb生成dssp文件和对应的csv文件，返回名称字典
def make_dssp_csv(structure_prediction_files_folder, dssp_folder, csv_folder):
    backbone_folders = os.listdir(structure_prediction_files_folder)
    backbone_folders = sorted(backbone_folders, key=natural_sort_key)
    ndx_dict = {}
    for backbone_folder in backbone_folders:
        backbone_folder_path = os.path.join(structure_prediction_files_folder, backbone_folder)
        ndx_dict[backbone_folder] = {}
        dssp_backbone_folder_path = os.path.join(dssp_folder, backbone_folder)
        if not os.path.exists(dssp_backbone_folder_path):
            os.makedirs(dssp_backbone_folder_path, exist_ok=True)

        csv_backbone_folder_path = os.path.join(csv_folder, backbone_folder)
        if not os.path.exists(csv_backbone_folder_path):
            os.makedirs(csv_backbone_folder_path, exist_ok=True)

        pdb_files = [file for file in os.listdir(backbone_folder_path) if file.endswith(".pdb")]
        pdb_files = sorted(pdb_files, key=natural_sort_key)

        for pdb_file in pdb_files:
            file_name = os.path.splitext(pdb_file)[0]
            ndx_dict[backbone_folder][file_name] = {}
            pdb_file_path = os.path.join(backbone_folder_path, pdb_file)
            dssp_file_path = os.path.join(dssp_backbone_folder_path, f'{file_name}.dssp')
            csv_file_path = os.path.join(csv_backbone_folder_path, f'{file_name}.csv')
            run_dssp(pdb_file_path, dssp_file_path)
            dssp_to_csv(dssp_file_path, csv_file_path)
    
    return ndx_dict

def read_ptm_and_plddt(structure_prediction_files_folder, ndx_dict):
    def read_json(json_file):
        with open(json_file, 'r') as f:
            data = json.load(f)
        plddt_l = data['plddt']
        plddt = np.mean(np.array(plddt_l))
        ptm = data['ptm']
        return ptm, plddt

    score_dict = {}
    for backbone, backbone_ndx_dict in ndx_dict.items():
        score_dict[backbone] = {}
        for file_name, _ in backbone_ndx_dict.items():
            json_file_path = os.path.join(structure_prediction_files_folder, backbone, f'{file_name}.json')
            ptm, plddt = read_json(json_file_path)
            score_dict[backbone][file_name] = {'plddt': round(plddt,2), 'ptm': ptm}

    return score_dict

def add_ptm_and_plddt(score_dict, esmfold_report_path):
    """根据行的backbone和index从字典提取数据（无匹配返回NaN）"""
    df = pd.read_csv(esmfold_report_path)
    # 2. 查看CSV结构（重点确认首行数据，是否为原始蛋白）
    # 2. 验证必要列（whether_pass用于检测，后续会删除）
    required_cols = ["backbone", "index", "whether_pass"]  # 仍需验证whether_pass存在
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"CSV文件缺少必要列：{missing_cols}，请检查列名是否正确！")

    # 3. 查看原始蛋白行（第一行）
    print("=== 原始蛋白行（第一行）关键信息 ===")
    print("=== Key information of original protein row (first row) ===")
    print(df.iloc[0][["backbone", "index", "whether_pass"]])

    

    dict_target_pairs = set()
    for backbone, index_dict in score_dict.items():
        for idx in index_dict.keys():
            dict_target_pairs.add((str(backbone), str(idx)))
    print(f"\n=== 字典目标行总数：{len(dict_target_pairs)} ===")
    print(f"\n=== Total number of dictionary target rows: {len(dict_target_pairs)} ===")

    # 2. 检测字典目标行与whether_pass=True行的对应关系
    df["whether_pass_bool"] = df["whether_pass"].apply(
        lambda x: x if isinstance(x, bool) else (x.lower() == "true")
    )
    pass_rows = df[(df["whether_pass_bool"]) & (df.index != 0)]  # 排除原始蛋白
    pass_target_pairs = set((str(row["backbone"]), str(row["index"])) for _, row in pass_rows.iterrows())

    # 输出检测报告
    dict_extra = dict_target_pairs - pass_target_pairs
    pass_extra = pass_target_pairs - dict_target_pairs
    print(f"\n=== 一一对应检测报告 ===")
    print(f"\n=== One-to-one correspondence detection report ===")
    print(f"匹配行数量：{len(dict_target_pairs & pass_target_pairs)}")
    print(f"Number of matched rows: {len(dict_target_pairs & pass_target_pairs)}")
    print(f"字典有但whether_pass≠True的行：{dict_extra if dict_extra else '无'}")
    print(f"Rows in dictionary but whether_pass≠True: {dict_extra if dict_extra else 'None'}")
    print(f"whether_pass=True但字典没有的行：{pass_extra if pass_extra else '无'}")
    print(f"Rows with whether_pass=True but not in dictionary: {pass_extra if pass_extra else 'None'}")
    print(f"是否完全一一对应：{'是' if len(dict_extra) == 0 and len(pass_extra) == 0 else '否'}")
    print(f"Completely one-to-one correspondence: {'Yes' if len(dict_extra) == 0 and len(pass_extra) == 0 else 'No'}")

    def is_keep_row(row):
        """判断行是否需要保留"""
        if row.name == 0:  # 原始蛋白行（index=0）
            return True
        # 非原始蛋白行：判断是否在字典目标对中
        row_pair = (str(row["backbone"]), str(row["index"]))
        return row_pair in dict_target_pairs

    # 筛选保留行
    df_filtered = df[df.apply(is_keep_row, axis=1)].copy()
    print(f"\n=== 筛选后的数据统计 ===")
    print(f"\n=== Statistics after filtering ===")
    print(f"原CSV总行数：{len(df)}")
    print(f"Total rows in original CSV: {len(df)}")
    print(f"筛选后总行数：{len(df_filtered)}（原始蛋白1行 + 字典目标行{len(df_filtered) - 1}行）")
    print(f"Total rows after filtering: {len(df_filtered)} (Original protein 1 row + Dictionary target rows {len(df_filtered) - 1})")

    # 2. 添加af2_ptm和af2_plddt列
    def get_ptm_plddt_val(row):
        """根据行提取af2_ptm/af2_plddt值"""
        if row.name == 0:  # 原始蛋白行：设为'-'
            return '-', '-'
        # 字典目标行：从字典中提取值
        backbone = str(row["backbone"])
        idx = str(row["index"])
        if backbone in score_dict and idx in score_dict[backbone]:
            return score_dict[backbone][idx]["ptm"], score_dict[backbone][idx]["plddt"]
        return np.nan, np.nan  # 理论上不会走到这里（因已筛选）

    # 应用函数生成新列
    new_cols = df_filtered.apply(get_ptm_plddt_val, axis=1)
    df_filtered["af2_ptm"] = [x[0] for x in new_cols]
    df_filtered["af2_plddt"] = [x[1] for x in new_cols]

    # 3. 删除临时列（whether_pass_bool），保持CSV整洁
    cols_to_drop = ["ss_filter", "whether_pass", "whether_pass_bool"]
    df_filtered = df_filtered.drop(columns=cols_to_drop, errors="ignore")  # errors=ignore避免列不存在报错

    # 1. 查看筛选后的数据（重点展示关键列）
    print(f"\n=== 筛选后的数据预览（前5行）===")
    print(f"\n=== Preview of filtered data (first 5 rows) ===")
    display_cols = ["backbone", "index", "af2_ptm", "af2_plddt"]
    print(df_filtered[display_cols].head())

    # 2. 保存新CSV文件
    work_dir = esmfold_report_path.rsplit("/", 1)[0]
    output_csv_path = os.path.join(work_dir, "alphafold2_report.csv")
    df_filtered.to_csv(output_csv_path, index=False, encoding="utf-8")
    #print(f"\n=== 新CSV文件已保存 ===")
    #print(f"保存路径：{output_csv_path}")
    #print(f"文件包含：原始蛋白行 + 字典目标行（共{len(df_filtered)}行）")
    #print(f"新增列：af2_ptm（原始蛋白为'-'）、af2_plddt（原始蛋白为'-'）")
    return output_csv_path


def add_ss_data(alphafold2_report_path, dssp_csv_path, start_res, end_res, ss_threshold, ss=None):
    def csv_column_ratio_with_list(csv_path, start_idx, end_idx, column):
        """
        提取CSV指定列、指定行范围的元素列表，并统计各元素占比
        参数：
            csv_path: CSV文件路径（绝对/相对）
            start_idx: 开始索引（从1开始）
            end_idx: 结束索引（从1开始，包含）
            column: 列名（如"性别"）或列索引（从0开始，如0）
        返回：
            字典：
                - 元素列表：指定范围的原始元素（含空值，Python原生列表）
                - 占比统计：{元素: 占比（保留2位小数）}，强制包含H/E/C（无则0.00）
        示例：
            >>> result = csv_column_ratio_with_list("test.csv", 1, 10, "类型")
            >>> print(result["占比统计"])  # {'H': 0.30, 'E': 0.00, 'C': 0.50, 'A': 0.20}
        """
        # 异常处理：检查文件是否存在
        if not os.path.exists(csv_path):
            raise FileNotFoundError(f"CSV文件不存在：{csv_path}")

        # 1. 读取CSV文件（跳过空行，保留原始索引）
        try:
            df = pd.read_csv(csv_path, skip_blank_lines=True)
        except Exception as e:
            raise ValueError(f"读取CSV失败：{str(e)}. csv_path: {csv_path}")

        # 异常处理：空DataFrame
        if df.empty:
            raise ValueError("CSV文件无有效数据（空文件/全是空行）")

        # 2. 验证索引范围（从1开始→转换为0开始的切片）
        total_rows = len(df)
        if start_idx < 1 or end_idx < 1:
            raise ValueError("开始/结束索引必须≥1")
        if start_idx > end_idx:
            raise ValueError("开始索引不能大于结束索引")
        if end_idx > total_rows:
            raise ValueError(f"结束索引超出文件总行数（总行数：{total_rows}）")

        # 转换为0开始的切片（左闭右开，end_idx无需-1）
        slice_start = start_idx - 1
        slice_end = end_idx

        # 3. 提取指定列的指定行范围数据
        try:
            # 按列名/列索引提取列
            if isinstance(column, int):
                if column >= len(df.columns):
                    raise IndexError(f"列索引{column}超出范围（文件列数：{len(df.columns)}）")
                target_column = df.iloc[:, column]
            else:
                if column not in df.columns:
                    raise KeyError(f"列名'{column}'不存在（文件列名：{df.columns.tolist()}）")
                target_column = df[column]

            # ① 原始元素列表（含空值，保留原生类型）
            raw_element_list = target_column.iloc[slice_start:slice_end].tolist()
            # ② 去空值的列表（用于统计占比）
            clean_element_list = [x for x in raw_element_list if pd.notna(x)]

        except (KeyError, IndexError) as e:
            raise e
        except Exception as e:
            raise ValueError(f"提取列数据失败：{str(e)}")

        # 4. 统计元素占比（强制包含H/E/C，保留2位小数）
        required_elements = ['H', 'E', 'C']  # 必须包含的元素
        ratio_result = {}
        total_clean = len(clean_element_list)

        if total_clean == 0:
            # 无有效数据时，H/E/C占比均为0
            ratio_result = {elem: 0.00 for elem in required_elements}
        else:
            # 步骤1：初始化占比字典，先将H/E/C设为0
            ratio_result = {elem: 0.00 for elem in required_elements}
            # 步骤2：计算所有元素的占比（包括H/E/C和其他元素）
            count_series = pd.Series(clean_element_list).value_counts()
            for elem, count in count_series.items():
                ratio_result[str(elem)] = round(count / total_clean, 4)

        # 5. 按文档要求返回字典
        str_element = ''.join(raw_element_list)
        return str_element, ratio_result

    # 首先读取esmfold_report.csv文件
    df = pd.read_csv(alphafold2_report_path)

    # 初始化结果列表
    design_ss3 = []
    design_ss8 = []
    h_prop = []
    e_prop = []
    c_prop = []

    # 构建dssp csv文件路径映射
    dssp_file_map = {}
    backbone_folders = os.listdir(dssp_csv_path)
    backbone_folders = sorted(backbone_folders, key=natural_sort_key)

    for backbone_folder in backbone_folders:
        backbone_folder_path = os.path.join(dssp_csv_path, backbone_folder)

        dssp_csv_files = os.listdir(backbone_folder_path)
        dssp_csv_files = sorted(dssp_csv_files, key=natural_sort_key)

        for dssp_csv_file in dssp_csv_files:
            # 从文件名中提取index（去掉.csv后缀）
            index = os.path.splitext(dssp_csv_file)[0]
            dssp_csv_file_path = os.path.join(backbone_folder_path, dssp_csv_file)
            dssp_file_map[index] = dssp_csv_file_path

    # 遍历df的每一行
    for _, row in df.iterrows():
        index = row['index']

        # 检查是否为原始蛋白（af2_ptm或af2_plddt为'-'）
        is_original = False
        if isinstance(row.get('af2_ptm'), str) and row['af2_ptm'] == '-':
            is_original = True
        if isinstance(row.get('af2_plddt'), str) and row['af2_plddt'] == '-':
            is_original = True

        if is_original:
            # 对于原始蛋白，使用'-'
            design_ss8.append('-')
            design_ss3.append('-')
            h_prop.append('-')
            e_prop.append('-')
            c_prop.append('-')
        else:
            # 对于设计序列，查找对应的dssp csv文件
            if index in dssp_file_map:
                dssp_csv_file_path = dssp_file_map[index]
                try:
                    design_seq8, useless = csv_column_ratio_with_list(
                        csv_path=dssp_csv_file_path,
                        start_idx=start_res,
                        end_idx=end_res,
                        column='SS_8'
                    )
                    design_seq, dict_prop = csv_column_ratio_with_list(
                        csv_path=dssp_csv_file_path,
                        start_idx=start_res,
                        end_idx=end_res,
                        column='SS_3'
                    )
                    design_ss8.append(design_seq8)
                    design_ss3.append(design_seq)
                    h_prop.append(dict_prop.get('H', 0.0))
                    e_prop.append(dict_prop.get('E', 0.0))
                    c_prop.append(dict_prop.get('C', 0.0))
                except Exception as e:
                    # 如果处理失败，使用默认值
                    print(f"处理{index}的dssp_csv文件时出错: {e}")
                    print(f"Error processing dssp_csv file for {index}: {e}")
                    design_ss8.append('')
                    design_ss3.append('')
                    h_prop.append(0.0)
                    e_prop.append(0.0)
                    c_prop.append(0.0)
            else:
                # 如果找不到对应的dssp csv文件，使用默认值
                design_ss8.append('')
                design_ss3.append('')
                h_prop.append(0.0)
                e_prop.append(0.0)
                c_prop.append(0.0)

    # 将结果赋值给df
    df['af2_ss8'] = design_ss8
    df['af2_ss3'] = design_ss3
    df['af2_H_prop'] = h_prop
    df['af2_E_prop'] = e_prop
    df['af2_C_prop'] = c_prop
    ss_filter = []
    if ss is not None and ss != 'None':
        if ss.lower() == 'helix':
            for prop in h_prop:
                if prop == '-':
                    ss_filter.append('-')
                else:
                    ss_filter.append(prop > ss_threshold)
            df['ss_filter'] = ss_filter
        elif ss.lower() == 'strand':
            for prop in e_prop:
                if prop == '-':
                    ss_filter.append('-')
                else:
                    ss_filter.append(prop > ss_threshold)
            df['ss_filter'] = ss_filter
        else:
            raise ValueError("输入的二级结构类型错误，正确的输入应为helix或strand，请仔细检查\n"
                             "The entered secondary structure type is incorrect. "
                             "The correct input should be 'helix' or 'strand'. Please check carefully.")

    df.to_csv(alphafold2_report_path, index=False)

    return

#处理alphafold2的CSV文件，根据阈值进行筛选，按backbone分类的字典返回
def process_alphafold2_csv(csv_input_path, fasta_output_path,
                        plddt_threshold=None, ptm_threshold=None):
    """
    处理alphafold2的CSV文件，按backbone分类的字典返回

    新增返回值:
    dict: 按backbone分组的字典，结构为 {backbone: {index: sequence, ...}, ...}
    """

    def clean_to_bool(val):
        # 处理布尔型：直接保留
        if isinstance(val, bool):
            return val
        # 处理字符串型：'True'→True，'False'→False，其他（如'-'）→True
        elif isinstance(val, str):
            val = val.strip()  # 去除首尾空格
            if val == 'False':
                return False
            elif val == 'True':
                return True
            else:  # '-'、空字符串等
                return True
        # 其他类型（如NaN）→True
        else:
            return True

    # 读取CSV并验证必要列（新增backbone列验证）
    csv_path = Path(csv_input_path)
    # 前置校验：文件是否存在
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV文件不存在：{csv_path}")

    # 前置校验：阈值合法性（0~1范围，非None时校验）
    def check_threshold(threshold, name):
        if threshold is not None:
            if not isinstance(threshold, (int, float)):
                raise ValueError(f"{name}阈值必须为数值型，当前类型：{type(threshold)}")

    check_threshold(ptm_threshold, "PTM")
    check_threshold(plddt_threshold, "pLDDT")

    # 【优化1】自动创建原文件备份（后缀.bak），避免覆盖丢失数据
    # bak_path = csv_path.with_suffix(".csv.bak")
    # shutil.copy2(csv_path, bak_path)
    # print(f"已创建原CSV文件备份：{bak_path}")

    # 1. 读取CSV文件并校验必要列
    print(f"\n正在读取CSV文件: {csv_path}")
    print(f"\nReading CSV file: {csv_path}")
    df = pd.read_csv(csv_path).copy()
    required_columns = ['af2_ptm', 'af2_plddt', 'index', 'sequence', 'backbone']
    missing_columns = [col for col in required_columns if col not in df.columns]

    if missing_columns:
        raise ValueError(f"CSV文件缺少必要的列: {', '.join(missing_columns)}")
    total_rows = len(df)
    print(f"CSV文件总数据条数：{total_rows}")
    print(f"Total number of data entries in CSV file: {total_rows}")

    # 2. 数据类型转换与筛选逻辑
    # 创建原始值的副本，用于判断是否为原始蛋白
    df['ptm_original'] = df['af2_ptm']
    df['plddt_original'] = df['af2_plddt']

    # 转换为数值类型，非法值转NaN
    df['af2_ptm'] = pd.to_numeric(df['af2_ptm'], errors='coerce')
    df['af2_plddt'] = pd.to_numeric(df['af2_plddt'], errors='coerce')

    # 【优化3】缺失值统计提示，知晓数据质量
    ptm_nan_num = df['af2_ptm'].isna().sum()
    plddt_nan_num = df['af2_plddt'].isna().sum()
    print(f"af2_ptm转换后缺失值(NaN)条数：{ptm_nan_num}")
    print(f"Number of missing values (NaN) after af2_ptm conversion: {ptm_nan_num}")
    print(f"af2_plddt转换后缺失值(NaN)条数：{plddt_nan_num}")
    print(f"Number of missing values (NaN) after af2_plddt conversion: {plddt_nan_num}")

    # 构建筛选条件（初始为False）
    pass_condition = pd.Series(False, index=df.index)

    # 多规则标记原始蛋白（任一条件满足即视为原始蛋白）
    is_original_protein = pd.Series(False, index=df.index)
    is_original_protein.iloc[0] = True  # 第一行视为原始蛋白
    is_original_protein |= df['ptm_original'] == '-'  # ptm为'-'
    is_original_protein |= df['plddt_original'] == '-'  # plddt为'-'
    is_original_protein |= df['af2_ptm'].isna()  # ptm为NaN
    is_original_protein |= df['af2_plddt'].isna()  # plddt为NaN
    original_protein_num = is_original_protein.sum()
    print(f"识别到原始蛋白条数：{original_protein_num}")
    print(f"Number of original protein entries identified: {original_protein_num}")

    # 应用不同筛选条件分支
    # print("df['ss_filter']: ", df['ss_filter'])
    # print("df['ss_filter'].apply(clean_to_bool): ",df['ss_filter'].apply(clean_to_bool))
    if plddt_threshold is not None and ptm_threshold is None:
        print(f"  - 仅plddt筛选，阈值: > {plddt_threshold}")
        print(f"  - Only plddt filtering, threshold: > {plddt_threshold}")
        if 'ss_filter' in df.columns:
            pass_condition = ((df['af2_plddt'] > plddt_threshold) &
                              df['ss_filter'].apply(clean_to_bool))
        else:
            pass_condition = df['af2_plddt'] > plddt_threshold

    elif ptm_threshold is not None and plddt_threshold is None:
        print(f"  - 仅ptm筛选，阈值: > {ptm_threshold}")
        print(f"  - Only ptm filtering, threshold: > {ptm_threshold}")
        if 'ss_filter' in df.columns:
            pass_condition = ((df['af2_ptm'] > ptm_threshold) &
                              df['ss_filter'].apply(clean_to_bool))
        else:
            pass_condition = df['af2_ptm'] > ptm_threshold

    elif plddt_threshold is not None and ptm_threshold is not None:
        print(f"  - 双条件筛选：ptm > {ptm_threshold} AND plddt > {plddt_threshold}")
        print(f"  - Dual condition filtering: ptm > {ptm_threshold} AND plddt > {plddt_threshold}")
        if 'ss_filter' in df.columns:
            pass_condition = ((df['af2_ptm'] > ptm_threshold) & (df['af2_plddt'] > plddt_threshold) &
                              df['ss_filter'].apply(clean_to_bool))
        else:
            pass_condition = (df['af2_ptm'] > ptm_threshold) & (df['af2_plddt'] > plddt_threshold)
    else:
        print(f"  - 不进行严格筛选，仅过滤plddt≤0的无效值")
        print(f"  - No strict filtering, only filter invalid values with plddt≤0")
        pass_condition = ((df['af2_plddt'] > 0) &
                          df['ss_filter'].apply(clean_to_bool))

    # 原始蛋白始终通过筛选
    pass_condition |= is_original_protein

    # 3. 更新CSV文件（添加whether_pass列）
    df['whether_pass'] = pass_condition
    df.loc[is_original_protein, 'whether_pass'] = '-'  # 原始蛋白设为'-'
    df.loc[is_original_protein, 'af2_ptm'] = '-'
    df.loc[is_original_protein, 'af2_plddt'] = '-'
    df.loc[is_original_protein, 'ss_filter'] = '-'
    df = df.drop(['ptm_original', 'plddt_original'], axis=1)  # 移除临时列

    # 【优化4】指定utf-8编码保存，避免编码混乱，消除潜在报错
    df.to_csv(csv_path, index=False, encoding='utf-8')
    print(f"\n已更新CSV文件，新增whether_pass列: {csv_path}")
    print(f"\nCSV file updated, added whether_pass column: {csv_path}")

    # 【优化5】详细筛选结果统计，无需手动查看CSV
    total_pass = (df['whether_pass'] != 'No').sum()  # '-'和True都算通过
    design_pass = ((df['whether_pass'] == True) & (~is_original_protein)).sum()
    design_total = total_rows - original_protein_num
    print("=" * 50)
    print("筛选结果统计：")
    print("Filtering result statistics:")
    print(f"  总通过条数（含原始蛋白）：{total_pass} / {total_rows}")
    print(f"  Total passed entries (including original protein): {total_pass} / {total_rows}")
    print(f"  设计结果总条数：{design_total}")
    print(f"  Total design results: {design_total}")
    print(f"  设计结果通过条数：{design_pass} / {design_total}")
    print(f"  Passed design results: {design_pass} / {design_total}")
    print(f"  原始蛋白条数（默认通过）：{original_protein_num}")
    print(f"  Original protein entries (default passed): {original_protein_num}")
    print("=" * 50)

    # 4. 生成FASTA文件（仅通过筛选的序列，排除原始蛋白）
    # 筛选通过的序列：whether_pass为True或'-'（原始蛋白）
    passed_df = df[(df['whether_pass'] == True) | (df['whether_pass'] == '-')].copy()
    # 排除原始蛋白
    non_original_df = passed_df[~is_original_protein].copy()

    print(f"\n筛选结果统计:")
    print(f"\nFiltering result statistics:")
    print(f"  - 总序列数: {len(df)}")
    print(f"  - Total sequences: {len(df)}")
    print(f"  - 通过筛选序列数: {len(passed_df)}")
    print(f"  - Number of sequences passed filtering: {len(passed_df)}")
    print(f"  - 其中原始蛋白数: {sum(is_original_protein)}")
    print(f"  - Number of original proteins: {sum(is_original_protein)}")
    print(f"  - 其中设计序列数: {len(non_original_df)}")
    print(f"  - Number of designed sequences: {len(non_original_df)}")

    if len(non_original_df) > 0:
        with open(fasta_output_path, 'w', encoding='utf-8') as f:
            for _, row in non_original_df.iterrows():
                f.write(f">{row['index']} pTM={row['af2_ptm']} pLDDT={row['af2_plddt']}\n{row['sequence']}\n")
        print(f"已生成FASTA文件: {fasta_output_path}")
        print(f"FASTA file generated: {fasta_output_path}")
    else:
        print("\n警告：无设计序列通过筛选，未生成FASTA文件")
        print("\nWarning: No designed sequences passed filtering, FASTA file not generated")

    # 5. 核心新增：生成按backbone分类的字典
    backbone_dict = {}
    # 遍历通过筛选的行，按backbone分组
    for _, row in passed_df.iterrows():
        backbone = row['backbone']
        seq_index = row['index']  # 保留原index作为键
        sequence = row['sequence']

        # 对于原始蛋白，使用None作为分数
        if is_original_protein.iloc[_]:
            af2_ptm = None
            af2_plddt = None
        else:
            af2_ptm = row['af2_ptm']
            af2_plddt = row['af2_plddt']

        # 若backbone未在字典中，初始化空字典；否则添加index-sequence键值对
        if backbone not in backbone_dict:
            backbone_dict[backbone] = {}
        if seq_index not in backbone_dict[backbone]:
            backbone_dict[backbone][seq_index] = {}
        backbone_dict[backbone][seq_index]['sequence'] = sequence
        backbone_dict[backbone][seq_index]['ptm'] = af2_ptm
        backbone_dict[backbone][seq_index]['plddt'] = af2_plddt

    print(f"\n已生成按backbone分类的字典，包含 {len(backbone_dict)} 个不同backbone")
    print(f"\nGenerated dictionary classified by backbone, containing {len(backbone_dict)} different backbones")
    return backbone_dict, len(non_original_df)  # 返回最终字典

#生成筛选后的文件夹
def filter_files(alphafold2_folder, filter_folder, backbone_dict):
    if not os.path.exists(filter_folder):
        os.makedirs(filter_folder, exist_ok=True)

    for backbone, values in backbone_dict.items():
        filter_backbone_folder = os.path.join(filter_folder, backbone)
        if not os.path.exists(filter_backbone_folder):
            os.makedirs(filter_backbone_folder, exist_ok=True)
        with open(f'{filter_backbone_folder}/{backbone}_filter.fa', 'a+', encoding='utf-8') as f:
            f.truncate(0)
            for index, values2 in values.items():
                # 只处理非原始蛋白，因为原始蛋白可能没有对应的PDB文件
                if values2['ptm'] is not None and values2['plddt'] is not None:
                    file_path = os.path.join(alphafold2_folder, 'structure_prediction_files', backbone, f'{index}.pdb')
                    copy_path = os.path.join(filter_folder, backbone, f'{index}.pdb')
                    if os.path.exists(file_path):
                        shutil.copy(file_path, copy_path)
                        f.write(f">{index} pTM={values2['ptm']} pLDDT={values2['plddt']}\n")
                        f.write(f"{values2['sequence']}\n")

    return


def get_start_end(input_str):
    """
    提取输入中的开始数字和结束数字
    """
    if " " in input_str:
        num_list = [int(num) for num in input_str.split()]
        return num_list[0], num_list[-1]
    elif "-" in input_str:
        match = re.match(r"^[A-Za-z]*(\d+)-(\d+)$", input_str)
        if match:
            start = int(match.group(1))
            end = int(match.group(2))
            return start, end
    else:
        match = re.match(r"^[A-Za-z]*(\d+)$", input_str)
        if match:
            num = int(match.group(1))
            return num, num
    return None, None


def main():
    args = parse_args()
    esmfold_report_path = args.esmfold_report_path
    alphafold2_folder = args.alphafold2_folder
    plddt_threshold = args.plddt_threshold
    ptm_threshold = args.ptm_threshold
    seq_range_str = args.seq_range_str
    ss = args.ss
    ss_threshold = args.ss_threshold
    ss_filter = args.ss_filter

    start, end = get_start_end(seq_range_str)

    structure_prediction_files_folder = os.path.join(alphafold2_folder, 'structure_prediction_files')
    dssp_folder = os.path.join(alphafold2_folder, 'dssp_files')
    csv_folder = os.path.join(alphafold2_folder, 'csv_files')

    ndx_dict = make_dssp_csv(
        structure_prediction_files_folder=structure_prediction_files_folder,
        dssp_folder=dssp_folder,
        csv_folder=csv_folder,
    )
    print('ndx_dict:', ndx_dict)

    score_dict = read_ptm_and_plddt(structure_prediction_files_folder, ndx_dict)
    print('score_dict:', score_dict)

    alphafold2_report_path = add_ptm_and_plddt(score_dict, esmfold_report_path)

    if not ss_filter:
        ss = None
    if ss is None or ss == 'None' or ss == '':
        print('是否开启二级结构过滤？\nNo')
        print('Is secondary structure filtering enabled?\nNo')
    else:
        print('是否开启二级结构过滤？\nYes')
        print('Is secondary structure filtering enabled?\nYes')
        print(f'目标二级结构类型：{ss}')
        print(f'Target secondary structure type: {ss}')
        print(f'过滤阈值：{ss_threshold}')
        print(f'Filtering threshold: {ss_threshold}')
    
    add_ss_data(
        alphafold2_report_path=alphafold2_report_path, 
        dssp_csv_path=csv_folder,
        start_res=start, 
        end_res=end, 
        ss_threshold=ss_threshold, 
        ss=ss
    )

    backbone_dict, seqs_num = process_alphafold2_csv(
        csv_input_path=alphafold2_report_path,
        fasta_output_path=f'{alphafold2_folder}/filter_result.fa',
        plddt_threshold=plddt_threshold,
        ptm_threshold=ptm_threshold
    )
    if seqs_num > 0:
        filter_files(
            alphafold2_folder=alphafold2_folder,
            filter_folder=f'{alphafold2_folder}/filter_files',
            backbone_dict=backbone_dict
        )
    else:
        print(f"没有蛋白质序列通过筛选，不进行后续操作。")
        print(f"No protein sequences passed filtering, no further operations will be performed.")
        sys.exit(102)

    return 


if __name__ == '__main__':
    main()



