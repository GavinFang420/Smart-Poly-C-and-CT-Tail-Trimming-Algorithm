# Import Libraries
import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import warnings
warnings.filterwarnings('ignore')

# The FatqReader class reads FASTQ files and filters out contaminated reads.
class FastqReader:
    def __init__(self, filename: str):
        self.filename = filename
        
    def read_records(self, max_records: int = 9999, filter_contaminated: bool = True):
        records = []
        contaminated = 0
        try:
            with open(self.filename, 'r') as f:
                while len(records) < max_records:
                    header = f.readline().strip()
                    if not header:
                        break
                    sequence = f.readline().strip()
                    plus = f.readline().strip()
                    quality = f.readline().strip()
                    if not all([header, sequence, plus, quality]):
                        break
                    if filter_contaminated and self.is_contaminated(sequence):
                        contaminated += 1
                        continue
                    records.append({'header': header, 'sequence': sequence, 'quality': quality})
        except Exception as e:
            print(f"Error reading {self.filename}: {e}")
            return []
        if contaminated > 0:
            print(f"Filtered {contaminated} contaminated reads from {self.filename}")
        return records

    def is_contaminated(self, sequence: str, threshold: float = 0.75) -> bool:
        if not sequence:
            return True
        length = len(sequence)
        counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
        for base in sequence.upper():
            if base in counts:
                counts[base] += 1
        
        # Check if any single base (A, C, G, T) exceeds threshold
        if any(counts[base]/length > threshold for base in ['A', 'C', 'G', 'T']):
            return True
        
        # Check if N bases exceed 20%
        if counts['N']/length >= 0.2:
            return True
        
        return False
    

class ConstantTrimOptimizer:
   def __init__(self, original_length: int = 150):
       self.original_length = original_length
       self.min_retention_rate = 0
       self.max_retention_rate = 0
       self.min_quality_score = 0
       self.max_quality_score = 0

   def rescale_value(self, value: float, min_value: float, max_value: float) -> float:
       """
       通用的 rescale 函数，接受一个值以及最小最大值范围，
       将其归一化到 [0, 1] 范围。
       """
       if max_value == min_value:
           return 0.5
       return (value - min_value) / (max_value - min_value)

   def calculate_data_retention_rate(self, trim_length: int) -> float:
       """
       计算留存率，给定修剪长度。修剪长度大于原始长度时，留存率为 0。
       """
       if trim_length >= self.original_length:
           return 0.0
       return ((self.original_length - trim_length) / self.original_length)
   
   def get_base_weights(self, from_start: str) -> Dict[str, float]:
       """
       获取不同碱基的权重，用于质量评估
       R1: C和T是关键碱基（CT tail），权重高
       R2: G和A是关键碱基（GG head，A是T的互补），权重高
       """
       if from_start == "tail":  # R1 tail trimming
           return {
               'C': 0.65,  # C最重要（CT tail中的C）
               'T': 0.15,  # T次重要
               'A': 0.05,  # 非关键碱基
               'G': 0.05,  # 非关键碱基
               'N': 0.10   # N碱基权重
           }
       else:  # R2 head trimming  
           return {
               'G': 0.65,  # G最重要（GG head中的G）
               'A': 0.15,  # A次重要（T的互补）
               'T': 0.05,  # 非关键碱基
               'C': 0.05,  # 非关键碱基
               'N': 0.10   # N碱基权重
           }

   def get_theoretical_baseline(self, from_start: str) -> Dict[str, float]:
       """
       使用理论上的亚硫酸盐测序期望分布
       """
       if from_start == "head":  # R2 trimming后期望
           return {'A': 0.425, 'T': 0.325, 'G': 0.025, 'C': 0.225}
       else:  # R1 trimming后期望  
           return {'A': 0.325, 'T': 0.425, 'G': 0.225, 'C': 0.025}

   def calculate_weighted_quality_score(self, sequence: str, from_start: Optional[str] = None) -> float:
       """
       使用加权方法的质量评估，重点关注CT tail相关碱基
       """
       if not sequence:
           return 0.0
       
       length = len(sequence)
       counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'N': 0}
       for base in sequence.upper():
           if base in counts:
               counts[base] += 1
       
       ratios = {base: count/length for base, count in counts.items()}
       
       # 获取理论基线和权重
       expected = self.get_theoretical_baseline(from_start)
       weights = self.get_base_weights(from_start)
       
       # 根据trimming程度调整期望
       if from_start == "head":  # R2 head trimming，期望G含量降低
           expected = expected.copy()
           expected['G'] *= 0.3  # G含量应该显著降低
           # 重新分配
           g_reduction = (0.025 - expected['G'])
           expected['A'] += g_reduction * 0.6
           expected['T'] += g_reduction * 0.4
           
       elif from_start == "tail":  # R1 tail trimming，期望C含量降低
           expected = expected.copy()
           expected['C'] *= 0.3  # C含量应该显著降低
           # 重新分配
           c_reduction = (0.025 - expected['C'])
           expected['T'] += c_reduction * 0.6
           expected['A'] += c_reduction * 0.4
       
       # 计算加权质量分数
       weighted_score = 0.0
       for base in ['A', 'T', 'G', 'C', 'N']:
           if base in expected:
               deviation = abs(ratios[base] - expected[base])
               # 偏差越小越好，转换为分数
               base_score = max(0, 1.0 - deviation * 4)  # 25%偏差 = 0分
               weighted_score += base_score * weights[base]
           else:  # N base
               # N含量越少越好
               n_score = max(0, 1.0 - ratios['N'] * 5)  # 20% N = 0分
               weighted_score += n_score * weights['N']
       
       return weighted_score

   def calculate_quality_improvement_rate(self, sequences: List[str], trim_length: int, from_start: Optional[str] = None) -> float:
       """
       针对亚硫酸盐测序数据计算质量改善率
       """
       if not sequences or trim_length <= 0:
           return 0.0
       
       # 修剪序列
       trimmed_seqs = []
       for seq in sequences:
           if from_start == "head":
               trimmed_seqs.append(seq[trim_length:] if trim_length < len(seq) else "")
           elif from_start == "tail":
               trimmed_seqs.append(seq[:-trim_length] if trim_length < len(seq) else "")
           else:
               trimmed_seqs.append(seq)
       
       valid_seqs = [seq for seq in trimmed_seqs if seq and len(seq) >= 50]
       if not valid_seqs:
           return 0.0
       
       total_score = 0.0
       for seq in valid_seqs:
           # 分析区域选择
           if from_start == "tail":
               analyze_region = seq[-80:] if len(seq) >= 80 else seq
           else:
               analyze_region = seq[:80] if len(seq) >= 80 else seq
           
           if not analyze_region:
               continue
               
           # 使用加权质量评估
           score = self.calculate_weighted_quality_score(analyze_region, from_start)
           total_score += score
       
       return total_score / len(valid_seqs) if valid_seqs else 0.0

   def find_optimal_trim_length(self, sequences: List[str], from_start: str = "head", trim_range: Tuple[int, int] = (5, 25)) -> Dict:
       """
       找到最佳修剪长度，优化留存率和质量得分的乘积。
       """
       min_trim, max_trim = trim_range
       results = []

       # 计算留存率的最小值和最大值
       all_retention_rates = [self.calculate_data_retention_rate(trim_len) for trim_len in range(min_trim, max_trim + 1)]
       self.max_retention_rate = max(all_retention_rates)
       self.min_retention_rate = min(all_retention_rates)

       # 计算质量得分的最小值和最大值
       all_quality_scores = [self.calculate_quality_improvement_rate(sequences, trim_len, from_start) for trim_len in range(min_trim, max_trim + 1)]
       self.max_quality_score = max(all_quality_scores)
       self.min_quality_score = min(all_quality_scores)

       # 打印不同算法的最佳修剪长度
       self.print_best_trim_lengths(sequences, from_start, trim_range)

       for trim_len in range(min_trim, max_trim + 1):
           retention_rate = self.calculate_data_retention_rate(trim_len)
           quality_rate = self.calculate_quality_improvement_rate(sequences, trim_len, from_start)

           # 使用通用的 rescale 函数分别对留存率和质量得分进行归一化
           rescaled_retention_rate = self.rescale_value(retention_rate, self.min_retention_rate, self.max_retention_rate)
           rescaled_quality_rate = self.rescale_value(quality_rate, self.min_quality_score, self.max_quality_score)

           product_score = rescaled_retention_rate * rescaled_quality_rate  # 乘积计算
           results.append({
               'trim_length': trim_len,
               'retention_rate': retention_rate,
               'quality_rate': quality_rate,
               'rescaled_retention_rate': rescaled_retention_rate,
               'rescaled_quality_rate': rescaled_quality_rate,
               'product_score': product_score
           })
       return {'results': results}

   def print_best_trim_lengths(self, sequences: List[str], from_start: str, trim_range: Tuple[int, int]):
       """
       打印三种不同算法的最佳修剪长度
       """
       min_trim, max_trim = trim_range
       
       # 预计算所有需要的数值
       all_retention_rates = [self.calculate_data_retention_rate(trim_len) for trim_len in range(min_trim, max_trim + 1)]
       all_quality_scores = [self.calculate_quality_improvement_rate(sequences, trim_len, from_start) for trim_len in range(min_trim, max_trim + 1)]
       
       max_retention = max(all_retention_rates)
       min_retention = min(all_retention_rates)
       max_quality = max(all_quality_scores)
       min_quality = min(all_quality_scores)
       
       # 方法1: quality score * retention rate (both scaled) - 推荐
       scaled_products = []
       for i, trim_len in enumerate(range(min_trim, max_trim + 1)):
           scaled_retention = self.rescale_value(all_retention_rates[i], min_retention, max_retention)
           scaled_quality = self.rescale_value(all_quality_scores[i], min_quality, max_quality)
           scaled_products.append(scaled_retention * scaled_quality)
       
       best_scaled_idx = np.argmax(scaled_products)
       best_scaled_trim = min_trim + best_scaled_idx
       
       # 方法2: quality score squared * retention rate (both scaled) - 偏好质量
       squared_quality_products = []
       for i, trim_len in enumerate(range(min_trim, max_trim + 1)):
           scaled_retention = self.rescale_value(all_retention_rates[i], min_retention, max_retention)
           scaled_quality = self.rescale_value(all_quality_scores[i], min_quality, max_quality)
           squared_quality_products.append(scaled_retention * (scaled_quality ** 2))
       
       best_quality_squared_idx = np.argmax(squared_quality_products)
       best_quality_squared_trim = min_trim + best_quality_squared_idx
       
       # 方法3: retention rate squared * quality score (both scaled) - 偏好保留率
       squared_retention_products = []
       for i, trim_len in enumerate(range(min_trim, max_trim + 1)):
           scaled_retention = self.rescale_value(all_retention_rates[i], min_retention, max_retention)
           scaled_quality = self.rescale_value(all_quality_scores[i], min_quality, max_quality)
           squared_retention_products.append((scaled_retention ** 2) * scaled_quality)
       
       best_retention_squared_idx = np.argmax(squared_retention_products)
       best_retention_squared_trim = min_trim + best_retention_squared_idx
       
       print(f"{from_start.upper()} Optimal Trim Lengths:")
       print(f"Quality*Retention (both scaled) - Recommended: {best_scaled_trim}bp")
       print(f"Quality^2*Retention (both scaled) - Prefer Quality: {best_quality_squared_trim}bp")
       print(f"Retention^2*Quality (both scaled) - Prefer Retention: {best_retention_squared_trim}bp")
       
       # 保存最佳结果以供后续使用
       self.best_scaled_trim = best_scaled_trim
       self.best_scaled_retention = all_retention_rates[best_scaled_idx]
       self.best_scaled_quality = all_quality_scores[best_scaled_idx]
       
       # 保存三个最佳结果用于绘图
       self.all_best_trims = {
           'Quality*Retention (both scaled) - Recommended': best_scaled_trim,
           'Quality^2*Retention (both scaled) - Prefer Quality': best_quality_squared_trim,
           'Retention^2*Quality (both scaled) - Prefer Retention': best_retention_squared_trim
       }


class ConstantTrimVisualizer:
    def __init__(self, output_dir: str = "constanttrim_results"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        plt.style.use('seaborn-v0_8' if 'seaborn-v0_8' in plt.style.available else 'default')

    def rescale_value(self, value: float, min_value: float, max_value: float) -> float:
        """
        通用的 rescale 函数，接受一个值以及最小最大值范围，
        将其归一化到 [0, 1] 范围。
        """
        if max_value == min_value:
            return 0.5
        return (value - min_value) / (max_value - min_value)

    def calculate_three_products(self, retention_rates, quality_scores, min_retention, max_retention, min_quality, max_quality):
        """计算三种方法的乘积分数"""
        products = {
            'Quality*Retention (both scaled) - Recommended': [],
            'Quality^2*Retention (both scaled) - Prefer Quality': [],
            'Retention^2*Quality (both scaled) - Prefer Retention': []
        }
        
        for i in range(len(retention_rates)):
            # Scaled values
            scaled_retention = self.rescale_value(retention_rates[i], min_retention, max_retention)
            scaled_quality = self.rescale_value(quality_scores[i], min_quality, max_quality)
            
            # Three methods
            products['Quality*Retention (both scaled) - Recommended'].append(scaled_retention * scaled_quality)
            products['Quality^2*Retention (both scaled) - Prefer Quality'].append(scaled_retention * (scaled_quality ** 2))
            products['Retention^2*Quality (both scaled) - Prefer Retention'].append((scaled_retention ** 2) * scaled_quality)
            
        return products

    def plot_optimization_curves(self, analysis: Dict, dataset: str):
        fig, axes = plt.subplots(2, 2, figsize=(20, 12))
        fig.suptitle(f'{dataset.upper()} - ConstantTrim Optimization (Quality vs Retention)', fontsize=16)

        r1_results = analysis['r1_analysis']['results']
        r2_results = analysis['r2_analysis']['results']

        # 获取最大和最小质量得分
        max_quality_score = max([r['quality_rate'] for r in r1_results + r2_results])
        min_quality_score = min([r['quality_rate'] for r in r1_results + r2_results])
        
        # 获取最大和最小留存率
        max_retention_rate = max([r['retention_rate'] for r in r1_results + r2_results])
        min_retention_rate = min([r['retention_rate'] for r in r1_results + r2_results])

        # R1 Optimization (Plot 1) - 移除Top3标记
        ax = axes[0, 0]
        r1_trims = [r['trim_length'] for r in r1_results]
        r1_ret = [r['retention_rate'] for r in r1_results]
        r1_qual = [r['quality_rate'] for r in r1_results]
        ax2 = ax.twinx()
        ax.plot(r1_trims, r1_ret, 'b-o', linewidth=2, label='Data Retention Rate')
        ax2.plot(r1_trims, r1_qual, 'r-s', linewidth=2, label='Quality Score')
        ax.set_xlabel('R1 Tail Trim Length (bp)')
        ax.set_ylabel('Data Retention Rate', color='b')
        ax2.set_ylabel('Quality Score', color='r')
        ax.set_title('R1 Tail Trimming Optimization')
        lines, labels = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax.legend(lines+lines2, labels+labels2, loc='center right')
        ax.grid(True, alpha=0.3)

        # R2 Optimization (Plot 2) - 移除Top3标记
        ax = axes[0, 1]
        r2_trims = [r['trim_length'] for r in r2_results]
        r2_ret = [r['retention_rate'] for r in r2_results]
        r2_qual = [r['quality_rate'] for r in r2_results]
        ax2 = ax.twinx()
        ax.plot(r2_trims, r2_ret, 'b-o', linewidth=2, label='Data Retention Rate')
        ax2.plot(r2_trims, r2_qual, 'r-s', linewidth=2, label='Quality Score')
        ax.set_xlabel('R2 Head Trim Length (bp)')
        ax.set_ylabel('Data Retention Rate', color='b')
        ax2.set_ylabel('Quality Score', color='r')
        ax.set_title('R2 Head Trimming Optimization')
        lines, labels = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax.legend(lines+lines2, labels+labels2, loc='center right')
        ax.grid(True, alpha=0.3)

        # Combined Score (Plot 3 & 4) - 使用原始分数，标记三个最佳值
        # R1 Combined Scores
        ax = axes[1, 0]
        r1_products = self.calculate_three_products(r1_ret, r1_qual, min_retention_rate, max_retention_rate, min_quality_score, max_quality_score)
        
        colors = ['red', 'blue', 'green']
        linestyles = ['-', '--', '-.']
        
        for i, (method, values) in enumerate(r1_products.items()):
            ax.plot(r1_trims, values, color=colors[i], linestyle=linestyles[i], 
                   linewidth=2, marker='o', markersize=4, label=method)
            # 找到最佳值并标记
            best_idx = np.argmax(values)
            best_trim = r1_trims[best_idx]
            ax.axvline(best_trim, color=colors[i], linestyle=':', alpha=0.7)
            ax.scatter([best_trim], [values[best_idx]], s=70, color=colors[i], zorder=5)
            
        ax.set_xlabel('R1 Trim Length (bp)')
        ax.set_ylabel('Combined Score (Original Scale)')
        ax.set_title('R1 Combined Score Comparison')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.grid(True, alpha=0.3)

        # R2 Combined Scores
        ax = axes[1, 1]
        r2_products = self.calculate_three_products(r2_ret, r2_qual, min_retention_rate, max_retention_rate, min_quality_score, max_quality_score)
        
        for i, (method, values) in enumerate(r2_products.items()):
            ax.plot(r2_trims, values, color=colors[i], linestyle=linestyles[i], 
                   linewidth=2, marker='o', markersize=4, label=method)
            # 找到最佳值并标记
            best_idx = np.argmax(values)
            best_trim = r2_trims[best_idx]
            ax.axvline(best_trim, color=colors[i], linestyle=':', alpha=0.7)
            ax.scatter([best_trim], [values[best_idx]], s=70, color=colors[i], zorder=5)

        ax.set_xlabel('R2 Trim Length (bp)')
        ax.set_ylabel('Combined Score (Original Scale)')
        ax.set_title('R2 Combined Score Comparison')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(self.output_dir / f"{dataset}_constanttrim_optimization.png", dpi=300, bbox_inches='tight')
        plt.close()


# 添加输出重定向类
class PrintLogger:
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.log = open(filename, "w", encoding='utf-8')

    def write(self, message):
        self.terminal.write(message)
        # 移除空行
        if message.strip():
            self.log.write(message.rstrip() + '\n')

    def flush(self):
        self.terminal.flush()
        self.log.flush()
        
    def close(self):
        self.log.close()


def main():
    parser = argparse.ArgumentParser(description='ConstantTrim Optimization')
    parser.add_argument('--data-dir', type=str, default='.', help='Directory containing FASTQ files')
    parser.add_argument('--datasets', nargs='+', default=['tumor', 'normal'], help='Dataset names to analyze')
    parser.add_argument('--max-reads', type=int, default=2000, help='Maximum number of reads to analyze')
    parser.add_argument('--output-dir', type=str, default='constanttrim_results', help='Output directory for results')
    parser.add_argument('--trim-range', nargs=2, type=int, default=[5, 25], help='Trim length range to test')
    args = parser.parse_args()

    # 设置输出重定向
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    logger = PrintLogger(output_dir / "constanttrim_results.txt")
    sys.stdout = logger

    try:
        optimizer = ConstantTrimOptimizer(original_length=150)
        visualizer = ConstantTrimVisualizer(args.output_dir)

        for dataset in args.datasets:
            r1_file = f"{dataset}_R1.fastq"
            r2_file = f"{dataset}_R2.fastq"
            if not (os.path.exists(r1_file) and os.path.exists(r2_file)):
                continue
            r1_records = FastqReader(r1_file).read_records(args.max_reads)
            r2_records = FastqReader(r2_file).read_records(args.max_reads)
            if not r1_records or not r2_records:
                continue
            r1_seqs = [r['sequence'] for r in r1_records]
            r2_seqs = [r['sequence'] for r in r2_records]

            print(f"=== Analysis for {dataset.upper()} ===")
            
            # 分别分析 R1 和 R2 序列
            r1_analysis = optimizer.find_optimal_trim_length(r1_seqs, from_start="tail", trim_range=tuple(args.trim_range))
            r2_analysis = optimizer.find_optimal_trim_length(r2_seqs, from_start="head", trim_range=tuple(args.trim_range))

            # 将 R1 和 R2 的分析结果整合到一起
            analysis = {'r1_analysis': r1_analysis, 'r2_analysis': r2_analysis}
            
            # 打印Overall Best结果（基于both scaled相乘）
            print(f"=== Overall Best (Quality*Retention both scaled) for {dataset.upper()} ===")
            print(f"R1 Tail: {optimizer.best_scaled_trim}bp (Retention: {optimizer.best_scaled_retention:.4f}, Quality: {optimizer.best_scaled_quality:.6f})")
            
            # 为R2重新计算一遍以获取最佳值
            r2_optimizer = ConstantTrimOptimizer(original_length=150)
            _ = r2_optimizer.find_optimal_trim_length(r2_seqs, from_start="head", trim_range=tuple(args.trim_range))
            print(f"R2 Head: {r2_optimizer.best_scaled_trim}bp (Retention: {r2_optimizer.best_scaled_retention:.4f}, Quality: {r2_optimizer.best_scaled_quality:.6f})")
            
            # 可视化结果
            print(f"Generating plots for dataset: {dataset}")
            visualizer.plot_optimization_curves(analysis, dataset)

    finally:
        # 恢复标准输出并关闭日志文件
        sys.stdout = logger.terminal
        logger.close()

if __name__ == "__main__":
    main()