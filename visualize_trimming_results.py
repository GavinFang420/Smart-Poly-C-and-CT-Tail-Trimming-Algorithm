#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

def plot_basic_stats(sample_name):
    """绘制基本统计图"""
    csv_file = f"{sample_name}_trimming_analysis.csv"
    if not os.path.exists(csv_file):
        print(f"File {csv_file} not found")
        return
    
    data = pd.read_csv(csv_file)
    stats = dict(zip(data['metric'], data['value']))
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f'{sample_name} - Poly-C Trimming Analysis', fontsize=14)
    
    # 1. Trimming rate pie chart
    trimmed = stats.get('trimmed_pairs', 0)
    total = stats.get('total_pairs', 1)
    not_trimmed = total - trimmed
    
    axes[0,0].pie([trimmed, not_trimmed], 
                  labels=[f'Trimmed ({trimmed:,})', f'Not Trimmed ({not_trimmed:,})'],
                  autopct='%1.1f%%', colors=['red', 'blue'])
    axes[0,0].set_title('Trimming Rate')
    
    # 2. Trimming patterns
    r1_only = stats.get('r1_only_trimmed', 0)
    r2_only = stats.get('r2_only_trimmed', 0) 
    both = stats.get('both_trimmed', 0)
    
    axes[0,1].bar(['R1 Only', 'R2 Only', 'Both'], [r1_only, r2_only, both])
    axes[0,1].set_title('Trimming Patterns')
    axes[0,1].set_ylabel('Count')
    
    # 3. C-tail types
    pure_c = stats.get('pure_c_tails', 0)
    mutations = stats.get('c_tails_with_mutations', 0)
    
    axes[1,0].bar(['Pure C-tails', 'C-tails with mutations'], [pure_c, mutations])
    axes[1,0].set_title('C-tail Types')
    axes[1,0].set_ylabel('Count')
    
    # 4. Summary text
    axes[1,1].axis('off')
    summary = f"""Summary Statistics:
    Total pairs: {total:,}
    Trimmed: {trimmed:,} ({100*trimmed/total:.1f}%)
    Pure C-tails: {pure_c:,}
    C-tails with mutations: {mutations:,}
    """
    axes[1,1].text(0.1, 0.5, summary, fontsize=12, verticalalignment='center')
    
    plt.tight_layout()
    plt.savefig(f'{sample_name}_stats.png', dpi=150, bbox_inches='tight')
    plt.show()

def plot_lengths(sample_name):
    """绘制长度分布"""
    length_file = f"{sample_name}_trimming_analysis_lengths.csv"
    if not os.path.exists(length_file):
        print(f"File {length_file} not found")
        return
    
    data = pd.read_csv(length_file)
    
    plt.figure(figsize=(10, 6))
    
    r1_data = data[data['read_type'] == 'R1']['trim_length']
    r2_data = data[data['read_type'] == 'R2']['trim_length']
    
    if len(r1_data) > 0:
        plt.hist(r1_data, bins=30, alpha=0.7, label=f'R1 (n={len(r1_data)})', color='red')
    if len(r2_data) > 0:
        plt.hist(r2_data, bins=30, alpha=0.7, label=f'R2 (n={len(r2_data)})', color='blue')
    
    plt.xlabel('Trim Length (bp)')
    plt.ylabel('Frequency')
    plt.title(f'{sample_name} - Trim Length Distribution')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{sample_name}_lengths.png', dpi=150, bbox_inches='tight')
    plt.show()

def plot_scores(sample_name):
    """绘制分数分布"""
    score_file = f"{sample_name}_trimming_analysis_scores.csv"
    if not os.path.exists(score_file):
        print(f"File {score_file} not found")
        return
    
    data = pd.read_csv(score_file)
    scores = data['trim_score']
    
    plt.figure(figsize=(10, 6))
    
    plt.hist(scores, bins=50, alpha=0.7, color='green', edgecolor='black')
    plt.axvline(scores.mean(), color='red', linestyle='--', label=f'Mean: {scores.mean():.1f}')
    plt.axvline(scores.median(), color='orange', linestyle='--', label=f'Median: {scores.median():.1f}')
    
    plt.xlabel('Trimming Score')
    plt.ylabel('Frequency')
    plt.title(f'{sample_name} - Score Distribution')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{sample_name}_scores.png', dpi=150, bbox_inches='tight')
    plt.show()

def main():
    # 自动找到所有样本
    samples = []
    for file in os.listdir('.'):
        if file.endswith('_trimming_analysis.csv'):
            sample = file.replace('_trimming_analysis.csv', '')
            samples.append(sample)
    
    if not samples:
        print("No analysis CSV files found. Run the C++ program first.")
        return
    
    print(f"Found samples: {samples}")
    
    for sample in samples:
        print(f"\nProcessing {sample}...")
        plot_basic_stats(sample)
        plot_lengths(sample)
        plot_scores(sample)
        print(f"Generated plots for {sample}")

if __name__ == "__main__":
    main()