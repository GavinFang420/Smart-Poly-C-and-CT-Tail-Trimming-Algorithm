#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def analyze_parameters():
    """分析参数优化结果"""
    
    # 读取参数优化结果
    try:
        df = pd.read_csv('parameter_optimization_results.csv')
    except FileNotFoundError:
        print("parameter_optimization_results.csv not found. Run C++ analysis first.")
        return
    
    print(f"Loaded {len(df)} parameter configurations")
    
    # 1. 最佳参数配置
    print("\nTOP 5 PARAMETER CONFIGURATIONS:")
    print("=" * 80)
    top5 = df.head(5)
    for i, row in top5.iterrows():
        print(f"Rank {row['rank']}: Window={row['window_size']}, InitScore={row['initial_score']:.1f}, "
              f"CScore={row['c_score']:.1f}, Penalty={row['penalty_score']:.1f}")
        print(f"           TrimRate={row['trim_rate']*100:.1f}%, AvgScore={row['avg_score']:.1f}")
        print()
    
    # 2. 参数vs性能关系图
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('Parameter Optimization Analysis', fontsize=16)
    
    # Window size vs performance
    axes[0,0].scatter(df['window_size'], df['avg_score'], alpha=0.6, c=df['trim_rate'], cmap='viridis')
    axes[0,0].set_xlabel('Window Size')
    axes[0,0].set_ylabel('Average Score')
    axes[0,0].set_title('Window Size vs Score')
    
    # Initial score vs performance
    axes[0,1].scatter(df['initial_score'], df['avg_score'], alpha=0.6, c=df['trim_rate'], cmap='viridis')
    axes[0,1].set_xlabel('Initial Score')
    axes[0,1].set_ylabel('Average Score')
    axes[0,1].set_title('Initial Score vs Performance')
    
    # Penalty score vs performance
    axes[0,2].scatter(df['penalty_score'], df['avg_score'], alpha=0.6, c=df['trim_rate'], cmap='viridis')
    axes[0,2].set_xlabel('Penalty Score')
    axes[0,2].set_ylabel('Average Score')
    axes[0,2].set_title('Penalty vs Performance')
    
    # Trim rate distribution
    axes[1,0].hist(df['trim_rate'] * 100, bins=30, alpha=0.7, color='skyblue')
    axes[1,0].set_xlabel('Trim Rate (%)')
    axes[1,0].set_ylabel('Frequency')
    axes[1,0].set_title('Trim Rate Distribution')
    
    # Score distribution
    axes[1,1].hist(df['avg_score'], bins=30, alpha=0.7, color='lightgreen')
    axes[1,1].set_xlabel('Average Score')
    axes[1,1].set_ylabel('Frequency')
    axes[1,1].set_title('Score Distribution')
    
    # Trim rate vs Score
    scatter = axes[1,2].scatter(df['trim_rate'] * 100, df['avg_score'], 
                               alpha=0.6, c=df['window_size'], cmap='plasma')
    axes[1,2].set_xlabel('Trim Rate (%)')
    axes[1,2].set_ylabel('Average Score')
    axes[1,2].set_title('Trim Rate vs Score')
    plt.colorbar(scatter, ax=axes[1,2], label='Window Size')
    
    plt.tight_layout()
    plt.savefig('parameter_optimization_analysis.png', dpi=150, bbox_inches='tight')
    plt.show()
    
    # 3. 推荐最佳参数
    print("\nRECOMMENDED PARAMETERS:")
    print("=" * 50)
    
    best = df.iloc[0]
    print(f"BEST OVERALL (Rank 1):")
    print(f"  Window Size: {best['window_size']}")
    print(f"  Initial Score: {best['initial_score']}")
    print(f"  C Score: {best['c_score']}")
    print(f"  Penalty Score: {best['penalty_score']}")
    print(f"  Performance: {best['trim_rate']*100:.1f}% trim rate, {best['avg_score']:.1f} avg score")
    
    # 寻找平衡的参数（中等trim rate但高分数）
    balanced = df[(df['trim_rate'] > 0.3) & (df['trim_rate'] < 0.7)].iloc[0]
    print(f"\nBALANCED OPTION (Moderate trim rate):")
    print(f"  Window Size: {balanced['window_size']}")
    print(f"  Initial Score: {balanced['initial_score']}")
    print(f"  C Score: {balanced['c_score']}")
    print(f"  Penalty Score: {balanced['penalty_score']}")
    print(f"  Performance: {balanced['trim_rate']*100:.1f}% trim rate, {balanced['avg_score']:.1f} avg score")

def suggest_next_steps():
    """根据结果推荐下一步"""
    print("\nNEXT STEPS:")
    print("=" * 40)
    print("1. Use the BEST OVERALL parameters for production runs")
    print("2. If trim rate is too high/low, try the BALANCED OPTION")
    print("3. Fine-tune around the best parameters:")
    print("   - Adjust window_size ±2")
    print("   - Adjust penalty_score ±0.5")
    print("4. Test on different samples to validate consistency")
    print("5. Monitor false positive/negative rates in real data")

if __name__ == "__main__":
    analyze_parameters()
    suggest_next_steps()