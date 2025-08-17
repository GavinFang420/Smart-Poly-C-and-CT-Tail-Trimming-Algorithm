# 项目介绍 - 开发中
本项目主要解决Poly C tail相关的问题，CT tail不在初期考虑范围内

目标是在 WGBS / 单链建库流程中，开发一个更加智能的 tail trimming 算法，避免传统方法带来的过度剪切或不足剪切问题。

同时，本方法可以避免 因不同 C tail 长度而误去重 的情况。
例如：148bp + 2C tail 和 145bp + 5C tail，虽然真实片段不同，但传统去重算法可能误判为重复。

## 问题与现状
针对单链建库，目前fastq的tail trimming有两个问题（2025/08/15未找到相关改进算法paper）：
1. Under Trimming - Some R2 heads that are longer than 15bp (the fastq recommended value), resulting in some short residues.
3. Over Trimming - Valuable R1 tail is trimmed along with some overkill R2 head. The waste is huge, especially when R1 and R2 do not have overlap (sum >= 300bp, meaning 15 bp info in R1 is wasted regardless of R2).

## 基本假设与原理
### 假设 - C/CT于R2末端富集 + 长度平均8bp
- PolyC/CT 尾富集：C/CT 序列主要集中在 R2 末端，平均长度约 8bp。
- C 碱基稀有性：经亚硫酸氢盐处理后，C 在基因组中极为稀有，仅占约 0.4%（≈99.5% 的 C 被转化为 T）。
- PolyC/CT 富集性：人工引入的 polyC/CT 段具有极高 C 含量（尾长度 ≥1bp，呈正态分布，均值 8bp）。

### 原理 - C tail位置衰减积分算法 应用于 R1+R2 Merged Sequence 
采用 位置衰减积分（Distance-Decay Integration）算法，在 R1+R2 拼接序列上识别 polyC 区域。

Scoring 规则：
- C 出现 → 加分
- A/G/T 出现 → 扣分（负分）

Distance Decay 原理：
- 越靠近 R2 末端，C 的分数权重越大，更可能属于 polyC 尾 → 应该被剪除。
- 越靠近 R1 尾部，权重越低，更可能是真实序列 → 应该被保留。

## 算法流程
1. 复制 & 合并：复制一份原始 reads → 构造 S = R2 头 K bp + R1 尾 K bp。
2. 智能打分：在 S 上用调参后的距离衰减打分，识别最优 polyC 区域。
3. 记录 index：在合并序列中找到 polyC 区域的起止坐标。
4. 映射回原始 reads：将 index 映射到原始 R1、R2 上的对应位置。
5. 精确剪切：直接对 原始 R1、R2 进行裁剪（保持 pair 结构，输出未合并的 reads）。
6. 验证（仅 CT 版本时启用）：校验剪切区域的 C/T 比例是否接近 85%C : 15%T。

## 使用方法
目前我们使用 fastp 先做 adapter trimming，随后运行本项目的 polyC tail trimming 模块。

## Features
- Preserves read pairing: Outputs separate R1/R2 files, not merged
- Position-aware: Higher weights for sequence ends
- Ratio validation: Confirms trimmed regions match expected C/T composition (for CT tailing Only, will develop)
- R1 preservation: Saves valuable R1 data that traditional tools waste (Single Chain)
- Fast: Theoretical Time Complexity O(n)

## 引用/参照/致谢
FASTQ处理模块基于 [fastp](https://github.com/OpenGene/fastp) 修改

Copyright (c) 2017 OpenGene - MIT License
