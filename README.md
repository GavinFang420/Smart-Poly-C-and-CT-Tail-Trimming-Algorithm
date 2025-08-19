# 项目介绍 - 开发中
本项目主要解决Poly C tail污染相关的问题，CT tail不在初期考虑范围内

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






## 适用人群与适用情境（Intended Users & Indications）

- 本项目当前适用于 WGBS 单链建库流程中，由文库化学步骤在片段 3′ 端引入的poly-C（在 R2 上通常表现为 poly-G） 合成尾和poly-CT 合成尾，特别针对 xGen Methyl-Sequencing DNA Library Prep Kit and Adaptase | IDT 试剂盒。
- 触发前提：R2 起始端存在 C/G 尾迹象（FastQC 的 per-base content 在 R2 前 5–12bp 出现 G 峰；或“R2 开头连续 G 的长度”呈集中分布，典型中位约 8–10bp）。若未观察到这些特征，算法会低触发或不触发，不会破坏数据。


# 适用人群
- WGBS 流水线维护者/分析员：需在**不合并双端读（PE）**前提下去除 R1 R2合成尾
- 低输入与临床相关样本的研究者（cfDNA/液体活检）：样本不可复采，片段短且常无重叠，对过度剪切与误去重高度敏感
- 重视可重复性与方法学呈现的研究团队：应用自适应 R2 截尾后，fastp 的 After-filtering 面板中，R2 的 G 尖峰消失；K-mer 异常热区清除，R1 指标保持稳定，证明仅去除了 R2 5′ 的合成尾而未伤及真实序列

# 适用判据（供参考，因为此方法运行效率可能比在R2头部，R1尾部直接剪切15碱基低）
- R2 的 1–10 bp 与 25–50 bp 的 G%/C% 差值 ≥ 10 个百分点；
- R1 的 倒数 1–10 bp 相对 倒数 25–50 bp 的 C/G% 差值 ≥6–8 个百分点，或 R1 末端连续 C/G 的中位 ≥4 bp → 可开启 R1 截尾（保守参数）。若无上述迹象，R1 不截更稳。

# 不适用/收益有限
- 无合成尾特征的文库；
- 双链酶转换（如 EM-seq）等非 Adaptase类流程；

## 重要文件一览
？？？
？？？
？？？？
？？？
？？？？

## 引用/参照/致谢
FASTQ处理模块基于 [fastp](https://github.com/OpenGene/fastp) 修改
感谢纪博帮助提供的测试数据
牢方当时的idea

Copyright (c) 2017 OpenGene - MIT License
