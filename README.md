# 项目介绍
本项目主要解决Poly C tail相关的问题，CT tail不在初期考虑范围内
## 问题与现状
针对单链建库，目前fastq的tail trimming有两个问题（2025/08/15未找到相关改进算法paper）：
1. Under Trimming - Some R2 heads that are longer than 15bp (the fastq recommended value), resulting in some short residues.
3. Over Trimming - Valuable R1 tail is trimmed along with some overkill R2 head. The waste is huge, especially when R1 and R2 do not have overlap (sum >= 300bp, meaning 15 bp info in R1 is wasted regardless of R2).

## 基本假设与原理
### 假设 - C/CT于R2末端富集 + 长度平均10bp
C碱基稀有性: 盐转化后的C碱基十分稀有，只占序列的0.4% （99.5%的C都转化为了T）
PolyC/CT尾C/CT富集性: 人工引入的polyC/CT段具有极高的C含量 （Tail长度最小为1，以10为中心呈正态分布）
oligo长度特性: CT oligo通常长度在10bp左右，集中在序列末端

### 原理 - C tail位置衰减积分算法 应用于 R1+R2 Merged Sequence 
Formula:
  threshold = POGRESSIVE_SUM(Position (start from R2) * Scoring) + Base （Will Tune）
  if threshold < 0 cut; else next
Scoring:
  C appearance - Add Unit Score
  A,G,T appearance - Penalize （Negative Score）
Distance Decay原理: 越靠近序列中段，越不应该被切除 （权重越小）
  i.e. 靠近R2 Head部分的C更应该被切除，靠近R1尾部的C相对其他的C更应该被切除

## 算法流程
1. Merge Read1 + Read2 → Focus on last 30bp (WILL TUNE: 20/21/22/23/25/30)
2. Intelligent scoring → Identify optimal polyC/CT regions (Distance-Decay Scoring Algo)
3. Position mapping → Find corresponding positions in original R1/R2 (Use recorded index to locate positions)
4. Precise trimming → Trim original R1/R2 separately (maintain pair structure, will not return merged Reads)
5. Validation → Verify 85%C + 15%T ratio in trimmed regions (Will Tune, ONLY FOR CT tailing OPTION)

## Features
Preserves read pairing: Outputs separate R1/R2 files, not merged
Position-aware: Higher weights for sequence ends
Ratio validation: Confirms trimmed regions match expected C/T composition (for CT tailing Only, will develop)
R1 preservation: Saves valuable R1 data that traditional tools waste (Single Chain)
Fast: Theoretical Time Complexity O(n)
