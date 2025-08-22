#include "mergeread.h"
#include "smarttrim.h"
#include <iostream>
#include <string>
#include <iomanip>

void debugMergeAndTrim(const std::string& r1_seq, const std::string& r2_seq, 
                       const std::string& test_name, SmartTrimmer& trimmer) {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "DEBUG: " << test_name << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    
    // 1. 显示原始序列
    std::cout << "\n1. ORIGINAL SEQUENCES:" << std::endl;
    std::cout << "R1 (" << r1_seq.length() << "bp): " << r1_seq << std::endl;
    std::cout << "R2 (" << r2_seq.length() << "bp): " << r2_seq << std::endl;
    
    // 2. 分析R2反向互补
    ReadMerger merger(10, 3, 0.3);
    std::string r2_rc = "";
    for (int i = r2_seq.length() - 1; i >= 0; i--) {
        char base = std::toupper(r2_seq[i]);
        switch (base) {
            case 'A': r2_rc += 'T'; break;
            case 'T': r2_rc += 'A'; break;
            case 'G': r2_rc += 'C'; break;
            case 'C': r2_rc += 'G'; break;
            case 'N': r2_rc += 'N'; break;
            default: r2_rc += 'N'; break;
        }
    }
    std::cout << "\n2. R2 REVERSE COMPLEMENT:" << std::endl;
    std::cout << "R2_RC (" << r2_rc.length() << "bp): " << r2_rc << std::endl;
    
    // 3. 尝试merge并显示详细信息
    std::cout << "\n3. MERGE ANALYSIS:" << std::endl;
    
    // 先分析overlap
    OverlapResult overlap = merger.analyzeOverlap(r1_seq, r2_seq);
    std::cout << "OVERLAP DETECTION:" << std::endl;
    std::cout << "  Overlapped: " << (overlap.overlapped ? "YES" : "NO") << std::endl;
    if (overlap.overlapped) {
        std::cout << "  Offset: " << overlap.offset << std::endl;
        std::cout << "  Overlap length: " << overlap.overlap_len << std::endl;
        std::cout << "  Diff count: " << overlap.diff_count << std::endl;
        std::cout << "  Diff percent: " << overlap.diff_percent << std::endl;
        
        // 显示overlap区域
        if (overlap.offset >= 0) {
            std::string r1_overlap = r1_seq.substr(overlap.offset, overlap.overlap_len);
            std::string r2_overlap = r2_rc.substr(0, overlap.overlap_len);
            std::cout << "  R1 overlap region: " << r1_overlap << std::endl;
            std::cout << "  R2 overlap region: " << r2_overlap << std::endl;
        } else {
            int r2_start = -(overlap.offset + 1);
            std::string r1_overlap = r1_seq.substr(0, overlap.overlap_len);
            std::string r2_overlap = r2_rc.substr(r2_start, overlap.overlap_len);
            std::cout << "  R1 overlap region: " << r1_overlap << std::endl;
            std::cout << "  R2 overlap region: " << r2_overlap << std::endl;
        }
    }
    
    MergeResult merge_result = merger.mergeReads(r1_seq, r2_seq);
    
    if (merge_result.merged) {
        std::cout << "MERGE SUCCESS" << std::endl;
        std::cout << "Merged (" << merge_result.sequence.length() << "bp): " << merge_result.sequence << std::endl;
        std::cout << "R1 contribution: " << merge_result.r1_bases << " bases" << std::endl;
        std::cout << "R2 contribution: " << merge_result.r2_bases << " bases" << std::endl;
        std::cout << "R1 start in merged: " << merge_result.r1_start_in_merged << std::endl;
        std::cout << "R2 end in merged: " << merge_result.r2_end_in_merged << std::endl;
        std::cout << "Overlap length: " << merge_result.overlap_length << std::endl;
        std::cout << "Merge offset: " << merge_result.merge_offset << std::endl;
        
        // 4. 分析merged序列的polyC tail
        std::cout << "\n4. POLYC TAIL ANALYSIS:" << std::endl;
        std::string merged_seq = merge_result.sequence;
        int merged_len = merged_seq.length();
        
        // 从末尾开始分析最多10个碱基
        std::cout << "Last 10 bases: ";
        int start_pos = std::max(0, merged_len - 10);
        for (int i = start_pos; i < merged_len; i++) {
            std::cout << merged_seq[i];
        }
        std::cout << std::endl;
        
        // 分析每个位置
        std::cout << "Position analysis (from end):" << std::endl;
        for (int i = 0; i < std::min(10, merged_len); i++) {
            int pos = merged_len - 1 - i;
            char base = std::toupper(merged_seq[pos]);
            std::cout << "  [" << i << "] pos=" << pos << " base=" << base;
            if (base == 'C' || base == 'N') {
                std::cout << " (polyC candidate)";
            } else {
                std::cout << " (NOT polyC)";
            }
            std::cout << std::endl;
        }
    } else {
        std::cout << " MERGE FAILED - will use fallback" << std::endl;
    }
    
    // 5. 运行trim算法并显示详细过程
    std::cout << "\n5. TRIM ALGORITHM:" << std::endl;
    TrimResult result = trimmer.findOptimalTrimPositions(r1_seq, r2_seq);
    
    if (result.is_valid) {
        std::cout << "TRIM SUCCESS" << std::endl;
        std::cout << "Score detail: " << result.score_detail << std::endl;
        std::cout << "R1 trim pos: " << result.r1_trim_pos << " (from end)" << std::endl;
        std::cout << "R2 trim pos: " << result.r2_trim_pos << " (from start)" << std::endl;
        
        // 6. 显示映射逻辑的详细过程
        std::cout << "\n6. MAPPING LOGIC DEBUG:" << std::endl;
        if (merge_result.merged) {
            int cut_length = result.r1_trim_pos + result.r2_trim_pos; // 假设这是总cut长度
            int merged_length = merge_result.sequence.length();
            int cut_start_pos = merged_length - cut_length;
            
            std::cout << "Cut length from end: " << cut_length << std::endl;
            std::cout << "Cut start position: " << cut_start_pos << std::endl;
            std::cout << "R1 start in merged: " << merge_result.r1_start_in_merged << std::endl;
            std::cout << "R2 end in merged: " << merge_result.r2_end_in_merged << std::endl;
            
            if (merge_result.r1_start_in_merged >= 0) {
                if (cut_start_pos >= merge_result.r1_start_in_merged) {
                    std::cout << "→ Cut region is in R1 part" << std::endl;
                } else if (cut_start_pos + cut_length <= merge_result.r2_end_in_merged) {
                    std::cout << "→ Cut region is in R2 part" << std::endl;
                } else {
                    std::cout << "→ Cut region spans R2 and R1" << std::endl;
                }
            }
        }
        
        // 7. 显示trim后的结果
        std::cout << "\n7. TRIM RESULTS:" << std::endl;
        auto trimmed = trimmer.trimReads(r1_seq, r2_seq, result);
        std::cout << "R1 before: " << r1_seq << " (" << r1_seq.length() << "bp)" << std::endl;
        std::cout << "R1 after:  " << trimmed.first << " (" << trimmed.first.length() << "bp)" << std::endl;
        std::cout << "R2 before: " << r2_seq << " (" << r2_seq.length() << "bp)" << std::endl;
        std::cout << "R2 after:  " << trimmed.second << " (" << trimmed.second.length() << "bp)" << std::endl;
        
        // 8. 验证trim是否合理
        std::cout << "\n8. TRIM VALIDATION:" << std::endl;
        if (result.r1_trim_pos > 0) {
            std::string r1_trimmed_part = r1_seq.substr(r1_seq.length() - result.r1_trim_pos);
            std::cout << "R1 trimmed part: " << r1_trimmed_part << std::endl;
            
            int c_count = 0;
            for (char base : r1_trimmed_part) {
                if (std::toupper(base) == 'C' || std::toupper(base) == 'N') c_count++;
            }
            double c_ratio = (double)c_count / r1_trimmed_part.length();
            std::cout << "R1 trimmed C+N ratio: " << std::fixed << std::setprecision(2) << c_ratio * 100 << "%" << std::endl;
        }
        
        if (result.r2_trim_pos > 0) {
            std::string r2_trimmed_part = r2_seq.substr(0, result.r2_trim_pos);
            std::cout << "R2 trimmed part: " << r2_trimmed_part << std::endl;
            
            int g_count = 0;
            for (char base : r2_trimmed_part) {
                if (std::toupper(base) == 'G' || std::toupper(base) == 'N') g_count++;
            }
            double g_ratio = (double)g_count / r2_trimmed_part.length();
            std::cout << "R2 trimmed G+N ratio: " << std::fixed << std::setprecision(2) << g_ratio * 100 << "%" << std::endl;
        }
        
    } else {
        std::cout << " TRIM FAILED" << std::endl;
    }
    
    std::cout << "\nDEBUG END\n" << std::endl;
}

int main() {
    std::cout << "=== DETAILED DEBUG TEST ===" << std::endl;
    
    // 使用默认参数
    TrimParams params(25, 0.0);
    SmartTrimmer trimmer(params);
    
    // 测试几个关键案例
    struct TestCase {
        std::string r1, r2, name;
    };
    
    std::vector<TestCase> tests = {
        {"ATCGATCGATCGATCGATCCCCCCC", "GGGGGATCGATCGATCGATCGA", "Basic PolyC"},
        {"ATCGATCGATCCCCCCCCCCCCC", "GGGGGGGGGGATCGATCGATC", "Long PolyC"},
        {"ATCGATCGATCGATCGATCGCACAC", "GTGTGATCGATCGATCGATCGA", "Mixed bases"}
    };
    
    // 只测试第一个案例，详细分析
    debugMergeAndTrim(tests[0].r1, tests[0].r2, tests[0].name, trimmer);
    
    return 0;
}
