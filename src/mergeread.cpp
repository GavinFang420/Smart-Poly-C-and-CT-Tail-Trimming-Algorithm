#include "mergeread.h"
#include <algorithm>
#include <cmath>
#include <iostream>

// Helper function to get reverse complement
std::string ReadMerger::reverseComplement(const std::string& seq) {
    std::string rc;
    rc.reserve(seq.length());
    
    for (int i = seq.length() - 1; i >= 0; i--) {
        rc += getComplementBase(seq[i]);
    }
    
    return rc;
}

bool ReadMerger::isValidBase(char base) {
    char upper_base = std::toupper(base);
    return (upper_base == 'A' || upper_base == 'T' || 
            upper_base == 'G' || upper_base == 'C' || upper_base == 'N');
}

char ReadMerger::getComplementBase(char base) {
    switch (std::toupper(base)) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'G': return 'C';
        case 'C': return 'G';
        case 'N': return 'N';  // N remains N, but could be context-dependent
        default: return 'N';
    }
}

// Modified overlap detection for PolyC tail processing
// R2 starts with GG tail (complement of polyC tail)
// Logic: complement R2 first, then merge with R1
// Result: R1 head at 3' end, R2 head at 5' end
OverlapResult ReadMerger::detectOverlap(const std::string& r1_seq, 
                                       const std::string& r2_seq,
                                       const std::string& r1_qual,
                                       const std::string& r2_qual) {
    OverlapResult result;
    
    if (r1_seq.empty() || r2_seq.empty()) {
        return result;
    }
    
    // For PolyC tail processing:
    // 1. R2 starts with GG tail (complement of polyC tail)
    // 2. Get reverse complement of R2
    std::string r2_rc = reverseComplement(r2_seq);
    std::string r2_qual_rev = r2_qual;
    if (!r2_qual_rev.empty()) {
        std::reverse(r2_qual_rev.begin(), r2_qual_rev.end());
    }
    
    // 3. Now we need to find overlap between R1 and R2_RC
    // The expected structure after merge should be: R2_head(5') + overlap + R1_head(3')
    // This means we're looking for R1's 3' end overlapping with R2_RC's 5' end
    
    int r1_len = r1_seq.length();
    int r2_rc_len = r2_rc.length();
    
    int best_offset = -1;
    int best_overlap_len = 0;
    int best_diff_count = 999999;
    
    // Look for overlap: R1's tail overlapping with R2_RC's head
    // R2_RC: [head]----------[tail]
    // R1:           [head]----------[tail]
    //              ^overlap region^
    
    for (int r1_start = 0; r1_start <= r1_len - min_overlap_len; r1_start++) {
        int max_possible_overlap = std::min(r1_len - r1_start, r2_rc_len);
        
        for (int overlap_len = min_overlap_len; overlap_len <= max_possible_overlap; overlap_len++) {
            // Check if R1 suffix matches R2_RC prefix
            std::string r1_region = r1_seq.substr(r1_start, overlap_len);
            std::string r2_region = r2_rc.substr(0, overlap_len);
            
            // Count mismatches
            int diff_count = 0;
            for (int i = 0; i < overlap_len; i++) {
                if (std::toupper(r1_region[i]) != std::toupper(r2_region[i])) {
                    diff_count++;
                }
            }
            
            double diff_percent = (double)diff_count / overlap_len;
            
            // Check if this overlap meets criteria
            bool valid_overlap = (diff_count <= max_diff_count) && 
                               (diff_percent <= max_diff_percent);
            
            if (valid_overlap && (diff_count < best_diff_count || 
                (diff_count == best_diff_count && overlap_len > best_overlap_len))) {
                best_offset = r1_start;  // Position in R1 where overlap starts
                best_overlap_len = overlap_len;
                best_diff_count = diff_count;
            }
        }
    }
    
    // Also check the reverse case: R2_RC extends beyond R1
    // R2_RC: [head]----------------[tail]
    // R1:          [head]----[tail]
    //              ^overlap^
    
    for (int r2_start = 0; r2_start <= r2_rc_len - min_overlap_len; r2_start++) {
        int max_possible_overlap = std::min(r2_rc_len - r2_start, r1_len);
        
        for (int overlap_len = min_overlap_len; overlap_len <= max_possible_overlap; overlap_len++) {
            // Check if R1 matches R2_RC suffix starting from r2_start
            std::string r1_region = r1_seq.substr(0, overlap_len);
            std::string r2_region = r2_rc.substr(r2_start, overlap_len);
            
            int diff_count = 0;
            for (int i = 0; i < overlap_len; i++) {
                char r1_base = std::toupper(r1_region[i]);
                char r2_base = std::toupper(r2_region[i]);
                
                // Special handling for N bases
                if (r1_base == 'N' || r2_base == 'N') {
                    if (!r1_qual.empty() && !r2_qual_rev.empty()) {
                        char r1_q = r1_qual[i];
                        char r2_q = r2_qual_rev[r2_start + i];
                        
                        if (r1_q < '#' || r2_q < '#') { // Phred < 3
                            continue;
                        }
                    }
                    diff_count += 0.5; // Half penalty for N
                } else if (r1_base != r2_base) {
                    diff_count++;
                }
            }
            
            double diff_percent = (double)diff_count / overlap_len;
            
            bool valid_overlap = (diff_count <= max_diff_count) && 
                               (diff_percent <= max_diff_percent);
            
            if (valid_overlap && (diff_count < best_diff_count || 
                (diff_count == best_diff_count && overlap_len > best_overlap_len))) {
                best_offset = -r2_start - 1; // Negative to indicate R2_RC extends beyond R1
                best_overlap_len = overlap_len;
                best_diff_count = diff_count;
            }
        }
    }
    
    if (best_offset != -1 && best_overlap_len >= min_overlap_len) {
        result.overlapped = true;
        result.offset = best_offset;
        result.overlap_len = best_overlap_len;
        result.diff_count = best_diff_count;
        result.diff_percent = (double)best_diff_count / best_overlap_len;
    }
    
    return result;
}

// Main merge function - sequence only
MergeResult ReadMerger::mergeReads(const std::string& r1_seq, const std::string& r2_seq) {
    return mergeReads(r1_seq, r2_seq, "", "");
}

// Modified merge function for PolyC tail processing with alignment info
// Result structure: R2_head(5') + overlap + R1_head(3')
MergeResult ReadMerger::mergeReads(const std::string& r1_seq, const std::string& r2_seq,
                                  const std::string& r1_qual, const std::string& r2_qual) {
    MergeResult result;
    
    OverlapResult overlap = detectOverlap(r1_seq, r2_seq, r1_qual, r2_qual);
    
    if (!overlap.overlapped) {
        return result; // No valid overlap found
    }
    
    // Prepare reverse complement of R2 (since R2 has GG tail at start)
    std::string r2_rc = reverseComplement(r2_seq);
    std::string r2_qual_rev = r2_qual;
    if (!r2_qual_rev.empty()) {
        std::reverse(r2_qual_rev.begin(), r2_qual_rev.end());
    }
    
    std::string merged_seq;
    std::string merged_qual;
    
    if (overlap.offset >= 0) {
        // Case 1: Standard overlap - R1 suffix overlaps with R2_RC prefix
        // 修正：不要丢弃R1的头部！
        // Structure: R1_prefix + [overlap region] + R1_suffix
        
        // 添加R1的前缀部分（overlap之前的部分）
        if (overlap.offset > 0) {
            merged_seq = r1_seq.substr(0, overlap.offset);
            if (!r1_qual.empty()) {
                merged_qual = r1_qual.substr(0, overlap.offset);
            }
        }
        
        // Add overlapped region (merge qualities, prefer higher quality base)
        for (int i = 0; i < overlap.overlap_len; i++) {
            char r1_base = r1_seq[overlap.offset + i];
            char r2_base = r2_rc[i];
            
            if (!r1_qual.empty() && !r2_qual_rev.empty()) {
                char r1_q = r1_qual[overlap.offset + i];
                char r2_q = r2_qual_rev[i];
                
                // Special handling for N bases - prefer non-N base
                if (std::toupper(r1_base) == 'N' && std::toupper(r2_base) != 'N') {
                    merged_seq += r2_base;
                    merged_qual += r2_q;
                } else if (std::toupper(r2_base) == 'N' && std::toupper(r1_base) != 'N') {
                    merged_seq += r1_base;
                    merged_qual += r1_q;
                } else if (std::toupper(r1_base) == 'N' && std::toupper(r2_base) == 'N') {
                    // Both are N, check if this is likely polyC region
                    // If in the tail region of merged sequence, assume C
                    if (i >= overlap.overlap_len - 5) { // Last 5bp of overlap
                        merged_seq += 'C'; // Assume N is miscalled C in polyC region
                        merged_qual += std::max(r1_q, r2_q);
                    } else {
                        merged_seq += 'N';
                        merged_qual += std::max(r1_q, r2_q);
                    }
                } else {
                    // Choose base with higher quality
                    if (r1_q >= r2_q) {
                        merged_seq += r1_base;
                        merged_qual += r1_q;
                    } else {
                        merged_seq += r2_base;
                        merged_qual += r2_q;
                    }
                }
            } else {
                // No quality scores available, prefer non-N base
                if (std::toupper(r1_base) == 'N' && std::toupper(r2_base) != 'N') {
                    merged_seq += r2_base;
                } else if (std::toupper(r2_base) == 'N' && std::toupper(r1_base) != 'N') {
                    merged_seq += r1_base;
                } else {
                    merged_seq += r1_base; // Default to R1
                }
                if (!r1_qual.empty()) {
                    merged_qual += r1_qual[overlap.offset + i];
                }
            }
        }
        
        // Add R1 suffix (after overlap region) 
        int r1_suffix_start = overlap.offset + overlap.overlap_len;
        if (r1_suffix_start < (int)r1_seq.length()) {
            merged_seq += r1_seq.substr(r1_suffix_start);
            if (!r1_qual.empty()) {
                merged_qual += r1_qual.substr(r1_suffix_start);
            }
        }
        
        // Add R2_RC suffix (parts that don't overlap with R1)
        if (overlap.overlap_len < (int)r2_rc.length()) {
            merged_seq += r2_rc.substr(overlap.overlap_len);
            if (!r2_qual_rev.empty()) {
                merged_qual += r2_qual_rev.substr(overlap.overlap_len);
            }
        }
        
        result.r1_bases = r1_seq.length();
        result.r2_bases = r2_rc.length() - overlap.overlap_len;
        
        // 修正对齐信息计算
        result.r1_start_in_merged = 0;
        result.r2_end_in_merged = overlap.offset + overlap.overlap_len + (r2_rc.length() - overlap.overlap_len);
        result.overlap_length = overlap.overlap_len;
        result.merge_offset = overlap.offset;
        
    } else {
        // Case 2: R2_RC extends beyond R1
        // R2_RC: [prefix]----[overlap]----[suffix]
        // R1:               [overlap]----[suffix]
        
        int r2_offset = -(overlap.offset + 1);
        
        // Add R2_RC prefix
        merged_seq = r2_rc.substr(0, r2_offset);
        if (!r2_qual_rev.empty()) {
            merged_qual = r2_qual_rev.substr(0, r2_offset);
        }
        
        // Add overlapped region
        for (int i = 0; i < overlap.overlap_len; i++) {
            char r1_base = r1_seq[i];
            char r2_base = r2_rc[r2_offset + i];
            
            if (!r1_qual.empty() && !r2_qual_rev.empty()) {
                char r1_q = r1_qual[i];
                char r2_q = r2_qual_rev[r2_offset + i];
                
                // Special handling for N bases - prefer non-N base
                if (std::toupper(r1_base) == 'N' && std::toupper(r2_base) != 'N') {
                    merged_seq += r2_base;
                    merged_qual += r2_q;
                } else if (std::toupper(r2_base) == 'N' && std::toupper(r1_base) != 'N') {
                    merged_seq += r1_base;
                    merged_qual += r1_q;
                } else if (std::toupper(r1_base) == 'N' && std::toupper(r2_base) == 'N') {
                    // Both are N, check if this is likely polyC region
                    if (i >= overlap.overlap_len - 5) { // Last 5bp of overlap
                        merged_seq += 'C'; // Assume N is miscalled C in polyC region
                        merged_qual += std::max(r1_q, r2_q);
                    } else {
                        merged_seq += 'N';
                        merged_qual += std::max(r1_q, r2_q);
                    }
                } else {
                    // Choose base with higher quality
                    if (r1_q >= r2_q) {
                        merged_seq += r1_base;
                        merged_qual += r1_q;
                    } else {
                        merged_seq += r2_base;
                        merged_qual += r2_q;
                    }
                }
            } else {
                // No quality scores, prefer non-N base
                if (std::toupper(r1_base) == 'N' && std::toupper(r2_base) != 'N') {
                    merged_seq += r2_base;
                } else if (std::toupper(r2_base) == 'N' && std::toupper(r1_base) != 'N') {
                    merged_seq += r1_base;
                } else {
                    merged_seq += r1_base; // Default to R1
                }
                if (!r1_qual.empty()) {
                    merged_qual += r1_qual[i];
                }
            }
        }
        
        // Add R1 suffix
        if (overlap.overlap_len < (int)r1_seq.length()) {
            merged_seq += r1_seq.substr(overlap.overlap_len);
            if (!r1_qual.empty()) {
                merged_qual += r1_qual.substr(overlap.overlap_len);
            }
        }
        
        // Add remaining R2_RC suffix if any
        int r2_remaining_start = r2_offset + overlap.overlap_len;
        if (r2_remaining_start < (int)r2_rc.length()) {
            merged_seq += r2_rc.substr(r2_remaining_start);
            if (!r2_qual_rev.empty()) {
                merged_qual += r2_qual_rev.substr(r2_remaining_start);
            }
        }
        
        result.r1_bases = r1_seq.length();
        result.r2_bases = r2_offset + overlap.overlap_len;
        
        // 计算对齐信息
        result.r2_end_in_merged = r2_offset + overlap.overlap_len;
        result.r1_start_in_merged = r2_offset;
        result.overlap_length = overlap.overlap_len;
        result.merge_offset = r2_offset;
    }
    
    result.merged = true;
    result.sequence = merged_seq;
    result.quality = merged_qual;
    
    return result;
}

OverlapResult ReadMerger::analyzeOverlap(const std::string& r1_seq, const std::string& r2_seq) {
    return detectOverlap(r1_seq, r2_seq);
}

std::string ReadMerger::getSubsequence(const std::string& seq, int start, int length) {
    if (start < 0 || start >= (int)seq.length()) {
        return "";
    }
    
    int actual_length = std::min(length, (int)seq.length() - start);
    return seq.substr(start, actual_length);
}

bool ReadMerger::hasValidOverlap(const OverlapResult& result, int min_len, 
                                int max_diff, double max_pct) {
    return result.overlapped && 
           result.overlap_len >= min_len &&
           result.diff_count <= max_diff &&
           result.diff_percent <= max_pct;
}

// Utility functions implementation
namespace MergeUtils {
    double qualToProb(char qual) {
        int q = qual - 33; // Phred+33 encoding
        return std::pow(10.0, -q / 10.0);
    }
    
    char probToQual(double prob) {
        if (prob <= 0) return 126; // Max quality
        double q = -10.0 * std::log10(prob);
        int qual_int = std::min(93, std::max(0, (int)(q + 0.5))); // Clamp to valid range
        return qual_int + 33;
    }
    
    char mergeQualityScores(char q1, char q2) {
        // Take the higher quality score (lower error probability)
        double p1 = qualToProb(q1);
        double p2 = qualToProb(q2);
        
        // Use the base with higher quality (lower error probability)
        // For consensus, we could use: combined_prob = p1 * p2 / (p1 * p2 + (1-p1) * (1-p2))
        // But for simplicity, just take the better quality
        return (p1 < p2) ? q1 : q2;
    }
    
    bool isValidDNASequence(const std::string& seq) {
        for (char base : seq) {
            char upper_base = std::toupper(base);
            if (upper_base != 'A' && upper_base != 'T' && 
                upper_base != 'G' && upper_base != 'C' && upper_base != 'N') {
                return false;
            }
        }
        return true;
    }
    
    SequenceStats getSequenceStats(const std::string& seq) {
        SequenceStats stats;
        stats.length = seq.length();
        stats.gc_count = 0;
        stats.n_count = 0;
        stats.polyG_start_length = 0;
        stats.polyC_end_length = 0;
        
        for (char base : seq) {
            char upper_base = std::toupper(base);
            if (upper_base == 'G' || upper_base == 'C') {
                stats.gc_count++;
            } else if (upper_base == 'N') {
                stats.n_count++;
            }
        }
        
        // Count polyG at start
        for (size_t i = 0; i < seq.length(); i++) {
            char upper_base = std::toupper(seq[i]);
            if (upper_base == 'G' || upper_base == 'N') { // N might be miscalled G
                stats.polyG_start_length++;
            } else {
                break;
            }
        }
        
        // Count polyC at end
        for (int i = seq.length() - 1; i >= 0; i--) {
            char upper_base = std::toupper(seq[i]);
            if (upper_base == 'C' || upper_base == 'N') { // N might be miscalled C
                stats.polyC_end_length++;
            } else {
                break;
            }
        }
        
        stats.gc_percent = stats.length > 0 ? 
            (double)stats.gc_count / stats.length * 100.0 : 0.0;
        
        return stats;
    }
    
    // Reverse complement utility (exposed for external use)
    std::string reverseComplement(const std::string& seq) {
        std::string rc;
        rc.reserve(seq.length());
        
        for (int i = seq.length() - 1; i >= 0; i--) {
            char base = std::toupper(seq[i]);
            switch (base) {
                case 'A': rc += 'T'; break;
                case 'T': rc += 'A'; break;
                case 'G': rc += 'C'; break;
                case 'C': rc += 'G'; break;
                case 'N': rc += 'N'; break;
                default: rc += 'N'; break;
            }
        }
        
        return rc;
    }
    
    // PolyC tail specific utilities
    bool hasPolyTail(const std::string& seq, char base, int min_length, bool count_N_as_base) {
        if (seq.length() < min_length) return false;
        
        char upper_base = std::toupper(base);
        int count = 0;
        
        for (int i = seq.length() - 1; i >= 0; i--) {
            char seq_base = std::toupper(seq[i]);
            if (seq_base == upper_base || (count_N_as_base && seq_base == 'N')) {
                count++;
                if (count >= min_length) return true;
            } else {
                break;
            }
        }
        
        return false;
    }
    
    int getPolyTailLength(const std::string& seq, char base, bool count_N_as_base) {
        char upper_base = std::toupper(base);
        int count = 0;
        
        for (int i = seq.length() - 1; i >= 0; i--) {
            char seq_base = std::toupper(seq[i]);
            if (seq_base == upper_base || (count_N_as_base && seq_base == 'N')) {
                count++;
            } else {
                break;
            }
        }
        
        return count;
    }
    
    std::string trimPolyTail(const std::string& seq, char base, int min_length, bool count_N_as_base) {
        int tail_length = getPolyTailLength(seq, base, count_N_as_base);
        if (tail_length >= min_length) {
            return seq.substr(0, seq.length() - tail_length);
        }
        return seq;
    }
    
    // Quality-aware versions
    bool hasPolyTailWithQuality(const std::string& seq, const std::string& qual,
                               char base, int min_length, char min_quality, bool count_N_as_base) {
        if (seq.length() < min_length || seq.length() != qual.length()) return false;
        
        char upper_base = std::toupper(base);
        int count = 0;
        
        for (int i = seq.length() - 1; i >= 0; i--) {
            char seq_base = std::toupper(seq[i]);
            char q = qual[i];
            
            if (q >= min_quality && (seq_base == upper_base || (count_N_as_base && seq_base == 'N'))) {
                count++;
                if (count >= min_length) return true;
            } else {
                break;
            }
        }
        
        return false;
    }
    
    int getPolyTailLengthWithQuality(const std::string& seq, const std::string& qual,
                                    char base, char min_quality, bool count_N_as_base) {
        if (seq.length() != qual.length()) return 0;
        
        char upper_base = std::toupper(base);
        int count = 0;
        
        for (int i = seq.length() - 1; i >= 0; i--) {
            char seq_base = std::toupper(seq[i]);
            char q = qual[i];
            
            if (q >= min_quality && (seq_base == upper_base || (count_N_as_base && seq_base == 'N'))) {
                count++;
            } else {
                break;
            }
        }
        
        return count;
    }
    
    std::string trimPolyTailWithQuality(const std::string& seq, const std::string& qual,
                                       char base, int min_length, char min_quality, bool count_N_as_base) {
        int tail_length = getPolyTailLengthWithQuality(seq, qual, base, min_quality, count_N_as_base);
        if (tail_length >= min_length) {
            return seq.substr(0, seq.length() - tail_length);
        }
        return seq;
    }
    
    // Quality validation
    bool isValidQualityString(const std::string& qual) {
        for (char q : qual) {
            if (q < '!' || q > '~') return false; // Phred+33 range
        }
        return true;
    }
    
    double getAverageQuality(const std::string& qual) {
        if (qual.empty()) return 0.0;
        
        double sum = 0.0;
        for (char q : qual) {
            sum += (q - 33); // Convert to Phred score
        }
        return sum / qual.length();
    }
    
    // Orientation analysis
    OrientationInfo analyzeOrientation(const std::string& seq) {
        OrientationInfo info;
        info.is_forward = true; // Default assumption
        
        // Check for polyG at start
        info.polyG_length = 0;
        info.N_count_in_polyG = 0;
        for (size_t i = 0; i < seq.length(); i++) {
            char base = std::toupper(seq[i]);
            if (base == 'G') {
                info.polyG_length++;
            } else if (base == 'N') {
                info.polyG_length++;
                info.N_count_in_polyG++;
            } else {
                break;
            }
        }
        info.has_polyG_start = (info.polyG_length >= 3);
        
        // Check for polyC at end
        info.polyC_length = 0;
        info.N_count_in_polyC = 0;
        for (int i = seq.length() - 1; i >= 0; i--) {
            char base = std::toupper(seq[i]);
            if (base == 'C') {
                info.polyC_length++;
            } else if (base == 'N') {
                info.polyC_length++;
                info.N_count_in_polyC++;
            } else {
                break;
            }
        }
        info.has_polyC_end = (info.polyC_length >= 3);
        
        return info;
    }
    
    OrientationInfo analyzeOrientationWithQuality(const std::string& seq, const std::string& qual) {
        OrientationInfo info = analyzeOrientation(seq);
        
        if (seq.length() == qual.length()) {
            // Calculate average quality in polyG region
            if (info.polyG_length > 0) {
                double sum = 0.0;
                for (int i = 0; i < info.polyG_length; i++) {
                    sum += (qual[i] - 33);
                }
                info.avg_quality_polyG = sum / info.polyG_length;
            }
            
            // Calculate average quality in polyC region
            if (info.polyC_length > 0) {
                double sum = 0.0;
                int start = seq.length() - info.polyC_length;
                for (int i = start; i < (int)seq.length(); i++) {
                    sum += (qual[i] - 33);
                }
                info.avg_quality_polyC = sum / info.polyC_length;
            }
        }
        
        return info;
    }
    
    // N base analysis
    NBaseAnalysis analyzeNBases(const std::string& seq, const std::string& qual) {
        NBaseAnalysis analysis;
        analysis.total_N_count = 0;
        analysis.N_in_first_20bp = 0;
        analysis.N_in_last_20bp = 0;
        analysis.avg_quality_at_N = 0.0;
        
        for (size_t i = 0; i < seq.length(); i++) {
            if (std::toupper(seq[i]) == 'N') {
                analysis.total_N_count++;
                analysis.N_positions.push_back(i);
                
                if (i < 20) analysis.N_in_first_20bp++;
                if (i >= seq.length() - 20) analysis.N_in_last_20bp++;
                
                if (!qual.empty() && i < qual.length()) {
                    analysis.N_qualities.push_back(qual[i]);
                }
            }
        }
        
        if (!analysis.N_qualities.empty()) {
            double sum = 0.0;
            for (char q : analysis.N_qualities) {
                sum += (q - 33);
            }
            analysis.avg_quality_at_N = sum / analysis.N_qualities.size();
        }
        
        return analysis;
    }
    
    // Quality-based base calling utilities
    bool isLowQualityBase(char quality, char threshold) {
        return quality < threshold;
    }
    
    bool isPotentialMiscall(char base, char quality, char threshold) {
        return quality < threshold;
    }
    
    char correctPotentialMiscall(char base, char expected_base, char quality, char quality_threshold) {
        if (quality < quality_threshold) {
            return expected_base;
        }
        return base;
    }
}
