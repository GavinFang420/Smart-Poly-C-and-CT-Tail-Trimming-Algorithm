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

// Modified overlap detection for PolyC tail processing with relaxed requirements
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

// New function: concatenate reads when no overlap found
MergeResult ReadMerger::concatenateReads(const std::string& r1_seq, 
                                        const std::string& r2_seq,
                                        const std::string& r1_qual,
                                        const std::string& r2_qual) {
    MergeResult result;
    
    // Get reverse complement of R2 for consistent orientation
    std::string r2_rc = reverseComplement(r2_seq);
    std::string r2_qual_rev = r2_qual;
    if (!r2_qual_rev.empty()) {
        std::reverse(r2_qual_rev.begin(), r2_qual_rev.end());
    }
    
    // Simple concatenation: R1 + R2_RC
    result.sequence = r1_seq + r2_rc;
    
    // Concatenate quality scores if available
    if (!r1_qual.empty() && !r2_qual_rev.empty()) {
        result.quality = r1_qual + r2_qual_rev;
    } else if (!r1_qual.empty()) {
        result.quality = r1_qual;
    }
    
    // Set result fields
    result.merged = true;        // Successfully processed
    result.overlapped = false;   // No overlap found, just concatenated
    result.r1_bases = r1_seq.length();
    result.r2_bases = r2_rc.length();
    result.r1_start_in_merged = 0;
    result.r2_end_in_merged = r1_seq.length() + r2_rc.length();
    result.overlap_length = 0;   // No overlap
    result.merge_offset = r1_seq.length();  // R2 starts right after R1
    
    return result;
}

// Main merge function - sequence only
MergeResult ReadMerger::mergeReads(const std::string& r1_seq, const std::string& r2_seq) {
    return mergeReads(r1_seq, r2_seq, "", "");
}

// Modified merge function for PolyC tail processing with concatenation fallback
// Result structure: R2_head(5') + overlap + R1_head(3') OR R1 + R2_RC (concatenated)
MergeResult ReadMerger::mergeReads(const std::string& r1_seq, const std::string& r2_seq,
                                  const std::string& r1_qual, const std::string& r2_qual) {
    MergeResult result;
    
    // First try to find overlap
    OverlapResult overlap = detectOverlap(r1_seq, r2_seq, r1_qual, r2_qual);
    
    if (!overlap.overlapped) {
        // No overlap found - check if concatenation is enabled
        if (enable_concatenation) {
            return concatenateReads(r1_seq, r2_seq, r1_qual, r2_qual);
        } else {
            return result; // Return empty result if concatenation disabled
        }
    }
    
    // Overlap found - proceed with merging
    result.overlapped = true;
    
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
        if (overlap.offset > 0) {
            merged_seq = r1_seq.substr(0, overlap.offset);
            if (!r1_qual.empty()) {
                merged_qual = r1_qual.substr(0, overlap.offset);
            }
        }
        
        for (int i = 0; i < overlap.overlap_len; i++) {
            char r1_base = r1_seq[overlap.offset + i];
            char r2_base = r2_rc[i];
            
            if (!r1_qual.empty() && !r2_qual_rev.empty()) {
                char r1_q = r1_qual[overlap.offset + i];
                char r2_q = r2_qual_rev[i];
                
                if (std::toupper(r1_base) == 'N' && std::toupper(r2_base) != 'N') {
                    merged_seq += r2_base;
                    merged_qual += r2_q;
                } else if (std::toupper(r2_base) == 'N' && std::toupper(r1_base) != 'N') {
                    merged_seq += r1_base;
                    merged_qual += r1_q;
                } else if (std::toupper(r1_base) == 'N' && std::toupper(r2_base) == 'N') {
                    if (i >= overlap.overlap_len - 5) {
                        merged_seq += 'C'; // Assume N is miscalled C in polyC region
                        merged_qual += std::max(r1_q, r2_q);
                    } else {
                        merged_seq += 'N';
                        merged_qual += std::max(r1_q, r2_q);
                    }
                } else {
                    if (r1_q >= r2_q) {
                        merged_seq += r1_base;
                        merged_qual += r1_q;
                    } else {
                        merged_seq += r2_base;
                        merged_qual += r2_q;
                    }
                }
            } else {
                if (std::toupper(r1_base) == 'N' && std::toupper(r2_base) != 'N') {
                    merged_seq += r2_base;
                } else if (std::toupper(r2_base) == 'N' && std::toupper(r1_base) != 'N') {
                    merged_seq += r1_base;
                } else {
                    merged_seq += r1_base;
                }
                if (!r1_qual.empty()) {
                    merged_qual += r1_qual[i];
                }
            }
        }
        
        if (overlap.overlap_len < (int)r1_seq.length()) {
            merged_seq += r1_seq.substr(overlap.overlap_len);
            if (!r1_qual.empty()) {
                merged_qual += r1_qual.substr(overlap.overlap_len);
            }
        }
        
        int r2_remaining_start = overlap.offset + overlap.overlap_len;
        if (r2_remaining_start < (int)r2_rc.length()) {
            merged_seq += r2_rc.substr(r2_remaining_start);
            if (!r2_qual_rev.empty()) {
                merged_qual += r2_qual_rev.substr(r2_remaining_start);
            }
        }
        
        result.r1_bases = r1_seq.length();
        result.r2_bases = r2_rc.length() - overlap.overlap_len;
        
        result.r1_start_in_merged = 0;
        result.r2_end_in_merged = overlap.offset + overlap.overlap_len + (r2_rc.length() - overlap.overlap_len);
        result.overlap_length = overlap.overlap_len;
        result.merge_offset = overlap.offset;
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
