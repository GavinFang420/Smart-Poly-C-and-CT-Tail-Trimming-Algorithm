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
        case 'N': return 'N';
        default: return 'N';
    }
}

// Core overlap detection algorithm (adapted from fastp)
OverlapResult ReadMerger::detectOverlap(const std::string& r1_seq, 
                                       const std::string& r2_seq,
                                       const std::string& r1_qual,
                                       const std::string& r2_qual) {
    OverlapResult result;
    
    if (r1_seq.empty() || r2_seq.empty()) {
        return result;
    }
    
    // Get reverse complement of read2 for overlap analysis
    std::string r2_rc = reverseComplement(r2_seq);
    std::string r2_qual_rev = r2_qual;
    if (!r2_qual_rev.empty()) {
        std::reverse(r2_qual_rev.begin(), r2_qual_rev.end());
    }
    
    int r1_len = r1_seq.length();
    int r2_len = r2_rc.length();
    
    int best_offset = -1;
    int best_overlap_len = 0;
    int best_diff_count = 999999;
    
    // Try different overlap positions
    // Case 1: R2 reverse complement overlaps with R1 tail
    for (int offset = 0; offset <= r1_len - min_overlap_len; offset++) {
        int max_overlap = std::min(r1_len - offset, r2_len);
        
        for (int overlap_len = min_overlap_len; overlap_len <= max_overlap; overlap_len++) {
            // Get overlapping regions
            std::string r1_region = r1_seq.substr(offset, overlap_len);
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
                best_offset = offset;
                best_overlap_len = overlap_len;
                best_diff_count = diff_count;
            }
        }
    }
    
    // Case 2: R1 overlaps with R2 reverse complement tail
    for (int offset = 0; offset <= r2_len - min_overlap_len; offset++) {
        int max_overlap = std::min(r2_len - offset, r1_len);
        
        for (int overlap_len = min_overlap_len; overlap_len <= max_overlap; overlap_len++) {
            std::string r1_region = r1_seq.substr(0, overlap_len);
            std::string r2_region = r2_rc.substr(offset, overlap_len);
            
            int diff_count = 0;
            for (int i = 0; i < overlap_len; i++) {
                if (std::toupper(r1_region[i]) != std::toupper(r2_region[i])) {
                    diff_count++;
                }
            }
            
            double diff_percent = (double)diff_count / overlap_len;
            
            bool valid_overlap = (diff_count <= max_diff_count) && 
                               (diff_percent <= max_diff_percent);
            
            if (valid_overlap && (diff_count < best_diff_count || 
                (diff_count == best_diff_count && overlap_len > best_overlap_len))) {
                best_offset = -offset - 1; // Negative to indicate R2 extends beyond R1
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

// Main merge function - with quality scores
MergeResult ReadMerger::mergeReads(const std::string& r1_seq, const std::string& r2_seq,
                                  const std::string& r1_qual, const std::string& r2_qual) {
    MergeResult result;
    
    OverlapResult overlap = detectOverlap(r1_seq, r2_seq, r1_qual, r2_qual);
    
    if (!overlap.overlapped) {
        return result; // No valid overlap found
    }
    
    // Prepare reverse complement of R2
    std::string r2_rc = reverseComplement(r2_seq);
    std::string r2_qual_rev = r2_qual;
    if (!r2_qual_rev.empty()) {
        std::reverse(r2_qual_rev.begin(), r2_qual_rev.end());
    }
    
    std::string merged_seq;
    std::string merged_qual;
    
    if (overlap.offset >= 0) {
        // Case 1: R1 prefix + overlapped region + R2 suffix
        // R1: ----[overlap]
        // R2:     [overlap]----
        
        // Add R1 prefix before overlap
        merged_seq = r1_seq.substr(0, overlap.offset);
        if (!r1_qual.empty()) {
            merged_qual = r1_qual.substr(0, overlap.offset);
        }
        
        // Add overlapped region (prefer R1 bases, merge qualities)
        for (int i = 0; i < overlap.overlap_len; i++) {
            char r1_base = r1_seq[overlap.offset + i];
            char r2_base = r2_rc[i];
            
            // Prefer R1 base (fastp behavior), but could be improved with quality consideration
            merged_seq += r1_base;
            
            if (!r1_qual.empty() && !r2_qual_rev.empty()) {
                char r1_q = r1_qual[overlap.offset + i];
                char r2_q = r2_qual_rev[i];
                merged_qual += MergeUtils::mergeQualityScores(r1_q, r2_q);
            } else if (!r1_qual.empty()) {
                merged_qual += r1_qual[overlap.offset + i];
            }
        }
        
        // Add R2 suffix after overlap
        if (overlap.overlap_len < r2_rc.length()) {
            merged_seq += r2_rc.substr(overlap.overlap_len);
            if (!r2_qual_rev.empty()) {
                merged_qual += r2_qual_rev.substr(overlap.overlap_len);
            }
        }
        
        result.r1_bases = overlap.offset + overlap.overlap_len;
        result.r2_bases = r2_rc.length() - overlap.overlap_len;
        
    } else {
        // Case 2: R2 prefix + overlapped region + R1 suffix  
        // R2: ----[overlap]
        // R1:     [overlap]----
        
        int r2_offset = -(overlap.offset + 1);
        
        // Add R2 prefix before overlap
        merged_seq = r2_rc.substr(0, r2_offset);
        if (!r2_qual_rev.empty()) {
            merged_qual = r2_qual_rev.substr(0, r2_offset);
        }
        
        // Add overlapped region
        for (int i = 0; i < overlap.overlap_len; i++) {
            char r1_base = r1_seq[i];
            char r2_base = r2_rc[r2_offset + i];
            
            merged_seq += r1_base; // Prefer R1
            
            if (!r1_qual.empty() && !r2_qual_rev.empty()) {
                char r1_q = r1_qual[i];
                char r2_q = r2_qual_rev[r2_offset + i];
                merged_qual += MergeUtils::mergeQualityScores(r1_q, r2_q);
            } else if (!r1_qual.empty()) {
                merged_qual += r1_qual[i];
            }
        }
        
        // Add R1 suffix after overlap
        if (overlap.overlap_len < r1_seq.length()) {
            merged_seq += r1_seq.substr(overlap.overlap_len);
            if (!r1_qual.empty()) {
                merged_qual += r1_qual.substr(overlap.overlap_len);
            }
        }
        
        result.r1_bases = r1_seq.length() - overlap.overlap_len;
        result.r2_bases = r2_offset + overlap.overlap_len;
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
    if (start < 0 || start >= seq.length()) {
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
        
        for (char base : seq) {
            char upper_base = std::toupper(base);
            if (upper_base == 'G' || upper_base == 'C') {
                stats.gc_count++;
            } else if (upper_base == 'N') {
                stats.n_count++;
            }
        }
        
        stats.gc_percent = stats.length > 0 ? 
            (double)stats.gc_count / stats.length * 100.0 : 0.0;
        
        return stats;
    }
}
