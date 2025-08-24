#ifndef MERGEREAD_H
#define MERGEREAD_H

#include <string>
#include <vector>
#include <algorithm>

struct OverlapResult {
    bool overlapped;        // Whether overlap was found
    int offset;            // Offset position where overlap starts (positive: in R1, negative: in R2_RC)
    int overlap_len;       // Length of the overlapped region
    int diff_count;        // Number of mismatches in overlap
    double diff_percent;   // Percentage of mismatches
    
    OverlapResult() : overlapped(false), offset(0), overlap_len(0), 
                     diff_count(0), diff_percent(0.0) {}
};

struct MergeResult {
    bool merged;           // Whether reads were successfully processed (merged or concatenated)
    bool overlapped;       // Whether actual overlap was found (true=merged, false=concatenated)
    std::string sequence;  // Final sequence (merged or concatenated)
    std::string quality;   // Final quality string (if provided)
    int r1_bases;         // Number of bases contributed from read1
    int r2_bases;         // Number of bases contributed from read2 (after reverse complement)
    
    // 对齐信息 - 用于正确的trim位置映射
    int r1_start_in_merged;  // R1在merged序列中的起始位置
    int r2_end_in_merged;    // R2在merged序列中的结束位置  
    int overlap_length;      // overlap的长度 (0 if concatenated)
    int merge_offset;        // merge偏移：R1头到R2原本尾巴新位置的距离（可正可负）
    
    MergeResult() : merged(false), overlapped(false), r1_bases(0), r2_bases(0), 
                   r1_start_in_merged(-1), r2_end_in_merged(-1), 
                   overlap_length(0), merge_offset(0) {}
};

class ReadMerger {
private:
    // Default parameters (can be adjusted)
    int min_overlap_len;      // Minimum overlap length required (reduced to 5)
    int max_diff_count;       // Maximum allowed mismatches
    double max_diff_percent;  // Maximum allowed mismatch percentage
    bool enable_concatenation; // Whether to concatenate when no overlap found
    
    // Internal helper functions
    std::string reverseComplement(const std::string& seq);
    bool isValidBase(char base);
    char getComplementBase(char base);
    
    // Core overlap detection algorithm (modified for PolyC tail processing)
    // For PolyC tail: R2 starts with GG tail (complement of polyC tail)
    // Logic: complement R2 first, then find overlap with R1
    // Result: R1 head at 3' end, R2 head at 5' end
    OverlapResult detectOverlap(const std::string& r1_seq, 
                               const std::string& r2_seq,
                               const std::string& r1_qual = "",
                               const std::string& r2_qual = "");
    
    // New function: concatenate reads when no overlap found
    MergeResult concatenateReads(const std::string& r1_seq, 
                                const std::string& r2_seq,
                                const std::string& r1_qual = "",
                                const std::string& r2_qual = "");
    
    // Score overlap region considering quality scores
    double scoreOverlapRegion(const std::string& r1_region,
                             const std::string& r2_region,
                             const std::string& r1_qual_region = "",
                             const std::string& r2_qual_region = "");

public:
    // Constructor with default parameters optimized for PolyC tail processing
    // min_overlap: minimum overlap length (reduced from 15 to 5)
    // max_diff: maximum allowed mismatches 
    // max_diff_pct: maximum allowed mismatch percentage
    // enable_concat: whether to concatenate when no overlap found
    ReadMerger(int min_overlap = 5, int max_diff = 2, double max_diff_pct = 0.20, 
               bool enable_concat = true) 
        : min_overlap_len(min_overlap), max_diff_count(max_diff), 
          max_diff_percent(max_diff_pct), enable_concatenation(enable_concat) {}
    
    // Main merge function - sequence only
    // Processes PolyC tail data: R2 starts with GG tail
    // Returns merged sequence with R2_head(5') + overlap + R1_head(3')
    // If no overlap found and concatenation enabled, returns concatenated sequence
    MergeResult mergeReads(const std::string& r1_seq, const std::string& r2_seq);
    
    // Main merge function - with quality scores
    // Merges quality scores by selecting higher quality base in overlap regions
    // If no overlap found and concatenation enabled, concatenates with quality scores
    MergeResult mergeReads(const std::string& r1_seq, const std::string& r2_seq,
                          const std::string& r1_qual, const std::string& r2_qual);
    
    // Get overlap analysis without merging
    // Useful for debugging and quality control
    OverlapResult analyzeOverlap(const std::string& r1_seq, const std::string& r2_seq);
    
    // Parameter setters
    void setMinOverlapLength(int len) { min_overlap_len = len; }
    void setMaxDiffCount(int count) { max_diff_count = count; }
    void setMaxDiffPercent(double percent) { max_diff_percent = percent; }
    void setConcatenationEnabled(bool enabled) { enable_concatenation = enabled; }
    
    // Parameter getters
    int getMinOverlapLength() const { return min_overlap_len; }
    int getMaxDiffCount() const { return max_diff_count; }
    double getMaxDiffPercent() const { return max_diff_percent; }
    bool isConcatenationEnabled() const { return enable_concatenation; }
    
    // Utility functions
    static std::string getSubsequence(const std::string& seq, int start, int length);
    static bool hasValidOverlap(const OverlapResult& result, int min_len, 
                               int max_diff, double max_pct);
    
    // PolyC-specific utility functions
    static bool hasPolyGStart(const std::string& seq, int min_length = 3, 
                             bool count_N_as_G = true);
    static int countPolyGLength(const std::string& seq, bool count_N_as_G = true);
    static std::string getOriginalR2Orientation(const std::string& merged_seq, 
                                               const MergeResult& result);
    
    // Quality-aware PolyC analysis
    static bool hasPolyGStartWithQuality(const std::string& seq, const std::string& qual,
                                        int min_length = 3, char min_quality = '!',
                                        bool count_N_as_G = true);
    static int countPolyGLengthWithQuality(const std::string& seq, const std::string& qual,
                                          char min_quality = '!', bool count_N_as_G = true);
};

// Standalone utility functions
namespace MergeUtils {
    // Convert quality score to probability
    double qualToProb(char qual);
    
    // Convert probability to quality score
    char probToQual(double prob);
    
    // Merge two quality scores (takes higher quality)
    char mergeQualityScores(char q1, char q2);
    
    // Check if sequence contains only valid DNA bases
    bool isValidDNASequence(const std::string& seq);
    
    // Simple sequence statistics
    struct SequenceStats {
        int length;
        int gc_count;
        int n_count;
        double gc_percent;
        int polyG_start_length;  // Length of polyG at start (for R2 analysis)
        int polyC_end_length;    // Length of polyC at end (for merged sequence analysis)
    };
    
    SequenceStats getSequenceStats(const std::string& seq);
    
    // PolyC tail specific utilities with quality awareness
    bool hasPolyTail(const std::string& seq, char base, int min_length = 3, 
                    bool count_N_as_base = true);
    int getPolyTailLength(const std::string& seq, char base, bool count_N_as_base = true);
    std::string trimPolyTail(const std::string& seq, char base, int min_length = 3,
                           bool count_N_as_base = true);
    
    // Quality-aware poly tail analysis
    bool hasPolyTailWithQuality(const std::string& seq, const std::string& qual,
                               char base, int min_length = 3, char min_quality = '!',
                               bool count_N_as_base = true);
    int getPolyTailLengthWithQuality(const std::string& seq, const std::string& qual,
                                    char base, char min_quality = '!', 
                                    bool count_N_as_base = true);
    std::string trimPolyTailWithQuality(const std::string& seq, const std::string& qual,
                                       char base, int min_length = 3, char min_quality = '!',
                                       bool count_N_as_base = true);
    
    // Reverse complement utility (exposed for external use)
    std::string reverseComplement(const std::string& seq);
    
    // Quality validation
    bool isValidQualityString(const std::string& qual);
    double getAverageQuality(const std::string& qual);
    
    // Sequence orientation utilities for PolyC processing
    struct OrientationInfo {
        bool is_forward;           // True if sequence is in forward orientation
        bool has_polyG_start;      // True if starts with polyG (R2 characteristic)
        bool has_polyC_end;        // True if ends with polyC (merged characteristic)
        int polyG_length;          // Length of polyG at start
        int polyC_length;          // Length of polyC at end
        int N_count_in_polyG;      // Number of N bases in polyG region (potential C bases)
        int N_count_in_polyC;      // Number of N bases in polyC region (potential C bases)
        double avg_quality_polyG;   // Average quality in polyG region
        double avg_quality_polyC;   // Average quality in polyC region
    };
    
    OrientationInfo analyzeOrientation(const std::string& seq);
    OrientationInfo analyzeOrientationWithQuality(const std::string& seq, 
                                                  const std::string& qual);
    
    // N base analysis utilities
    struct NBaseAnalysis {
        int total_N_count;                    // Total N bases in sequence
        int N_in_first_20bp;                  // N bases in first 20bp (R2 head region)
        int N_in_last_20bp;                   // N bases in last 20bp (potential polyC region)
        std::vector<int> N_positions;         // Positions of all N bases
        std::vector<char> N_qualities;        // Quality scores at N positions
        double avg_quality_at_N;              // Average quality at N positions
    };
    
    NBaseAnalysis analyzeNBases(const std::string& seq, const std::string& qual = "");
    
    // Quality-based base calling utilities
    bool isLowQualityBase(char quality, char threshold = '#');  // Phred score < 3
    bool isPotentialMiscall(char base, char quality, char threshold = '('); // Phred score < 8
    char correctPotentialMiscall(char base, char expected_base, char quality, 
                               char quality_threshold = '(');
}

#endif