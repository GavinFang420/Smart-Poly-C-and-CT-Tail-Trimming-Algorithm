#ifndef MERGEREAD_H
#define MERGEREAD_H

#include <string>
#include <algorithm>

struct OverlapResult {
    bool overlapped;        // Whether overlap was found
    int offset;            // Offset position in read1 where overlap starts
    int overlap_len;       // Length of the overlapped region
    int diff_count;        // Number of mismatches in overlap
    double diff_percent;   // Percentage of mismatches
    
    OverlapResult() : overlapped(false), offset(0), overlap_len(0), 
                     diff_count(0), diff_percent(0.0) {}
};

struct MergeResult {
    bool merged;           // Whether reads were successfully merged
    std::string sequence;  // Merged sequence
    std::string quality;   // Merged quality string (if provided)
    int r1_bases;         // Number of bases from read1
    int r2_bases;         // Number of bases from read2
    
    MergeResult() : merged(false), r1_bases(0), r2_bases(0) {}
};

class ReadMerger {
private:
    // Default parameters (can be adjusted)
    int min_overlap_len;      // Minimum overlap length required
    int max_diff_count;       // Maximum allowed mismatches
    double max_diff_percent;  // Maximum allowed mismatch percentage
    
    // Internal helper functions
    std::string reverseComplement(const std::string& seq);
    bool isValidBase(char base);
    char getComplementBase(char base);
    
    // Core overlap detection algorithm (extracted from fastp)
    OverlapResult detectOverlap(const std::string& r1_seq, 
                               const std::string& r2_seq,
                               const std::string& r1_qual = "",
                               const std::string& r2_qual = "");
    
    // Score overlap region considering quality scores
    double scoreOverlapRegion(const std::string& r1_region,
                             const std::string& r2_region,
                             const std::string& r1_qual_region = "",
                             const std::string& r2_qual_region = "");

public:
    // Constructor with default fastp parameters
    ReadMerger(int min_overlap = 30, int max_diff = 5, double max_diff_pct = 0.2) 
        : min_overlap_len(min_overlap), max_diff_count(max_diff), 
          max_diff_percent(max_diff_pct) {}
    
    // Main merge function - sequence only
    MergeResult mergeReads(const std::string& r1_seq, const std::string& r2_seq);
    
    // Main merge function - with quality scores
    MergeResult mergeReads(const std::string& r1_seq, const std::string& r2_seq,
                          const std::string& r1_qual, const std::string& r2_qual);
    
    // Get overlap analysis without merging
    OverlapResult analyzeOverlap(const std::string& r1_seq, const std::string& r2_seq);
    
    // Parameter setters
    void setMinOverlapLength(int len) { min_overlap_len = len; }
    void setMaxDiffCount(int count) { max_diff_count = count; }
    void setMaxDiffPercent(double percent) { max_diff_percent = percent; }
    
    // Parameter getters
    int getMinOverlapLength() const { return min_overlap_len; }
    int getMaxDiffCount() const { return max_diff_count; }
    double getMaxDiffPercent() const { return max_diff_percent; }
    
    // Utility functions
    static std::string getSubsequence(const std::string& seq, int start, int length);
    static bool hasValidOverlap(const OverlapResult& result, int min_len, 
                               int max_diff, double max_pct);
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
    };
    
    SequenceStats getSequenceStats(const std::string& seq);
}

#endif
