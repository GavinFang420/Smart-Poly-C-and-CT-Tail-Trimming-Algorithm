#include "mergeread.h"
#include "smarttrim.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include <functional>

// FASTQ record structure
struct FastqRecord {
    std::string header;
    std::string sequence;
    std::string plus_line;
    std::string quality;
    
    bool is_valid() const {
        return !header.empty() && !sequence.empty() && 
               !quality.empty() && sequence.length() == quality.length();
    }
};

// Processing statistics
struct ProcessingStats {
    int total_reads = 0;
    int merged_reads = 0;
    int trimmed_reads = 0;
    int r1_trimmed = 0;
    int r2_trimmed = 0;
    int both_trimmed = 0;
    
    // Quality metrics
    double avg_r1_length = 0.0;
    double avg_r2_length = 0.0;
    double avg_merged_length = 0.0;
    double avg_r1_trim_length = 0.0;
    double avg_r2_trim_length = 0.0;
    
    // Tail analysis - extended for CT detection
    int r1_with_polyc = 0;
    int r1_with_polyt = 0;
    int r1_with_ct_tail = 0;
    int r2_with_polyg = 0;
    int r2_with_polya = 0;
    int r2_with_ag_tail = 0;
    
    double avg_polyc_length = 0.0;
    double avg_polyt_length = 0.0;
    double avg_ct_tail_length = 0.0;
    double avg_polyg_length = 0.0;
    double avg_polya_length = 0.0;
    double avg_ag_tail_length = 0.0;
    
    // Score distribution
    std::vector<double> trim_scores;
    std::map<std::string, int> score_ranges;
    
    // CT specific metrics
    int ct_pattern_detected = 0;
    int c_only_detected = 0;
    int mixed_tail_detected = 0;
    
    void updateStats(const FastqRecord& r1, const FastqRecord& r2, 
                    const MergeResult& merge_result, const TrimResult& trim_result) {
        total_reads++;
        
        // Basic length stats
        avg_r1_length = (avg_r1_length * (total_reads - 1) + r1.sequence.length()) / total_reads;
        avg_r2_length = (avg_r2_length * (total_reads - 1) + r2.sequence.length()) / total_reads;
        
        if (merge_result.merged) {
            merged_reads++;
            avg_merged_length = (avg_merged_length * (merged_reads - 1) + 
                               merge_result.sequence.length()) / merged_reads;
        }
        
        if (trim_result.is_valid) {
            trimmed_reads++;
            trim_scores.push_back(trim_result.final_score);
            
            if (trim_result.r1_trim_pos > 0) {
                r1_trimmed++;
                avg_r1_trim_length = (avg_r1_trim_length * (r1_trimmed - 1) + 
                                    trim_result.r1_trim_pos) / r1_trimmed;
            }
            
            if (trim_result.r2_trim_pos > 0) {
                r2_trimmed++;
                avg_r2_trim_length = (avg_r2_trim_length * (r2_trimmed - 1) + 
                                    trim_result.r2_trim_pos) / r2_trimmed;
            }
            
            if (trim_result.r1_trim_pos > 0 && trim_result.r2_trim_pos > 0) {
                both_trimmed++;
            }
            
            // Score range classification
            int score_int = (int)trim_result.final_score;
            if (score_int < 0) score_ranges["negative"]++;
            else if (score_int < 10) score_ranges["0-10"]++;
            else if (score_int < 20) score_ranges["10-20"]++;
            else if (score_int < 50) score_ranges["20-50"]++;
            else score_ranges["50+"]++;
        }
        
        // Enhanced tail content analysis
        analyzeTailContent(r1.sequence, r2.sequence, trim_result);
    }
    
private:
    void analyzeTailContent(const std::string& r1_seq, const std::string& r2_seq, 
                           const TrimResult& trim_result) {
        // Analyze R1 tail (last 20bp or trimmed region)
        int analyze_length = 20;
        if (trim_result.is_valid && trim_result.r1_trim_pos > 0) {
            analyze_length = std::min(trim_result.r1_trim_pos, 20);
        }
        
        if (analyze_length > 0 && r1_seq.length() >= analyze_length) {
            std::string r1_tail = r1_seq.substr(r1_seq.length() - analyze_length);
            
            // Count different bases
            int c_count = 0, t_count = 0, total_ct = 0;
            for (char base : r1_tail) {
                char upper_base = std::toupper(base);
                if (upper_base == 'C') { c_count++; total_ct++; }
                else if (upper_base == 'T') { t_count++; total_ct++; }
                else if (upper_base == 'N') { total_ct++; } // N could be C or T
            }
            
            double ct_ratio = (double)total_ct / analyze_length;
            double c_ratio = (double)c_count / analyze_length;
            double t_ratio = (double)t_count / analyze_length;
            
            // Classify tail type
            if (ct_ratio >= 0.7) {
                if (c_count >= 3 && t_count >= 2) {
                    // Mixed CT tail
                    r1_with_ct_tail++;
                    ct_pattern_detected++;
                    avg_ct_tail_length = (avg_ct_tail_length * (r1_with_ct_tail - 1) + 
                                        total_ct) / r1_with_ct_tail;
                } else if (c_ratio >= 0.7) {
                    // Mostly C
                    r1_with_polyc++;
                    c_only_detected++;
                    avg_polyc_length = (avg_polyc_length * (r1_with_polyc - 1) + 
                                      c_count) / r1_with_polyc;
                } else if (t_ratio >= 0.7) {
                    // Mostly T
                    r1_with_polyt++;
                    avg_polyt_length = (avg_polyt_length * (r1_with_polyt - 1) + 
                                      t_count) / r1_with_polyt;
                }
            } else if (ct_ratio >= 0.4) {
                // Mixed tail with some CT content
                mixed_tail_detected++;
            }
        }
        
        // Analyze R2 head (first 20bp or trimmed region)
        analyze_length = 20;
        if (trim_result.is_valid && trim_result.r2_trim_pos > 0) {
            analyze_length = std::min(trim_result.r2_trim_pos, 20);
        }
        
        if (analyze_length > 0 && r2_seq.length() >= analyze_length) {
            std::string r2_head = r2_seq.substr(0, analyze_length);
            
            int g_count = 0, a_count = 0, total_ag = 0;
            for (char base : r2_head) {
                char upper_base = std::toupper(base);
                if (upper_base == 'G') { g_count++; total_ag++; }
                else if (upper_base == 'A') { a_count++; total_ag++; }
                else if (upper_base == 'N') { total_ag++; }
            }
            
            double ag_ratio = (double)total_ag / analyze_length;
            double g_ratio = (double)g_count / analyze_length;
            double a_ratio = (double)a_count / analyze_length;
            
            if (ag_ratio >= 0.7) {
                if (g_count >= 3 && a_count >= 2) {
                    r2_with_ag_tail++;
                    avg_ag_tail_length = (avg_ag_tail_length * (r2_with_ag_tail - 1) + 
                                        total_ag) / r2_with_ag_tail;
                } else if (g_ratio >= 0.7) {
                    r2_with_polyg++;
                    avg_polyg_length = (avg_polyg_length * (r2_with_polyg - 1) + 
                                      g_count) / r2_with_polyg;
                } else if (a_ratio >= 0.7) {
                    r2_with_polya++;
                    avg_polya_length = (avg_polya_length * (r2_with_polya - 1) + 
                                      a_count) / r2_with_polya;
                }
            }
        }
    }
};

// Dataset configuration
struct DatasetConfig {
    std::string name;
    std::string description;
    std::string r1_file;
    std::string r2_file;
    TrimParams trim_params;
    std::string output_prefix;
};

class ComprehensiveProcessor {
private:
    bool readFastqRecord(std::ifstream& file, FastqRecord& record) {
        std::string line;
        
        if (!std::getline(file, record.header) || record.header.empty() || record.header[0] != '@') {
            return false;
        }
        
        if (!std::getline(file, record.sequence)) {
            return false;
        }
        
        if (!std::getline(file, record.plus_line) || record.plus_line.empty() || record.plus_line[0] != '+') {
            return false;
        }
        
        if (!std::getline(file, record.quality)) {
            return false;
        }
        
        return record.is_valid();
    }
    
    void writeFastqRecord(std::ofstream& file, const FastqRecord& record) {
        file << record.header << "\n"
             << record.sequence << "\n"
             << record.plus_line << "\n"
             << record.quality << "\n";
    }

public:
    ProcessingStats processDataset(const DatasetConfig& config, int max_reads = 0, bool verbose = false) {
        std::ifstream r1_file(config.r1_file);
        std::ifstream r2_file(config.r2_file);
        
        if (!r1_file.is_open()) {
            std::cerr << "Error: Cannot open " << config.r1_file << std::endl;
            return ProcessingStats();
        }
        
        if (!r2_file.is_open()) {
            std::cerr << "Error: Cannot open " << config.r2_file << std::endl;
            return ProcessingStats();
        }
        
        // Output files
        std::ofstream r1_out(config.output_prefix + "_R1.fastq");
        std::ofstream r2_out(config.output_prefix + "_R2.fastq");
        std::ofstream merged_out(config.output_prefix + "_merged.fastq");
        std::ofstream stats_out(config.output_prefix + "_stats.txt");
        
        SmartTrimmer trimmer(config.trim_params);
        ReadMerger merger(10, 3, 0.3);
        ProcessingStats stats;
        
        FastqRecord r1_record, r2_record;
        int processed = 0;
        
        std::cout << "\nProcessing " << config.description << "..." << std::endl;
        std::cout << "Files: " << config.r1_file << " + " << config.r2_file << std::endl;
        std::cout << "Config: " << config.trim_params.target_bases 
                  << " (C:" << config.trim_params.c_score 
                  << ", T:" << config.trim_params.t_score << ")" << std::endl;
        
        while (readFastqRecord(r1_file, r1_record) && readFastqRecord(r2_file, r2_record)) {
            processed++;
            
            if (max_reads > 0 && processed > max_reads) {
                break;
            }
            
            if (verbose && processed % 1000 == 0) {
                std::cout << "Processed " << processed << " read pairs..." << std::endl;
            }
            
            // Try to merge reads
            MergeResult merge_result = merger.mergeReads(r1_record.sequence, r2_record.sequence,
                                                        r1_record.quality, r2_record.quality);
            
            // Apply trimming
            TrimResult trim_result = trimmer.findOptimalTrimPositions(r1_record.sequence, r2_record.sequence);
            
            // Update statistics
            stats.updateStats(r1_record, r2_record, merge_result, trim_result);
            
            // Apply trimming to sequences
            FastqRecord trimmed_r1 = r1_record;
            FastqRecord trimmed_r2 = r2_record;
            
            if (trim_result.is_valid) {
                auto trimmed_seqs = trimmer.trimReads(r1_record.sequence, r2_record.sequence, trim_result);
                
                if (trim_result.r1_trim_pos > 0) {
                    trimmed_r1.sequence = trimmed_seqs.first;
                    trimmed_r1.quality = r1_record.quality.substr(0, trimmed_seqs.first.length());
                }
                
                if (trim_result.r2_trim_pos > 0) {
                    trimmed_r2.sequence = trimmed_seqs.second;
                    trimmed_r2.quality = r2_record.quality.substr(trim_result.r2_trim_pos);
                }
            }
            
            // Write output
            writeFastqRecord(r1_out, trimmed_r1);
            writeFastqRecord(r2_out, trimmed_r2);
            
            if (merge_result.merged) {
                FastqRecord merged_record;
                merged_record.header = r1_record.header;
                merged_record.sequence = merge_result.sequence;
                merged_record.plus_line = "+";
                merged_record.quality = merge_result.quality.empty() ? 
                    std::string(merge_result.sequence.length(), 'I') : merge_result.quality;
                writeFastqRecord(merged_out, merged_record);
            }
            
            // Detailed logging for first few reads
            if (verbose && processed <= 3) {
                std::cout << "\n--- Read " << processed << " (" << config.name << ") ---" << std::endl;
                std::cout << "R1: " << r1_record.sequence << std::endl;
                std::cout << "R2: " << r2_record.sequence << std::endl;
                
                if (merge_result.merged) {
                    std::cout << "Merged: " << merge_result.sequence << std::endl;
                }
                
                if (trim_result.is_valid) {
                    std::cout << "Trim: R1[" << trim_result.r1_trim_pos << "], R2[" 
                              << trim_result.r2_trim_pos << "] (score: " << trim_result.score_detail << ")" << std::endl;
                    std::cout << "After trim - R1: " << trimmed_r1.sequence << std::endl;
                    std::cout << "After trim - R2: " << trimmed_r2.sequence << std::endl;
                }
            }
        }
        
        r1_file.close();
        r2_file.close();
        r1_out.close();
        r2_out.close();
        merged_out.close();
        
        // Write statistics
        writeStatistics(stats_out, stats, config);
        stats_out.close();
        
        std::cout << "Completed: " << processed << " read pairs processed." << std::endl;
        return stats;
    }
    
    void writeStatistics(std::ofstream& file, const ProcessingStats& stats, const DatasetConfig& config) {
        file << "CT Tail Trimming Statistics - " << config.description << "\n";
        file << "=================================================\n\n";
        
        file << "Dataset:\n";
        file << "  Name: " << config.name << "\n";
        file << "  R1 file: " << config.r1_file << "\n";
        file << "  R2 file: " << config.r2_file << "\n\n";
        
        file << "Parameters:\n";
        file << "  Target bases: " << config.trim_params.target_bases << "\n";
        file << "  Window size: " << config.trim_params.window_size << "\n";
        file << "  C score: " << config.trim_params.c_score << "\n";
        file << "  T score: " << config.trim_params.t_score << "\n";
        file << "  A score: " << config.trim_params.a_score << "\n";
        file << "  G score: " << config.trim_params.g_score << "\n";
        file << "  N score: " << config.trim_params.n_score << "\n";
        file << "  Initial score: " << config.trim_params.initial_score << "\n\n";
        
        file << "Processing Summary:\n";
        file << "  Total read pairs: " << stats.total_reads << "\n";
        file << "  Merged reads: " << stats.merged_reads 
             << " (" << (double)stats.merged_reads / stats.total_reads * 100 << "%)\n";
        file << "  Trimmed reads: " << stats.trimmed_reads 
             << " (" << (double)stats.trimmed_reads / stats.total_reads * 100 << "%)\n";
        file << "  R1 trimmed: " << stats.r1_trimmed 
             << " (" << (double)stats.r1_trimmed / stats.total_reads * 100 << "%)\n";
        file << "  R2 trimmed: " << stats.r2_trimmed 
             << " (" << (double)stats.r2_trimmed / stats.total_reads * 100 << "%)\n";
        file << "  Both trimmed: " << stats.both_trimmed 
             << " (" << (double)stats.both_trimmed / stats.total_reads * 100 << "%)\n\n";
        
        file << "Length Statistics:\n";
        file << "  Avg R1 length: " << std::fixed << std::setprecision(1) << stats.avg_r1_length << "\n";
        file << "  Avg R2 length: " << stats.avg_r2_length << "\n";
        file << "  Avg merged length: " << stats.avg_merged_length << "\n";
        file << "  Avg R1 trim length: " << stats.avg_r1_trim_length << "\n";
        file << "  Avg R2 trim length: " << stats.avg_r2_trim_length << "\n\n";
        
        file << "Tail Analysis:\n";
        file << "  R1 with polyC: " << stats.r1_with_polyc 
             << " (" << (double)stats.r1_with_polyc / stats.total_reads * 100 << "%)\n";
        file << "  R1 with polyT: " << stats.r1_with_polyt 
             << " (" << (double)stats.r1_with_polyt / stats.total_reads * 100 << "%)\n";
        file << "  R1 with CT tail: " << stats.r1_with_ct_tail 
             << " (" << (double)stats.r1_with_ct_tail / stats.total_reads * 100 << "%)\n";
        file << "  R2 with polyG: " << stats.r2_with_polyg 
             << " (" << (double)stats.r2_with_polyg / stats.total_reads * 100 << "%)\n";
        file << "  R2 with polyA: " << stats.r2_with_polya 
             << " (" << (double)stats.r2_with_polya / stats.total_reads * 100 << "%)\n";
        file << "  R2 with AG head: " << stats.r2_with_ag_tail 
             << " (" << (double)stats.r2_with_ag_tail / stats.total_reads * 100 << "%)\n\n";
        
        file << "Pattern Detection:\n";
        file << "  CT pattern detected: " << stats.ct_pattern_detected 
             << " (" << (double)stats.ct_pattern_detected / stats.total_reads * 100 << "%)\n";
        file << "  C only detected: " << stats.c_only_detected 
             << " (" << (double)stats.c_only_detected / stats.total_reads * 100 << "%)\n";
        file << "  Mixed tail detected: " << stats.mixed_tail_detected 
             << " (" << (double)stats.mixed_tail_detected / stats.total_reads * 100 << "%)\n\n";
        
        file << "Average Tail Lengths:\n";
        file << "  PolyC length: " << stats.avg_polyc_length << "\n";
        file << "  PolyT length: " << stats.avg_polyt_length << "\n";
        file << "  CT tail length: " << stats.avg_ct_tail_length << "\n";
        file << "  PolyG length: " << stats.avg_polyg_length << "\n";
        file << "  PolyA length: " << stats.avg_polya_length << "\n";
        file << "  AG head length: " << stats.avg_ag_tail_length << "\n\n";
        
        file << "Score Distribution:\n";
        for (const auto& range : stats.score_ranges) {
            file << "  " << range.first << ": " << range.second 
                 << " (" << (double)range.second / stats.total_reads * 100 << "%)\n";
        }
    }

    void printComparison(const std::vector<std::pair<std::string, ProcessingStats>>& all_stats) {
        std::cout << "\n" << std::string(100, '=') << std::endl;
        std::cout << "COMPREHENSIVE COMPARISON SUMMARY" << std::endl;
        std::cout << std::string(100, '=') << std::endl;
    
        auto print_metric = [&](const std::string& name, std::function<double(const ProcessingStats&)> getter) {
            std::cout << std::left << std::setw(25) << name;
            for (const auto& stat_pair : all_stats) {
                std::cout << std::setw(15) << std::fixed << std::setprecision(1) << getter(stat_pair.second);
            }
            std::cout << std::endl;
        };
    
        // Print header
        std::cout << std::left << std::setw(25) << "Metric";
        for (const auto& stat_pair : all_stats) {
            std::cout << std::setw(15) << stat_pair.first;
        }
        std::cout << std::endl;
        std::cout << std::string(100, '-') << std::endl;
        
        print_metric("Trim rate (%)", [](const ProcessingStats& s) { 
            return (double)s.trimmed_reads / s.total_reads * 100; 
        });
        print_metric("R1 trim rate (%)", [](const ProcessingStats& s) { 
            return (double)s.r1_trimmed / s.total_reads * 100; 
        });
        print_metric("R2 trim rate (%)", [](const ProcessingStats& s) { 
            return (double)s.r2_trimmed / s.total_reads * 100; 
        });
        print_metric("Both trim rate (%)", [](const ProcessingStats& s) { 
            return (double)s.both_trimmed / s.total_reads * 100; 
        });
        print_metric("CT pattern (%)", [](const ProcessingStats& s) { 
            return (double)s.ct_pattern_detected / s.total_reads * 100; 
        });
        print_metric("C only (%)", [](const ProcessingStats& s) { 
            return (double)s.c_only_detected / s.total_reads * 100; 
        });
        print_metric("Mixed tail (%)", [](const ProcessingStats& s) { 
            return (double)s.mixed_tail_detected / s.total_reads * 100; 
        });
        print_metric("Avg R1 trim len", [](const ProcessingStats& s) { 
            return s.avg_r1_trim_length; 
        });
        print_metric("Avg R2 trim len", [](const ProcessingStats& s) { 
            return s.avg_r2_trim_length; 
        });
        print_metric("Merge rate (%)", [](const ProcessingStats& s) { 
            return (double)s.merged_reads / s.total_reads * 100; 
        });
    }
};

int main() {
    std::cout << "Comprehensive CT Tail Analysis - Three Datasets" << std::endl;
    std::cout << "===============================================" << std::endl;
    
    // Configure datasets
    std::vector<DatasetConfig> datasets = {
        {
            "tumor", 
            "Tumor CT Data", 
            "tumor_R1.fastq", 
            "tumor_R2.fastq",
            TrimParams::createCTConfig(25),
            "tumor_trimmed"
        },
        {
            "normal", 
            "Normal CT Data", 
            "normal_R1.fastq", 
            "normal_R2.fastq",
            TrimParams::createCTConfig(25),
            "normal_trimmed"
        },
        {
            "abclonal", 
            "ABclonal PolyC Data", 
            "ABclonal_polyC_2407_R1.fastq", 
            "ABclonal_polyC_2407_R2.fastq",
            TrimParams::createPolyCConfig(25),
            "abclonal_trimmed"
        }
    };
    
    // Adjust CT parameters for tumor and normal
    datasets[0].trim_params.t_score = 8.0;  // Tumor: T = 80% of C
    datasets[0].trim_params.consecutive_decay_rate = 1.5;
    
    datasets[1].trim_params.t_score = 8.0;  // Normal: T = 80% of C
    datasets[1].trim_params.consecutive_decay_rate = 1.5;
    
    ComprehensiveProcessor processor;
    std::vector<std::pair<std::string, ProcessingStats>> all_results;
    
    for (const auto& dataset : datasets) {
        std::cout << "\n" << std::string(80, '-') << std::endl;
        std::cout << "PROCESSING: " << dataset.description << std::endl;
        std::cout << std::string(80, '-') << std::endl;
        
        // Process first 3000 reads for comparison
        ProcessingStats stats = processor.processDataset(dataset, 3000, (dataset.name == "tumor"));
        
        all_results.push_back({dataset.name, stats});
        
        std::cout << "\nQuick Summary for " << dataset.description << ":" << std::endl;
        std::cout << "  Total reads: " << stats.total_reads << std::endl;
        std::cout << "  Trim rate: " << std::fixed << std::setprecision(1) 
                  << (double)stats.trimmed_reads / stats.total_reads * 100 << "%" << std::endl;
        std::cout << "  CT patterns: " << stats.ct_pattern_detected 
                  << " (" << (double)stats.ct_pattern_detected / stats.total_reads * 100 << "%)" << std::endl;
        std::cout << "  C only: " << stats.c_only_detected 
                  << " (" << (double)stats.c_only_detected / stats.total_reads * 100 << "%)" << std::endl;
        std::cout << "  Both trimmed: " << stats.both_trimmed 
                  << " (" << (double)stats.both_trimmed / stats.total_reads * 100 << "%)" << std::endl;
    }
    
    // Print comprehensive comparison
    processor.printComparison(all_results);
    
    std::cout << "\n" << std::string(100, '=') << std::endl;
    std::cout << "ANALYSIS COMPLETE" << std::endl;
    std::cout << "Output files generated:" << std::endl;
    std::cout << "  tumor_trimmed_[R1|R2|merged|stats].*" << std::endl;
    std::cout << "  normal_trimmed_[R1|R2|merged|stats].*" << std::endl;
    std::cout << "  abclonal_trimmed_[R1|R2|merged|stats].*" << std::endl;
    std::cout << "Check individual stats files for detailed analysis." << std::endl;
    std::cout << std::string(100, '=') << std::endl;
    
    return 0;
}
