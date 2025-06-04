#include <iostream>    
#include <fstream>      
#include <string>       
#include <vector>              
#include <algorithm>    
#include <cmath>        
#include <cstdint>      
#include <filesystem>   
#include <chrono> 

const int K_TUPLE_LENGTH = 20;

const size_t HASH_TABLE_LENGTH = 1 << 25;

enum class Base : uint8_t {
    A = 0, C = 1, G = 2, T = 3, UNKNOWN = 4
};

struct TargetAuxInfo {

    std::string id;
    std::vector<int> line_lengths;
    std::vector<std::pair<size_t, char>> non_acgt_chars;
};

struct MatchRecord {
    
    bool is_match;
    size_t ref_start_pos;
    int length;
    std::vector<Base> mismatch_sequence;

    MatchRecord(size_t r_pos, int l) : is_match(true), ref_start_pos(r_pos), length(l) {}

    MatchRecord(const std::vector<Base>& seq) : is_match(false), ref_start_pos(0), length(static_cast<int>(seq.size())), mismatch_sequence(seq) {}
};

Base char_to_base(char ch) {;
    if (ch == 'A') return Base::A;
    if (ch == 'C') return Base::C;
    if (ch == 'G') return Base::G;
    if (ch == 'T') return Base::T;
    return Base::UNKNOWN;
}


void parse_fasta(const std::string& file_path, TargetAuxInfo& aux_info, std::vector<Base>& out_pure_acgt_sequence, bool is_target_file) {

    std::ifstream file(file_path);

    if (!file.is_open()) {
        throw std::runtime_error("Cannot open FASTA file: " + file_path);
    }

    out_pure_acgt_sequence.clear(); 
    if (is_target_file) {
        aux_info = TargetAuxInfo(); 
    }

    std::string line; 
    bool first_header = true; 
    size_t original_char_position = 0;

    
    while (std::getline(file, line)) {

        if (line.empty() || line[0] == '\r' || line[0] == '\n') continue;

        if (line[0] == '>') { 

            if (first_header) {

                if (is_target_file) {
                    size_t id_end = line.find_first_of(" \t\r\n");
                    aux_info.id = line.substr(1, id_end != std::string::npos ? id_end - 1 : std::string::npos);
                }
                first_header = false;

            }
        } else { 
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.empty()) continue;

            if (is_target_file) {
                aux_info.line_lengths.push_back(line.length()); 
            }

            for (char c : line) { 

                Base b = char_to_base(c); 
                if (b != Base::UNKNOWN) { 
                    out_pure_acgt_sequence.push_back(b);
                } else { 
                    if (is_target_file) {
                        aux_info.non_acgt_chars.push_back({original_char_position, c});
                    }
                }
                original_char_position++;

            }
        }
    }
}

uint64_t calculate_kmer_value(const std::vector<Base>& seq, size_t start_index, int k_t_l) {

    uint64_t value = 0;
    uint64_t multiplier = 1; 
    for (int i = 0; i < k_t_l; ++i) {

        if (start_index + i >= seq.size()) { 
             throw std::out_of_range("k-mer calculation out of bounds");
        }

        value += static_cast<uint64_t>(seq[start_index + i]) * multiplier;

        if (i < k_t_l - 1) {  
            multiplier *= 4;
        }
    }
    return value;
}

using HashTable = std::vector<std::vector<size_t>>;

void build_reference_hash_table(const std::vector<Base>& ref_pure_acgt_sequence, HashTable& hash_table) {

    hash_table.assign(HASH_TABLE_LENGTH, std::vector<size_t>()); 
    if (ref_pure_acgt_sequence.size() < static_cast<size_t>(K_TUPLE_LENGTH)) { 
        return;
    }

    for (size_t i = 0; i <= ref_pure_acgt_sequence.size() - K_TUPLE_LENGTH; ++i) {
        uint64_t kmer_val = calculate_kmer_value(ref_pure_acgt_sequence, i, K_TUPLE_LENGTH); 
        size_t hash_table_index = kmer_val % HASH_TABLE_LENGTH; 
        hash_table[hash_table_index].push_back(i); 
    }
}

void greedy_match_target(const std::vector<Base>& target_pure_acgt_sequence, const std::vector<Base>& ref_pure_acgt_sequence, const HashTable& ref_hash_table, std::vector<MatchRecord>& out_match_records) {
   
    out_match_records.clear(); 
    if (target_pure_acgt_sequence.empty()) return; 

    size_t current_target_pos = 0; 
    std::vector<Base> current_mismatch_buffer; 


    while (current_target_pos < target_pure_acgt_sequence.size()) {

        if (target_pure_acgt_sequence.size() - current_target_pos < static_cast<size_t>(K_TUPLE_LENGTH)) {

            for (size_t i = current_target_pos; i < target_pure_acgt_sequence.size(); ++i) {
                current_mismatch_buffer.push_back(target_pure_acgt_sequence[i]);
            }

            current_target_pos = target_pure_acgt_sequence.size(); 
            break; 
        }


        uint64_t target_kmer_val = calculate_kmer_value(target_pure_acgt_sequence, current_target_pos, K_TUPLE_LENGTH);
        size_t hash_table_index = target_kmer_val % HASH_TABLE_LENGTH; 

        int best_match_len = 0; 
        size_t best_match_ref_start= 0; 

        if (hash_table_index < ref_hash_table.size()) { 

            for (size_t ref_kmer_start_pos : ref_hash_table[hash_table_index]) { 
                
                bool direct_kmer_match = true;
                
                if (ref_kmer_start_pos + K_TUPLE_LENGTH > ref_pure_acgt_sequence.size()) continue;

                for (int k_offset = 0; k_offset < K_TUPLE_LENGTH; ++k_offset) {

                    if (ref_pure_acgt_sequence[ref_kmer_start_pos + k_offset] != target_pure_acgt_sequence[current_target_pos + k_offset]) {

                        direct_kmer_match = false; 
                        break;

                    }
                }

                if (direct_kmer_match) { 
                    int current_match_len = K_TUPLE_LENGTH;
                    
                    while (ref_kmer_start_pos + current_match_len < ref_pure_acgt_sequence.size() &&
                           current_target_pos + current_match_len < target_pure_acgt_sequence.size() &&
                           ref_pure_acgt_sequence[ref_kmer_start_pos + current_match_len] ==
                           target_pure_acgt_sequence[current_target_pos + current_match_len]) {

                        current_match_len++;
                    }
                    
                    if (current_match_len > best_match_len) {

                        best_match_len = current_match_len;
                        best_match_ref_start= ref_kmer_start_pos;

                    }
                }
            }
        }

        
        if (best_match_len >= K_TUPLE_LENGTH) { 
            if (!current_mismatch_buffer.empty()) {
                out_match_records.emplace_back(current_mismatch_buffer);
                current_mismatch_buffer.clear();
            }
            
            out_match_records.emplace_back(best_match_ref_start, best_match_len);
            current_target_pos += best_match_len;
            
        } else { 
            current_mismatch_buffer.push_back(target_pure_acgt_sequence[current_target_pos]);
            current_target_pos++; 
        }
    }

    if (!current_mismatch_buffer.empty()) {
        out_match_records.emplace_back(current_mismatch_buffer);
    }
}

void write_rle_vector(std::ostream& out, const std::vector<int>& data) {

    if (data.empty()) {
        out << "0";
        return;
    }

    std::vector<int> rle_sequence; 
    rle_sequence.push_back(data[0]); 
    int current_run_count = 1;      

    
    for (size_t i = 1; i < data.size(); ++i) {
        if (data[i] == data[i-1]) { 
            current_run_count++;    
        } else { 
            rle_sequence.push_back(current_run_count); 
            rle_sequence.push_back(data[i]);         
            current_run_count = 1;                   
        }
    }
    rle_sequence.push_back(current_run_count); 

    out << rle_sequence.size(); 
    for (int val : rle_sequence) { 
        out << " " << val;
    }
}


void serialize_compressed_data(std::ostream& out, const TargetAuxInfo& aux_info, const std::vector<MatchRecord>& match_records) {
    
    out << ">" << aux_info.id << "\n";

    write_rle_vector(out, aux_info.line_lengths);
    out << "\n";

    out << aux_info.non_acgt_chars.size(); 
    size_t last_non_acgt_pos = 0;
    for (const auto& non_acgt_item : aux_info.non_acgt_chars) {
        out << " " << (non_acgt_item.first - last_non_acgt_pos); 
        out << " " << static_cast<int>(non_acgt_item.second);   
        last_non_acgt_pos = non_acgt_item.first;
    }
    out << "\n";


    size_t last_ref_match_pos = 0; 
    for (const auto& rec : match_records) {
        if (rec.is_match) {
            out << "M " << (rec.ref_start_pos - last_ref_match_pos); 
            out << " " << rec.length << "\n";                                     
            last_ref_match_pos = rec.ref_start_pos; 
        } else {
            out << "S " << rec.length << " "; 
            for (Base b : rec.mismatch_sequence) {
                out << static_cast<int>(b); 
            }
            out << "\n";
        }
    }
}

int main(int argc, char* argv[]) {
    
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <reference.fasta> <target.fasta>\n";
        return 1; 
    }

    auto start_time = std::chrono::high_resolution_clock::now();

    std::string ref_fasta_path = argv[1];  
    std::string target_fasta_path = argv[2]; 

    try {
        std::cout << "Parsing reference: " << ref_fasta_path << "..." << std::endl;
        TargetAuxInfo ref_aux_info;
        std::vector<Base> ref_pure_acgt_sequence; 
        parse_fasta(ref_fasta_path, ref_aux_info, ref_pure_acgt_sequence, false);
        
        if (ref_pure_acgt_sequence.empty()) { 
            throw std::runtime_error("Reference sequence is empty or contains no ACGT characters.");
        }
        std::cout << "Building hash table from reference..." << std::endl;
        HashTable reference_hash_table; 
        build_reference_hash_table(ref_pure_acgt_sequence, reference_hash_table);
        

        std::cout << "Parsing target: " << target_fasta_path << "..." << std::endl;
        TargetAuxInfo target_aux_info; 
        std::vector<Base> target_pure_acgt_sequence;
        parse_fasta(target_fasta_path, target_aux_info, target_pure_acgt_sequence, true);
       
        std::cout << "Performing greedy matching..." << std::endl;
        std::vector<MatchRecord> match_records; 
        greedy_match_target(target_pure_acgt_sequence, ref_pure_acgt_sequence, reference_hash_table, match_records);

        std::filesystem::path target_p(target_fasta_path); 
        std::string output_filename_str = (target_p.parent_path() / ("compressed_" + target_p.stem().string() + ".txt")).string();
        
        std::cout << "Writing compressed data to: " << output_filename_str << "..." << std::endl;
        std::ofstream out_file(output_filename_str); 
        if (!out_file.is_open()) { 
            throw std::runtime_error("Cannot open output file: " + output_filename_str);
        }

        serialize_compressed_data(out_file, target_aux_info, match_records);
        out_file.close(); 

        auto end_time = std::chrono::high_resolution_clock::now(); 
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        std::cout << "Compression successful." << std::endl;
        std::cout << "Total execution time: " << duration.count() << " ms" << std::endl;
        

    } catch (const std::exception& e) { 
        std::cerr << "Error: " << e.what() << std::endl; 
        return 1; 
    }

    

    return 0; 
}