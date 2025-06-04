#include <iostream>     
#include <fstream>     
#include <string>       
#include <vector>       
#include <map>          
#include <algorithm>    
#include <sstream>      
#include <stdexcept>    
#include <filesystem>   
#include <cctype>       
#include <cstdint>  
#include <chrono>     


//enum class for representing DNA bases A, C, G, T
enum class Base : uint8_t {
    A = 0, C = 1, G = 2, T = 3, UNKNOWN = 4
};

//struct for auxiliary information from target fasta file
struct TargetAuxInfo {

    std::string id;
    std::vector<int> line_lengths;
    std::vector<std::pair<size_t, char>> non_acgt_chars;
};

//struct for representing a found refrence sequence or a mismatched literal
struct MatchRecord {
    
    bool is_match;
    size_t ref_start_pos;
    int length;
    std::vector<Base> mismatch_sequence;

};

//converts ints to base enum value
Base int_to_base(int val) {
    if (val >= 0 && val <= 3) return static_cast<Base>(val);
    return Base::UNKNOWN; 
}

//converts chars to base enum value
Base char_to_base(char ch) {
    if (ch == 'A') return Base::A;
    if (ch == 'C') return Base::C;
    if (ch == 'G') return Base::G;
    if (ch == 'T') return Base::T;
    return Base::UNKNOWN;
}

//converts base enum value to chars
char base_to_char(Base b) {
    if (b == Base::A) return 'A';
    if (b == Base::C) return 'C';
    if (b == Base::G) return 'G';
    if (b == Base::T) return 'T';
    return '?';
}

//function for fasta file parsing which takes file_path as an input and outputs ACGT sequence
void parse_reference_fasta(const std::string& file_path, std::vector<Base>& out_pure_acgt_sequence) {
    
    std::ifstream file(file_path);

    if (!file.is_open()) {
        throw std::runtime_error("Cannot open FASTA file: " + file_path);
    }

    out_pure_acgt_sequence.clear(); 

    std::string line;
    bool first_header_processed = false; 

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '\r' || line[0] == '\n') continue; 

        if (line[0] == '>') {

            if (!first_header_processed) {
                first_header_processed = true;
            }
            continue; 
            
        }

        if (!line.empty() && line.back() == '\r') line.pop_back(); 
        if (line.empty()) continue;

        for (char c : line) { 

            Base b = char_to_base(c);
            if(b != Base::UNKNOWN) {
                out_pure_acgt_sequence.push_back(b);
            }
            
        }
    }
}

//function that deserializes Run-Length Encoding format input 
void deserialize_rle_vector(std::istream& in, std::vector<int>& data) {

    data.clear();
    int rle_count;
    in >> rle_count; 

    if (rle_count == 0) { 
        return;
    }

    int value, count;
    for (int i = 0; i < rle_count / 2; ++i) { 

        in >> value >> count; 
        for (int k = 0; k < count; ++k) {
            data.push_back(value);
        }
        
    }

    std::string rest_of_the_line; 
    std::getline(in, rest_of_the_line);
}


//function for deserializing compressed data from input file
void deserialize_data_compact(std::istream& in,TargetAuxInfo& aux_info, std::vector<MatchRecord>& match_records) {
    
    aux_info = TargetAuxInfo(); 
    match_records.clear();      
    std::string line;          

    //header line
    std::getline(in, line);
    aux_info.id = line.substr(1);

    //RLE line lengths
    deserialize_rle_vector(in, aux_info.line_lengths);

    //deserialize Non-ACGT characters
    int num_non_acgt_occurrences;
    in >> num_non_acgt_occurrences;

    size_t last_non_acgt_pos = 0;
    for (int i = 0; i < num_non_acgt_occurrences; ++i) {
        size_t delta_pos;
        int char_as_int;

         if (!(in >> delta_pos >> char_as_int)) {
            throw std::runtime_error("Error reading non-ACGT char position or value at occurrence " + std::to_string(i+1));
        }

        last_non_acgt_pos += delta_pos;
        aux_info.non_acgt_chars.push_back({last_non_acgt_pos, static_cast<char>(char_as_int)});
    }
    std::getline(in, line); 
    

    //match/mismatch data
    size_t last_ref_match_pos = 0; 
    while (std::getline(in, line)) { 

        if (line.empty()) continue; 

        std::istringstream iss(line);
        char type_char;
        iss >> type_char; 

        MatchRecord rec; 
        if (type_char == 'M') {

            //match 
            rec.is_match = true;
            size_t ref_pos;
            iss >> ref_pos >> rec.length; 
            rec.ref_start_pos = last_ref_match_pos+ ref_pos; 
            last_ref_match_pos = rec.ref_start_pos; 

        } else if (type_char == 'S') {
            
            //mismatch
            rec.is_match = false;
            iss >> rec.length; 
            std::string base_digits_str;
            iss >> base_digits_str; 
            
            rec.mismatch_sequence.reserve(rec.length);
            for (char digit : base_digits_str) { 
                rec.mismatch_sequence.push_back(int_to_base(digit - '0'));
            
            }
        }
        match_records.push_back(rec); 
    }
}

//function for reconstructing fasta file
void reconstruct_and_write_fasta(std::ostream& out_stream, const TargetAuxInfo& aux_info, const std::vector<MatchRecord>& match_records, const std::vector<Base>& ref_pure_acgt_sequence) {
    
    //reconstruct ACGT sequence
    std::vector<Base> target_pure_acgt_sequence;
    size_t estimated_target_pure_len = 0;
    for(const auto& rec : match_records) estimated_target_pure_len += rec.length;
    target_pure_acgt_sequence.reserve(estimated_target_pure_len);

    for (const auto& rec : match_records) {
        if (rec.is_match) {
            for (int i = 0; i < rec.length; ++i) {
                target_pure_acgt_sequence.push_back(ref_pure_acgt_sequence[rec.ref_start_pos + i]);
            }
        } else {
            target_pure_acgt_sequence.insert(target_pure_acgt_sequence.end(), rec.mismatch_sequence.begin(), rec.mismatch_sequence.end());
        }
    }

    std::string final_sequence_str;
    size_t final_estimated_length = target_pure_acgt_sequence.size() + aux_info.non_acgt_chars.size();
    final_sequence_str.reserve(final_estimated_length);

    //create a map of original positions for non-ACGT Chars
    std::map<size_t, char> non_acgt_insertions_map;
    for (const auto& item : aux_info.non_acgt_chars) {
        non_acgt_insertions_map[item.first] = item.second;
    }

    size_t pure_acgt_idx = 0;         
    size_t original_pos_tracker = 0; 

    size_t total_original_length_from_lines = 0;
    for(int len : aux_info.line_lengths) total_original_length_from_lines += len;

    size_t total_length_to_build = total_original_length_from_lines > 0 ? total_original_length_from_lines : final_estimated_length;

    if (total_length_to_build == 0 && !target_pure_acgt_sequence.empty() && aux_info.non_acgt_chars.empty()){
        total_length_to_build = target_pure_acgt_sequence.size();
    }

    //build the string character by character
    for (original_pos_tracker = 0; original_pos_tracker < total_length_to_build; ++original_pos_tracker) {

        auto it = non_acgt_insertions_map.find(original_pos_tracker);
        if (it != non_acgt_insertions_map.end()) { 
            
            final_sequence_str += it->second;

        } else { 

            if (pure_acgt_idx < target_pure_acgt_sequence.size()) {
                final_sequence_str += base_to_char(target_pure_acgt_sequence[pure_acgt_idx++]);
            }
        }
    }
    
    //write to output
    out_stream << ">" << aux_info.id << "\n";
    size_t current_char_idx_output = 0;
    if (!aux_info.line_lengths.empty()) {
        for (int line_len : aux_info.line_lengths) {
            if (current_char_idx_output >= final_sequence_str.length()) break;
            size_t len_to_write = std::min(static_cast<size_t>(line_len), final_sequence_str.length() - current_char_idx_output);
            out_stream << final_sequence_str.substr(current_char_idx_output, len_to_write) << "\n";
            current_char_idx_output += len_to_write;
        }
        if (current_char_idx_output < final_sequence_str.length()){
             out_stream << final_sequence_str.substr(current_char_idx_output) << "\n";
        }
    }
}

//main program function that takes refrence fasta file argument and compressed txt argument
int main(int argc, char* argv[]) {
    
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <reference.fasta> <compressed_file.txt>\n";
        return 1; 
    }

    auto start_time = std::chrono::high_resolution_clock::now();

    std::string ref_fasta_path = argv[1];        
    std::string compressed_file_path = argv[2];
    std::filesystem::path comp_p(compressed_file_path);  
    std::string output_fasta_path = (comp_p.parent_path() / ("decompressed.fna")).string();

    try {

        //parse reference fasta file
        std::cout << "Parsing reference: " << ref_fasta_path << "..." << std::endl;
        std::vector<Base> ref_pure_acgt_sequence; 
        parse_reference_fasta(ref_fasta_path, ref_pure_acgt_sequence);
       
        if (ref_pure_acgt_sequence.empty()) { 
            throw std::runtime_error("Reference sequence is empty or contains no ACGT characters for decompression.");
        }

        //deserialize the data from the compact compressed file
        std::cout << "Reading and deserializing compressed file: " << compressed_file_path << "..." << std::endl;
        std::ifstream in_compressed_file(compressed_file_path); 

        if (!in_compressed_file.is_open()) { 
            throw std::runtime_error("Cannot open compressed file: " + compressed_file_path);
        }
        TargetAuxInfo target_aux_info; 
        std::vector<MatchRecord> match_records; 
        deserialize_data_compact(in_compressed_file, target_aux_info, match_records);
        in_compressed_file.close(); 
        std::cout << "Deserialization complete. Target ID: " << target_aux_info.id << std::endl;

        //Reconstruct the original target fasta sequence and write it to the output file
        std::cout << "Reconstructing and writing FASTA to: " << output_fasta_path << "..." << std::endl;
        std::ofstream out_fasta_file(output_fasta_path); 
        
        if (!out_fasta_file.is_open()) { 
            throw std::runtime_error("Cannot open output FASTA file: " + output_fasta_path);
        }
        reconstruct_and_write_fasta(out_fasta_file, target_aux_info, match_records, ref_pure_acgt_sequence);
        out_fasta_file.close(); 

        auto end_time = std::chrono::high_resolution_clock::now(); 
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        std::cout << "Decompression successful." << std::endl;
        std::cout << "Total execution time: " << duration.count() << " ms" << std::endl;

    } catch (const std::exception& e) { 
        std::cerr << "Error: " << e.what() << std::endl; 
        return 1; 
    }

    return 0; 
}