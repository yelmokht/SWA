#pragma once
#include <vector>
#include <string>
#include <iostream>

#include "protein.hpp"

class Database
{
private:
    std::string file_path;
    int NUMBER_PROTEIN;
    int header_offset_start = 0;
    int sequence_offset_start;
    std::vector<Protein> protein_vector;
    
public:
    Database(std::string file_path);
    void searchNumberProteins();
    uint32_t conversion_endian_btl(uint32_t value);
    void searchStartOffset();
    int getOffset(int index, std::string offset_type);
    std::string getProteinId(int header_offset);
    unsigned char decode_char(unsigned char c);
    std::string getProteinSequence(int sequence_offset);
    void alignProteins(Protein protein, std::vector<std::vector<int>> &blosum, int gap_open_penalty, int gap_extension_penalty);
    static void align_thread(Protein &protein, Protein &protein_query, std::vector<std::vector<int>> &blosum, int gap_open_penalty, int gap_extension_penalty);
    std::vector<Protein> getMostAlignedProteins();
};