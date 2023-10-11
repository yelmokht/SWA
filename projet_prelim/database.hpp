#pragma once
#include <vector>
#include <string>
#include <iostream>

#include "protein.hpp"

class Database
{
private:
    std::string file_path;
    std::string seq_vect;
    std::vector<unsigned char> hex_vect;
    const static int NUMBER_PROTEIN = 568363;

public:
    Database(std::string file_path);
    void updateNUMBEROFPROTEIN();
    int getIndex(std::string protein_sequence);
    void decode_char(unsigned char *buffer);
    int getOffset(int index);
    uint32_t conversion_endian_btl(uint32_t value);
    Protein getProtein(int header_offset);
};