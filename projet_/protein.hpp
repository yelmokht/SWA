#ifndef _PROTEIN_HPP
#define _PROTEIN_HPP
#include <string>
#include <vector>

class Protein
{
private:
    std::string id, sequence;
    int score = 0;

public:
    Protein(std::string id, std::string sequence);
    std::string getId();
    std::string getSequence();
    int getScore();
    int substitution_matrix(std::vector<std::vector<int>> &blosum, char a, char b);
    int gap_penalty(int gap_open_penalty, int gap_extension_penalty, int length);
    void calculate_score(std::string sequence, std::vector<std::vector<int>> &blosum, int gap_open_penalty, int gap_extension_penalty);
};
#endif