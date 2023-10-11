#pragma once
#include <string>

class Protein
{
private:
    std::string id, name, choice, sequence;

public:
    Protein(std::string id, std::string name, std::string choice, std::string sequence);
    std::string getId();
    std::string getName();
    std::string getChoice();
    std::string getSequence();
};