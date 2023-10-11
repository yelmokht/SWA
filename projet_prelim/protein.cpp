#include <string>

#include "protein.hpp"

using namespace std;

Protein::Protein(string id, string name, string choice, string sequence)
{
    this->id = id;             // sp|P00533|EGFR_HUMAN
    this->name = name;         // Epidermal growth factor receptor
    this->choice = choice;     // OS=Homo sapiens OX=9606 GN=EGFR PE=1 SV=2
    this->sequence = sequence; // MRPSGTAG ...
}

string Protein::getId()
{
    return this->id;
}

string Protein::getName()
{
    return this->name;
}

string Protein::getChoice()
{
    return this->choice;
}

string Protein::getSequence()
{
    return this->sequence;
}