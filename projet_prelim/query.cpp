#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstring>

#include "query.hpp"

using namespace std;

Query::Query(string file_path)
{
    this->file_path = file_path;
}

Protein Query::makeProteinFromQuery()
{
    ifstream file_stream(file_path, ios::binary);
    string id, name, choice, sequence;
    if (file_stream.is_open())
    {
        int linepos = 0;
        string buffer;
        while (getline(file_stream, buffer))
        {
            if (linepos == 0)
            { // Si on se trouve à la première ligne
                stringstream ss(buffer);
                string token;
                int pos = 0;
                while (getline(ss, token, ' '))
                { // On split la ligne par un espace et on stocke dans token
                    if (token.find('>') != string::npos)
                    {
                        id = token;
                    }
                    /* Pour partie projet si besoin
                    else
                    {
                        if (token.find('=') == string::npos)
                        {
                            name += token + " ";
                        }
                        else
                        {
                            choice += token + " ";
                        }
                    }
                    */
                }
            }
            else
            {
                sequence += buffer;
            }
            linepos++;
        }
    }
    id = id.substr(1, id.size() - 1); // On enlève '>' de l'id
    /*
    name.pop_back();                  // On enlève l'espace restant
    choice.pop_back();                // On enlève l'espace restant
    */
    Protein protein = Protein(id, name, choice, sequence);
    return protein;
}