#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstring>
#include <vector>
#include <algorithm>

#include "query.hpp"

using namespace std;

/**
 * Constructeur de Query
 * @param file_path Le chemin de du fichier à traiter
 */
Query::Query(string file_path)
{
    this->file_path = file_path;
}

/**
 * Cette fonction retourne une matrice BLOSUM à partir du fichier BLOSUM62 contenant les valeurs de la matrice
 * @return La matrice BLOSUM
 */
vector<vector<int>> Query::makeMatrixBlosum()
{
    ifstream file_stream(file_path);
    vector<vector<int>> blosum;
    if (file_stream.is_open())
    {
        int line = 0;
        string buffer;
        while (getline(file_stream, buffer))
        {
            if (buffer.find('#') == string::npos)
            {
                if (line > 6)
                {
                    stringstream ss(buffer);
                    string value;
                    vector<int> vector;
                    int column = 0;
                    while (getline(ss, value, ' '))
                    { // On split la ligne par un espace et on stocke dans value
                        if (column > 0 & value != "")
                        { // On enlève les espaces restants
                            vector.push_back(stoi(value));
                        }
                        column++;
                    }
                    blosum.push_back(vector);
                }
            }
            line++;
        }
    }
    return blosum;
}

/**
 * Cette fonction crée un objet Protein à partir d'un fichier .fasta contenant les infos d'une protéine
 * @return La protéine de requête
 */
Protein Query::makeProteinFromQuery()
{
    ifstream file_stream(file_path);
    string id, sequence;
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
    Protein protein = Protein(id, sequence);
    return protein;
}