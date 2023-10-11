#include <string>
#include <vector>
#include <iostream>

#include "protein.hpp"

using namespace std;

Protein::Protein(string id, string sequence)
{
    this->id = id;             // sp|P00533|EGFR_HUMAN
    this->sequence = sequence; // MRPSGTAG ...
}

string Protein::getId()
{
    return this->id;
}

string Protein::getSequence()
{
    return this->sequence;
}

int Protein::getScore()
{
    return this->score;
}

/**
 * Cette fonction retourne la valeur de la matrice BLOSUM en fonction de chaque acide aminé de deux protéines alignés
 * @param blosum La matrice BLOSUM
 * @param a L'acide aminé de la protéine 1
 * @param b L'acide aminé de la protéine 2
 * @return La valeur de la matrice BLOSUM en fonction de chaque acide aminé de deux protéines alignés
 */
int Protein::substitution_matrix(vector<vector<int>> &blosum, char a, char b)
{
    string amino_acids = "ARNDCQEGHILKMFPSTWYVBZX*";
    return blosum[amino_acids.find(a)][amino_acids.find(b)];
}

/**
 * Cette fonction retourne le gap penalty en fonction du gap open penalty et du gap extension penalty
 * @param gap_open_penalty Le gap open penalty
 * @param gap_extension_penalty Le gap extension penalty
 * @param length Longueur du décalage
 * @return Le gap penalty en fonction du gap open penalty, du gap extension penalty et de la longueur du décalage
 */
int Protein::gap_penalty(int gap_open_penalty, int gap_extension_penalty, int length)
{
    return gap_extension_penalty * length + gap_open_penalty;
}

/**
 * Cette fonction calcule le score d'une protéine aligné avec une protéine de requête via l'algorithme de Smith-Waterman
 * @param sequence La séquence de la protéine de requête
 * @param blosum La matrice BLOSUM
 * @param gap_open_penalty Le gap open penalty
 * @param gap_extension_penalty Le gap extension penalty
 */
void Protein::calculate_score(string sequence, vector<vector<int>> &blosum, int gap_open_penalty, int gap_extension_penalty)
{
    vector<vector<int>> matrix(sequence.size() + 1, vector<int>(this->sequence.size() + 1, 0));
    vector<int> max_values_lines(sequence.size() + 1, 0);
    vector<int> max_values_columns(this->sequence.size() + 1, 0);
    for (int line = 0; line < sequence.size() + 1; line++)
    {
        for (int column = 0; column < this->sequence.size() + 1; column++)
        {
            if (line == 0 | column == 0)
            {
                matrix[line][column] = 0;
            }
            else
            {
                int h1 = matrix[line - 1][column - 1] + substitution_matrix(blosum, sequence[line - 1], this->sequence[column - 1]);
                int h2 = matrix[line][max_values_lines[line]] - gap_penalty(gap_open_penalty, gap_extension_penalty, column - max_values_lines[line]);
                int h3 = matrix[max_values_columns[column]][column] - gap_penalty(gap_open_penalty, gap_extension_penalty, line - max_values_columns[column]);
                int h4 = 0;
                matrix[line][column] = max(max(h1, h2), max(h3, h4));

                if (matrix[line][column] > score)
                {
                    score = matrix[line][column];
                }

                if (matrix[line][column] - gap_penalty(gap_open_penalty, gap_extension_penalty, 1) >= h2)
                {
                    max_values_lines[line] = column;
                }

                if (matrix[line][column] - gap_penalty(gap_open_penalty, gap_extension_penalty, 1) >= h3)
                {
                    max_values_columns[column] = line;
                }
            }
        }
    }
}