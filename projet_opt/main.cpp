#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "database.hpp"
#include "query.hpp"
#include "protein.hpp"

using namespace std;

int main(int argc, char *argv[])
{
    // Arguments
    Query query_protein(argv[1]);
    Database database(argv[2]);
    Query query_blosum(argv[3]);
    int gap_open_penalty = stoi(argv[4]);
    int gap_extension_penalty = stoi(argv[5]);

    // Création de la matrice BLOSUM et de la protéine de requête
    vector<vector<int>> blosum = query_blosum.makeMatrixBlosum();
    Protein protein_query = query_protein.makeProteinFromQuery();

    // Alignement de la protéine de requête avec toutes les protéines de la base de données
    database.alignProteins(protein_query, blosum, gap_open_penalty, gap_extension_penalty);
    
    // Affichage des 20 protéines les plus alignées avec la protéine de requête
    vector<Protein> most_aligned_proteins = database.getMostAlignedProteins();
    for (int i = 0; i < most_aligned_proteins.size(); i++)
    {
        cout << most_aligned_proteins[i].getId() << " " << most_aligned_proteins[i].getScore() << endl;
    }
    return 0;
}