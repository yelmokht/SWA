#include <iostream>
#include <fstream>
#include <string>

#include "database.hpp"
#include "query.hpp"
#include "protein.hpp"

using namespace std;

int main(int argc, char *argv[])
{
    Query query(argv[1]);
    Database db(argv[2]);
    Protein protein_query = query.makeProteinFromQuery();
    int index_protein = db.getIndex(protein_query.getSequence());
    int offset_protein = db.getOffset(index_protein);
    Protein protein_database = db.getProtein(offset_protein);
    cout << protein_database.getId();
    /*
    cout << "Proteine de requête : " << protein_query.getId() << endl << "Protéine de la base de données : " << protein_database.getId() << endl;
    */
    return 0;
}