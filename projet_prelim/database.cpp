#include <fstream>
#include <iostream>
#include <sstream>
#include <cstring>
#include <arpa/inet.h>

#include "database.hpp"

using namespace std;

Database::Database(string file_path) {
    this->file_path = file_path;
}
/*
Cette fonction récupère la position de la protéine parmi N protéines
*/
int Database::getIndex(string protein_sequence) {
    string protein_sequence_file_path = file_path + ".psq";
    ifstream sequence_stream(protein_sequence_file_path, ios::binary);
    unsigned char buffer;
    unsigned int index = 0; //Indice de la proteine dans le fichier .psq
    bool store_sequence = true; //Permission de stockage
    if (sequence_stream.is_open())
    {
        while (!sequence_stream.eof())
        {
            sequence_stream.read((char *)&buffer, sizeof(unsigned char));
            decode_char(&buffer); 
            // Only store the good sequence.
            if (store_sequence)
            {
                if (buffer != '-')  //Rejet du '-'
                {
                    seq_vect.push_back(buffer); 
                    //Subvector de la query_sequence qui permet la compairason de petit morceau de code
                    string subquery_sequence(protein_sequence.begin(), protein_sequence.begin() + seq_vect.size());
                    
                    if (subquery_sequence != seq_vect)
                    {   //Si le seq_vect et subquery divergent, il faut plus stocker et il faut nettoyer seq_vect
                        store_sequence = false;
                        seq_vect.clear();
                    }
                    else if (seq_vect == protein_sequence)
                    {   //Si le seq_vect correspond à la query_sequence
                        break;
                    }
                }
            }
            else
            {
                if (buffer == '-')
                {   //Si on se trouve à la fin de la séquence courant, alors l'indice incrémente et on peut à nouveau
                    //stocker les valeurs dans seq_vect
                    index++;
                    store_sequence = true;
                }
            }
        }
        sequence_stream.close();
    }
    else
    {
        cerr << "Couldn't open the .psq file" << endl;
        exit(1);
    }
    return index;
}

/*
Cette fonction retourne l'acide aminé selon une certaine valeur définie par le fichier .psq
*/
void Database::decode_char(unsigned char *buffer) {
    const char encode_sequence[29] = "-ABCDEFGHIKLMNPQRSTVWXYZU*OJ";
    *buffer = encode_sequence[*buffer];
}

/*
Cette fonction retourne l'offset (= position) de la protéine dans le fichier .pin
*/
int Database::getOffset(int index) {
    string index_file_path = file_path + ".pin";
    ifstream index_file_stream(index_file_path, ios::binary);
    int offset;
    if (index_file_stream.is_open()) {
        uint32_t buffer;
        index_file_stream.seekg((22 + index) * sizeof(uint32_t)); //On change de position : on lit à partir de 22 (début du header offset table) + index
        index_file_stream.read((char *)&buffer, sizeof(uint32_t));
        offset = conversion_endian_btl(buffer);
    } else {
        cerr << "Couldn't open the file" << endl;
    }
    return offset;
}

/*
Cette fonction convertit une valeur big endian en little endian
*/
uint32_t Database::conversion_endian_btl(uint32_t value) {
    return ntohl(value);
}

/*
Cette fonction retourne la protéine se trouvant dans la base de données
à partir du offset du fichier .pin
*/
Protein Database::getProtein(int header_offset) {
    // Reading .phr file
    string header_file_path = file_path + ".phr";
    ifstream header_stream(header_file_path, ios::binary);
    if (header_stream.is_open())
    {
        int current_offset = 0; //Donne la position actuelle dans le fichier
        int store_hex = 0;  //Permission de stocker ou pas le byte lu
        unsigned char buffer;   
        while (!header_stream.eof())
        {
            header_stream.read((char *)&buffer, sizeof(unsigned char));
            if (current_offset >= header_offset) //Si on se trouve au bon endroit
            {
                if (buffer == 0x1A) // S'il s'agit d'un string
                {
                    store_hex = 1; // Commencez à stocker les hex
                }
                else if (buffer == 0x20 && store_hex) // S'il s'agit d'un espace et qu'on stockait des données
                {
                    break;
                }
                if (store_hex) //Permission de stocker des données
                {
                    hex_vect.push_back(buffer);
                }
            }
            current_offset++; // Avancez dans la lecture
        }
        header_stream.close();
    }
    else
    {
        cerr << "Couldn't read the .phr file" << endl;
    }
    // Parsing hex vector
    if (hex_vect.size() > 0)    //Pour éviter des erreurs de segmentation
    {
        hex_vect.erase(hex_vect.begin(), hex_vect.begin() + 2); // Rejet du token '1A' et de la taille
    }
    string ret_string;
    for (int i = 0; i < hex_vect.size(); i++)
    {
        ret_string += hex_vect[i]; // Concatenation de tous les hex pour former un seul string
    }
    Protein protein = Protein(ret_string, "","",""); //TODO
    return protein;
}