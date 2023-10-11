#include <fstream>
#include <iostream>
#include <sstream>
#include <cstring>
#include <algorithm>
#include <arpa/inet.h>
#include <thread>

#include "database.hpp"

using namespace std;

/**
 * Constructeur de Database
 * @param file_path Le chemin de la base de données FASTA
 */
Database::Database(string file_path)
{
    this->file_path = file_path;
}

/**
 * Cette fonction convertit une valeur big endian en little endian.
 * @param value Une valeur en big endian
 * @return Cette même valeur en little endian
 */
uint32_t Database::conversion_endian_btl(uint32_t value)
{
    return ntohl(value);
}

/**
 * Cette fonction cherche le nombre de protéines contenu dans la base de données à partir du fichier .pjs.
 * Une fois trouvée, la valeur NUMBER_PROTEIN est modifiée.
 */
void Database::searchNumberProteins()
{
    string info_file_path = file_path + ".pjs";
    ifstream info_file_stream(info_file_path, ios::binary);
    if (info_file_stream.is_open())
    {
        string buffer;
        string number;
        while (getline(info_file_stream, buffer, ','))
        {
            if (buffer.find("number-of-sequences") != string::npos)
            {
                number = buffer.substr(buffer.find(":") + 2, buffer.size() - buffer.find(":"));
                break;
            }
        }
        NUMBER_PROTEIN = atoi(number.c_str());
    }
    info_file_stream.close();
}

/**
 * Cette fonction cherche le header_offset_start et sequence_offset_start. Une fois calculé, le header/sequence_offset_start
 * peut être réutilisé pour lire directement à partir du début du header/sequence offset table.
 */
void Database::searchStartOffset()
{
    string index_file_path = file_path + ".pin";
    ifstream index_file_stream(index_file_path, ios::binary);
    if (index_file_stream.is_open())
    {
        int i = 0;
        bool checkzero = false;
        uint32_t buffer;
        while (index_file_stream.read((char *)&buffer, sizeof(uint32_t)))
        {
            buffer = conversion_endian_btl(buffer);
            if (buffer == NUMBER_PROTEIN)
            {
                checkzero = true;
            }
            if (checkzero && buffer == 0)
            {
                i++;
                if (i == 2)
                {
                    header_offset_start = (index_file_stream.tellg() / sizeof(uint32_t)) - 1;
                    sequence_offset_start = header_offset_start + NUMBER_PROTEIN + 1;
                    break;
                }
            }
        }
    }
    index_file_stream.close();
}

/**
 * Cette fonction retourne le header/sequence offset de la protéine se trouvant le fichier .pin
 * à partir de l'index (= numéro de la protéine) et du type (header/sequence).
 * @param index L'indice de la protéine dans la base de données
 * @param offset_type Le type d'offset (header ou sequence)
 * @return Le header/sequence offset de la protéine
 */
int Database::getOffset(int index, string offset_type)
{
    // Si on lit le fichier pour la première fois, on cherche le début de la table header/sequence
    // Ceci permet de stocker le position du début de deux tables pour les prochaines lectures
    if (header_offset_start == 0)
    {
        searchStartOffset();
    }
    string index_file_path = file_path + ".pin";
    ifstream index_file_stream(index_file_path, ios::binary);
    int offset;
    uint32_t buffer;
    if (index_file_stream.is_open())
    {
        if (offset_type == "header")
        { // On change de position : on lit à partir du début du header offset table + index
            index_file_stream.seekg((header_offset_start + index) * sizeof(uint32_t));
        }
        else if (offset_type == "sequence")
        { // On change de position : on lit à partir du début du sequence offset table + index
            index_file_stream.seekg((sequence_offset_start + index) * sizeof(uint32_t));
        }
        else
        {
            cerr << "Type is incorrect: must be header or sequence" << endl;
            exit(1);
        }
        index_file_stream.read((char *)&buffer, sizeof(uint32_t));
        offset = conversion_endian_btl(buffer);
        index_file_stream.close();
    }
    else
    {
        cerr << "Couldn't open the file" << endl;
        exit(1);
    }
    return offset;
}

/**
 * Cette fonction retourne l'id d'une protéine se trouvant dans le fichier .phr pour un header_offset donné.
 * @param header_offset Le header offset de la protéine
 * @return L'id de la protéine
 */
string Database::getProteinId(int header_offset)
{
    string header_file_path = file_path + ".phr";
    ifstream header_stream(header_file_path, ios::binary);
    string id;
    if (header_stream.is_open())
    {
        unsigned char buffer;
        header_stream.seekg(header_offset * sizeof(unsigned char)); // On se déplace directement à la position pour le header_offset donné
        while (header_stream.read((char *)&buffer, sizeof(unsigned char)))
        {
            if (buffer == ' ') // Un espace signifie qu'on se trouve à la fin de l'id -> on break
            {
                break;
            }
            id.push_back(buffer);
        }
        header_stream.close();
    }
    else
    {
        cerr << "Couldn't read the .phr file" << endl;
        exit(1);
    }

    char c = id[id.find('|') - 1];
    if (c == 'p')
    {
        id = id.substr(id.find('|') - 2, id.size() - 1); // On parse l'id pour enlever les premières valeurs et garder qu'à partir du "sp|..."
    }
    else
    {
        id = id.substr(id.find('|') + 1, id.size() - 1); // Si | se trouve devant sp, on parse l'id pour garder qu'à partir du "sp|..."
    }
    return id;
}

/**
 * Cette fonction retourne l'acide aminé correspondant à un certain caractère défini par le fichier .psq.
 * @param c Le caractère à décoder
 * @return L'acide aminé correspondant
 */
unsigned char Database::decode_char(unsigned char c)
{
    const char encode_sequence[29] = "-ABCDEFGHIKLMNPQRSTVWXYZU*OJ";
    return encode_sequence[c];
}

/**
 * Cette fonction retourne la séquence d'une protéine se trouvant dans le fichier .psq pour un sequence_offset donné.
 * @param sequence_offset Le sequence offset de la protéine
 * @return La séquence de la protéine
 */
string Database::getProteinSequence(int sequence_offset)
{
    string protein_sequence_file_path = file_path + ".psq";
    ifstream sequence_stream(protein_sequence_file_path, ios::binary);
    string sequence;
    if (sequence_stream.is_open())
    {
        unsigned char buffer;
        sequence_stream.seekg(sequence_offset); // On se déplace directement à la position pour le sequence_offset donné
        while (sequence_stream.read((char *)&buffer, sizeof(unsigned char)))
        {
            buffer = decode_char(buffer);
            if (buffer == '-') // Un - signifie qu'on se trouve à la fin de la séquence -> on break
            {
                break;
            }
            sequence.push_back(buffer);
        }
        sequence_stream.close();
    }
    else
    {
        cerr << "Couldn't open the .psq file" << endl;
        exit(1);
    }
    return sequence;
}

/**
 * Cette fonction aligne la protéine de requête avec toutes les protéines se trouvant dans la base de données.
 * Un score est caculée pour chaque protéine en fonction de la matrice de substitution BLOSUM, du gap open penalty et du gap extension penalty.
 * Chaque protéine est enfin stockée dans un vecteur (protein_vector).
 * @param protein_query La protéine de requête
 * @param blosum La matrice BLOSUM
 * @param gap_open_penalty Le gap open penalty
 * @param gap_extension_penalty Le gap extension penalty
 */
void Database::alignProteins(Protein protein_query, vector<vector<int>> &blosum, int gap_open_penalty, int gap_extension_penalty)
{
    searchNumberProteins();
    thread threads[NUMBER_PROTEIN];
    for (int index = 0; index < NUMBER_PROTEIN; index++)
    {
        int header_offset = getOffset(index, "header");
        int sequence_offset = getOffset(index, "sequence");
        string id = getProteinId(header_offset);
        string sequence = getProteinSequence(sequence_offset);
        Protein protein = Protein(id, sequence);
        protein_vector.push_back(protein);
    }

    for (int i = 0; i < NUMBER_PROTEIN; i++)
    {
        threads[i] = thread(align_thread, ref(protein_vector[i]), ref(protein_query), ref(blosum), gap_open_penalty, gap_extension_penalty);
    }

    for (int i = 0; i < NUMBER_PROTEIN ; i++) {
        if (threads[i].joinable()) {
            threads[i].join();
        }
    }
}

void Database::align_thread(Protein &protein, Protein &protein_query, vector<vector<int>> &blosum, int gap_open_penalty, int gap_extension_penalty) {
    protein.calculate_score(protein_query.getSequence(), blosum, gap_open_penalty, gap_extension_penalty);
}

/**
 * Cette fonction trie tous les protéines pour déterminer les 20 plus alignés par rapport à la protéine de requête.
 * Le tri se fait de manière à prendre les 20 protéines ayant le plus gros score.
 * @return Un vecteur contenant les 20 protéines les plus alignés par rapport à la protéine de requête
 */
vector<Protein> Database::getMostAlignedProteins()
{
    stable_sort(protein_vector.begin(), protein_vector.end(), [](Protein a, Protein b) {return a.getScore() > b.getScore();});
    vector<Protein> most_aligned_proteins;
    most_aligned_proteins.assign(protein_vector.begin(), protein_vector.begin() + 20);
    return most_aligned_proteins;
}