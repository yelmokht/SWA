#ifndef _QUERY_HPP
#define _QUERY_HPP
#include <string>
#include <vector>

#include "protein.hpp"

class Query
{
private:
    std::string file_path;

public:
    Query(std::string file_path);
    std::vector<std::vector<int>> makeMatrixBlosum();
    Protein makeProteinFromQuery();
};
#endif