#pragma once
#include <string>
#include <vector>

#include "protein.hpp"

class Query
{
private:
    std::string file_path;

public:
    Query(std::string file_path);
    Protein makeProteinFromQuery();
};