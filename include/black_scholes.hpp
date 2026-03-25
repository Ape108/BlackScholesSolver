#pragma once

#include <vector>
#include <string>
#include <map>
#include <fstream> 
#include <sstream>
#include <iostream>

// I'll write a data loader to grab the inputs from the CSV 
// Then I will discretize the grid,
// and initialize the boundary conditions.

/// @brief Loads key-value pairs from a CSV file into a map.
/// @details Parses the CSV file assuming the first row contains column headers
///          and the second row contains corresponding values. Each header becomes
///          a key in the returned map, mapped to its value from the second row.
/// @param file Input file stream to read the CSV data from.
/// @param print If true, prints the parsed CSV data to standard output.
/// @return A map containing header-value pairs from the CSV file.
/// @throws std::runtime_error If the CSV has fewer than 2 rows or if the data row
///         has fewer columns than the header row.
std::map<std::string, std::string> data_loader(std::ifstream &file, const bool& print=false);