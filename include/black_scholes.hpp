#pragma once

#include <vector>
#include <string>
#include <map>
#include <fstream> 
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <iomanip>

struct GridParams {
    double price_ceiling;
    double time_to_maturity;
    size_t num_price_steps; // determines δS
    size_t num_time_steps; // determines δT 

};  

struct MarketParams {
    double volatility;
    double risk_free_interest;
    double strike_price;
};

struct Coefficients {
    double alpha;
    double beta;
    double gamma;
};

/// @brief Loads key-value pairs from a CSV file into a map.
/// @details Parses the CSV file assuming the first row contains column headers
///          and the second row contains corresponding values. Each header becomes
///          a key in the returned map, mapped to its value from the second row.
/// @param file Input file stream to read the CSV data from.
/// @param print If true, prints the parsed CSV data to standard output.
/// @return A map containing header-value pairs from the CSV file.
/// @throws std::runtime_error If the CSV has fewer than 2 rows or if the data row
///         has fewer columns than the header row.
std::map<std::string, std::string> data_loader(std::ifstream &file, bool print=false);

std::map<std::string, std::string> open_file(const std::string& filename);

double price_option(
    std::vector<double>& V, 
    const size_t& num_price_steps,
    const std::map<std::string, std::string>& params
);

void print_option_diff(const double& option_price, const std::map<std::string, std::string>& params);


Coefficients calculate_coeffs(
    const double& vol, 
    const double& r, 
    const double& time_to_maturity, 
    const size_t& time_steps, 
    const size_t& i
);

std::vector<double> evaluate_rhs(
    const std::vector<double>& V_known,
    const std::vector<double>& alpha,   
    const std::vector<double>& beta,
    const std::vector<double>& gamma,
    double V_bound_lower_j, double V_bound_lower_j_plus_1,
    double V_bound_upper_j, double V_bound_upper_j_plus_1);

std::vector<double> formulate_black_scholes(const GridParams& grid, const MarketParams& market);

std::vector<double> evaluate_system(
    const size_t& num_time_steps, 
    const size_t& num_price_steps, 
    const std::map<std::string, std::string>& params
);