#include "../include/lower_upper.hpp"
#include "../include/forward_sub.hpp"
#include "../include/back_sub.hpp"

#include <iostream>


int main() {

    // Matrix A
    std::vector<float> a = {1, 7, 5};
    std::vector<float> b = {1, 8};
    std::vector<float> c = {2, 3};
    
    // Right-hand side
    std::vector<float> b = {6.0, 9.0, 6.0};

    Decomposed LU = lu_decomposition(a, b, c);

    std::vector<float> l = LU.lower;
    std::vector<float> y = forward_substitution(l, b);

    std::vector<float> u = LU.upper;
    std::vector<float> x = backward_substitution(u, b, y);

    std::cout << x[0] << " " << x[1] << " " << x[2] << std::endl;

    return 0;
}