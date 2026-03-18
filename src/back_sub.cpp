#include "../include/back_sub.hpp"

std::vector<float> backward_substitution(const std::vector<float>& u, const std::vector<float>& b, const std::vector<float>& y) {

    std::size_t n = u.size();
    std::vector<float> x = {static_cast<float>(y[n]) / u[n]};
    
    for (std::size_t i=n-1; i >= 0; i--) {
        float x_i = y[i] - (b[i+1] * y[i+1]);
        // x.push_front() DONT PUSH FRONT because it's O(n) to shift everything. 
        // Ill use a stack then pop everything from the stack into the vector result
        // at the end.
    }

    return x;
}