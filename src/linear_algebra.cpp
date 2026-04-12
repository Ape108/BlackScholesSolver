#include "linear_algebra.hpp"

Decomposed lu_decomposition(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c) {
    
    Decomposed LU;
    LU.upper.push_back(a[0]);
    size_t n = a.size();

    for (size_t i=1; i<n; i++) {
        double l_i = c[i-1] / LU.upper[i-1];
        LU.lower.push_back(l_i);
        double u_i = a[i] - (l_i * b[i-1]);
        LU.upper.push_back(u_i); 
    }

    return LU;
}

std::vector<double> forward_substitution(const std::vector<double>& lower, const std::vector<double>& b) {
    
    std::vector<double> y = {b[0]};
    size_t n = lower.size();
    
    for (size_t i=1; i<=n; i++) {
        double y_i = b[i] - (lower[i-1] * y[i-1]);
        y.push_back(y_i);
    }

    return y;
}

std::vector<double> backward_substitution(const std::vector<double>& u, const std::vector<double>& b, const std::vector<double>& y) {

    size_t n = u.size();
    std::stack<double> x;

    x.push(y[n-1] / u[n-1]);
    
    for (size_t i = n-1; i > 0; i--) {
        double x_i = y[i-1] - (b[i-1] * x.top());
        x_i /= u[i-1];
        x.push(x_i);
    }

    std::vector<double> solution;

    // Pop all values from stack into the vector in correct order
    while (!x.empty()) {
        double top = x.top();
        solution.push_back(top);
        x.pop();
    }

    return solution;
}