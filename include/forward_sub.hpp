/*
Given  L = [[1,               ],
            [l_2, 1,          ],
            [    l_3, 1,      ],
            ...
            [           l_n, 1]],

and    b = [b_1, b_2, ..., b_n]^T,

Since we know L is a lower-diagonal matrix with 1's on the main diagonal,
we can just store all of the information in a vector l.

We will evaluate the linear system: Ly = b
because we will use the values of y to solve Ux = y.

Forward substitution on an upper diagonal matrix of this
specific form is an O(n) operation.
*/

#pragma once

#include <vector>

/// @brief Solves a system of linear equations using forward substitution.
/// @param lower The lower triangular matrix (flattened into a 1D vector).
/// @param b The right-hand side constant vector.
/// @return A vector containing the solution to the system.
std::vector<float> forward_substitution(const std::vector<float>& lower, const std::vector<float>& b);
