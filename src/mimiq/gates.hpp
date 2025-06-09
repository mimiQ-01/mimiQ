
#pragma once

#include <vector>

#include "state_vector.hpp"

#define R3(x) (std::round((x) * 1000.0) / 1000.0)

namespace GATES 
{
    std::vector < std::vector < Coeff >> unitaryGate(double theta, double psi, double lam);
    void controlled_pauli_X(struct state_vector & sv, int cbit, int tbit);
    void toffoli(struct state_vector & sv, int cbit1, int cbit2, int tbit);
    void controlled_hadamard(struct state_vector & sv, int cbit, int tbit);
    void matrixmultiply2x2(const std::vector<std::vector<Coeff>>& g1, const std::vector<std::vector<Coeff>>& g2, std::vector<std::vector<Coeff>>& g3) ;
    void controlled_pauli_Z(struct state_vector & sv, int cbit, int tbit) ;
    void printMatrix(const std::vector < std::vector < Coeff >> & matrix) ;
    void kroneckerSUB(const std::vector < std::vector < Coeff >> & m1,
    const std::vector < std::vector < Coeff >> & m2, std::vector < std::vector < Coeff >> & m3);
    void kroneckerBUFF(const std::vector < std::vector < Coeff >> & gate, std::vector < std::vector < Coeff >> & res, int n_qbits, int tbit) ;
    std::vector < std::vector < Coeff >> inverse2x2(const std::vector < std::vector < Coeff >> & matrix);

    enum gateNumber
    {
        h = 1, 
        x = 2,
        y = 3,
        z = 4,
        s = 5,
        sd = 6,
        t = 7,
        td = 8,
        rx = 9,
        ry = 10,
        u3 = 11, 
        iu3 = 12,
        u1 = 13,
        u2 = 14,
    };
};
