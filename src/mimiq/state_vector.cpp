#include "state_vector.hpp"

#include <iostream>
#include <cmath>

// overloading operators for simplifying arithmetic operations with Coeff type
// structures
Coeff Coeff::operator*(const Coeff& other) const {
    Coeff result;
    result.real = real * other.real - complex * other.complex;
    result.complex = real * other.complex + complex * other.real;
    return result;
}
Coeff Coeff::operator+(const Coeff& other) const {
    Coeff result;
    result.real = real + other.real;
    result.complex = complex + other.complex;
    return result;
}
Coeff Coeff::operator-(const Coeff& other) const {
    Coeff result;
    result.real = real - other.real;
    result.complex = complex - other.complex;
    return result;
}
Coeff Coeff::operator/(const Coeff& other) const {
    Coeff result;
    result.real = ((real * other.real) + (complex * other.complex)) /
                  ((other.real * other.real) + (other.complex * other.complex));
    result.complex =
        ((other.complex * real) - (real * other.complex)) /
        ((other.real * other.real) + (other.complex * other.complex));
    return result;
}

// will return the square of amplitude though the name of the function is
// amplitude.
double Coeff::amp_sq() const {  // amp squared actually
    return (real * real) + (complex * complex);
}

// our prime basis states as we told are {|0>, |1>}
// so the number of basis states from n qbits is 2^n

// the state_vector is the array of the Coefficients of all the basis states
// from a system of n_qbits - number of qbits

state_vector::state_vector(int n) {
    n_qbits = n;
    coeffs.resize(1 << n, {0, 0});
}

state_vector::state_vector(const struct state_vector& sv) {
    n_qbits = sv.n_qbits;
    coeffs = sv.coeffs;
}

// appropriate overloading for ease in arithmetic operations
struct state_vector& state_vector::operator=(const struct state_vector& other) {
    // Check for self-assignment
    if (this == &other) {
        return *this;
    }

    // Copy the number of qubits
    n_qbits = other.n_qbits;

    // Deep copy the vector of coefficients
    coeffs = other.coeffs;

    return *this;
}
struct state_vector state_vector::operator*(
    const struct state_vector&
        other)  // returns a tensor product of 2 state_ vectors
{
    struct state_vector sv3;  //(*this).n_qbits == n_qbits
    sv3.n_qbits = n_qbits + other.n_qbits;
    int cnd1 = 1 << n_qbits, cnd2 = 1 << other.n_qbits;
    for (int i = 0; i < cnd1; i++)
        for (int j = 0; j < cnd2; j++)
            sv3.coeffs.push_back(coeffs[i] * other.coeffs[j]);
    return sv3;
}
void state_vector::print() {
    handler->latex_writer << "\\text{The state_ vector for the last shot is as follows: }";
    handler->latex_writer << "\\[\n\\begin{array}{@{}llll@{}}\n";
    std::cout << std::endl << "number of qbits = " << n_qbits << std::endl;
    int cnd = 1 << n_qbits;
    for (int i = 0; i < cnd; i++) {
        std::string binaryString;
        for (int j = n_qbits - 1; j >= 0; --j) {
            int bit = (i >> j) & 1;
            binaryString += (bit == 0) ? '0' : '1';
        }
        handler->latex_writer << "\\text{" << binaryString << ":} & "
                    << coeffs[i].amp_sq() * 100 << "\\% & " << coeffs[i].real
                    << " |0\\rangle &  ";
        handler->latex_writer << coeffs[i].complex << " |1\\rangle \\\\" << std::endl;
        std::cout << binaryString << " =>" << coeffs[i].amp_sq() * 100 << "% ( "
                  << coeffs[i].real << " |0> + " << coeffs[i].complex
                  << " i |1> )\n";
    }
    handler->latex_writer << "\\end{array}\n\\]\n";
    std::cout << std::endl;
}

std::pair<Coeff, Coeff> state_vector::measure_along(int bit) {
    int cnd = 1 << n_qbits;
    Coeff z0, z1;
    auto tmask = 1 << (n_qbits - 1 - bit);
    // std::cout << "mask: "<<tmask<<std::endl;
    for (int i = 0; i < cnd; i++) {
        // std::cout << "i: "<< i << " tmask: "<< tmask << "__ i | tmask: " <<
        // (i | tmask ) << "__ i & tmask: "<< (i & tmask) << "__i &~ tmask: "<<
        // (i & (~tmask)) <<std::endl;
        if (i & tmask) {
            z1.real += (coeffs[i].real * coeffs[i].real);
            z1.complex += (coeffs[i].complex * coeffs[i].complex);
        } else {
            z0.real += (coeffs[i].real * coeffs[i].real);
            z0.complex += (coeffs[i].complex * coeffs[i].complex);
        }
    }
    z0.real = std::sqrt(z0.real);
    z0.complex = std::sqrt(z0.complex);
    z1.real = std::sqrt(z1.real);
    z1.complex = std::sqrt(z1.complex);
    return {z0, z1};
}
void state_vector::printprobs() {
    handler->latex_writer << "\\begin{align*}\n";
    int cnd = 1 << n_qbits;
    for (int bit = 0; bit < n_qbits; bit++) {
        auto tmask = 1 << (n_qbits - 1 - bit);
        double prob = 0.0;
        for (int i = 0; i < cnd; i++)
            if (i & tmask) prob += coeffs[i].amp_sq();
        std::cout << "for qbit " << bit << ": " << "Prob of |1> : " << prob
                  << ",  ";
        std::cout << "Prob of |0> : " << 1 - prob << std::endl;
        handler->latex_writer << "\\text{For qbit " << bit
                    << "} \\quad & \\text{Probability of} |1\\rangle: " << prob
                    << ", \\text{Probability of} |0\\rangle: " << 1 - prob
                    << " \\\\\n";
    }
    std::cout << std::endl;
    handler->latex_writer << "\\end{align*}\n";
}