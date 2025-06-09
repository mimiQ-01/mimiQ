#pragma once 

#include <vector>
#include "mimiqhandler.hpp"

#define R3(x) (std::round((x) * 1000.0) / 1000.0)
struct Coeff 
{
    // prime basis state: |0> , |1>
    // for other basis sets: {|+>, |->}, {|i>, |-i>}, spontaneous computation is done as they are not often used
    double real;
    double complex;

    // default constructors
    Coeff(): real(0), complex(0) {}
    Coeff(double r, double c): real(r), complex(c) {}

    // overloading operators for simplifying arithmetic operations with Coeff type structures
    Coeff operator * (const Coeff & other) const ;
    Coeff operator + (const Coeff & other) const ;
    Coeff operator - (const Coeff & other) const ;
    Coeff operator / (const Coeff & other) const ;
    // will return the square of amplitude though the name of the function is amplitude. 
    double amp_sq() const ;
};
struct state_vector {
    int n_qbits; 
    std::vector < Coeff > coeffs; // 16 * (1<< nq) // for 2 qbits -> 544 cbits, for 3qbits -> 1056 cbits
    mimiqHandler* handler; 
    // default constuctor
    state_vector(): n_qbits(0) {}
    state_vector(int n);
    
    state_vector(const struct state_vector & sv);

    // appropriate overloading for ease in arithmetic operations
    struct state_vector& operator=(const struct state_vector& other) ;
    struct state_vector operator * (const struct state_vector & other);
    void print() ;
    std::pair<Coeff,Coeff> measureAlong(int bit);
    void printprobs();
};

