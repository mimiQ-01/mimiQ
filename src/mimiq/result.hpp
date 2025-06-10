#pragma once
#include "mimiqhandler.hpp"
#include "qcircuit.hpp"
#include <map>
#include <vector>
#include <cstdint>
#include <functional>

#define R3(x) (std::round((x) * 1000.0) / 1000.0)
struct result 
{
    uint64_t n_cbits;
    std::map < int, int > m; // TODO rename, restructure 
    std::vector<bool> creg;
    struct state_vector* state_;
    MimiqHandler* handler; 

    void generate_openqasm();
    void print_counts() const;
    std::pair<int,int> get_counts_of(int cbit) ;
};

struct result simulate(MimiqHandler* handler, std::function<Qcircuit::experiment(MimiqHandler* )> func = nullptr , int shots = 1 );
