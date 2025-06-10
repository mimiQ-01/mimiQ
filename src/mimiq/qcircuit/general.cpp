#include "../qcircuit.hpp"
#include "../gates.hpp"
#include <iostream>
#include <cmath>
#include <array>
#include <algorithm>

void Qcircuit::print_vector_() const {
    std::cout << "order_: \n";
    for (const auto& row : (*order_)) {
        for (int element : row) std::cout << element << " ";
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void Qcircuit::initialize(MimiqHandler* h, int nQ, int nC, std::string n) {
    name = n;
    pdf_done_ = false;
    n_qbits = nQ;
    n_cbits = nC;
    c_reg = 0;
    handler = h;
    
    std::vector<struct state_vector> sv;
    sv.resize(nQ);
    for (int i = 0; i < nQ; i++) {
        sv[i].n_qbits = 1;
        sv[i].coeffs.push_back({1, 0});
        sv[i].coeffs.push_back({0, 0});
    }
    if (nQ > 1)
        for (int i = 0; i < nQ - 1; i++) sv[i + 1] = sv[i] * sv[i + 1];
    
    state_ = new state_vector;
    (*state_) = sv[nQ - 1];
    (*state_).handler = handler;
    if(handler->qasm_gen) h->oqsm += ("qreg q[" + std::to_string(nQ) + "];\n");
    if(handler->qasm_gen) h->oqsm += ("creg c[" + std::to_string(nC) + "];\n");

    order_ = new std::vector<std::vector<int>>;
    euler_container_ = new std::vector<std::array<double, 3>>;
}
Qcircuit::Qcircuit(MimiqHandler* h, int nQ, std::string name = "") {
    initialize(h, nQ, 0, name);
}
Qcircuit::Qcircuit(MimiqHandler* h, int nQ, int nC, std::string name = "") {
    initialize(h, nQ, nC, name);
}

void Qcircuit::print_creg() const 
{ std::cout << "creg: " << c_reg << std::endl; }


void Qcircuit::setCbit1(int cbit) {
    c_reg = c_reg | (1 << (n_cbits - 1 - cbit));
}
void Qcircuit::setCbit0(int cbit) {
    c_reg = c_reg & (~(1 << (n_cbits - 1 - cbit)));
}


bool Qcircuit::accessCreg(int cbit) {
    return (c_reg & (1 << (n_cbits - 1 - cbit)));
}

// std::vector<struct experiment> simulation()
Experiment Qcircuit::simulation() {
    if (handler->circuit_drawn == false) {
        draw_circuit_();
        handler->circuit_drawn = true;
    }
    handler->qasm_gen=false;
    Experiment exp1;
    exp1.c_reg_value = c_reg;
    exp1.final_state = state_;
    exp1.n_cbits = static_cast<uint64_t>(n_cbits);
    if (!handler->latex_writer.is_open())
        std::cerr << "already closed before returning shot";
    // handler->latex_writer = NULL;
    return exp1;
}
