#pragma once 

#include "mimiqhandler.hpp"
#include "state_vector.hpp"
#include <vector>
#include <stdint.h>

enum Operation 
{
    singleGateOPeration = 1 ,
    controlledGateOperation = 2,
    twoControlledGateOperation = 3, 
    measureOperation = 4,
};

struct Qcircuit
{
private:
    std::vector<std::vector<int> > *order_; // for drawing circuits
    std::vector<std::array<double, 3> > *euler_container_; // for drawing circuits
    int n_qbits, n_cbits;
    bool pdf_done_, describe_ = false; // TODO deal with describe_
    uint64_t c_reg;
    struct state_vector* state_;

    void print_vector_() const;
    void draw_circuit_(int num=0);
    std::pair<std::string, int> draw_circuit_utility_(int key, int eulerIndex );
    void apply_utility_(const std::string &gate, int& G, int& esize,std::vector<double> arr , std::vector<std::vector<Coeff> >& matrix);
    void apply_gate_ (const std::string &gate, int t_qbit, std::vector<double> arr= std::vector<double>(), bool ORDERinc = true); // overloading apply_gate_

public:

    std::string name; 
    MimiqHandler* handler;    
    struct experiment
    {
        uint64_t c_reg_value,n_cbits ;
        struct state_vector* final_state;
    };

    void initialize(MimiqHandler* h, int nQ, int nC, std::string name );
    Qcircuit ( MimiqHandler* h, int nQ, int nC, std::string );
    Qcircuit ( MimiqHandler* h, int nQ, std::string );
    ~Qcircuit()
    {
       // delete handler; TODO analyse
    }
    void apply_controlled_gate (const std::string &gate, const std::string c_qbit_str,int t_qbit, std::vector<double> arr = std::vector<double>() );
    void u( int t_qbit, std::vector<double> arr= std::vector<double>(), bool ORDERinc = true);
    void x( int t_qbit, std::vector<double> arr= std::vector<double>(), bool ORDERinc = true);
    void y( int t_qbit, std::vector<double> arr= std::vector<double>(), bool ORDERinc = true);
    void z( int t_qbit, std::vector<double> arr= std::vector<double>(), bool ORDERinc = true);
    void h( int t_qbit, std::vector<double> arr= std::vector<double>(), bool ORDERinc = true);
    void u2( int t_qbit, std::vector<double> arr= std::vector<double>(), bool ORDERinc = true);
    void u1( int t_qbit, std::vector<double> arr= std::vector<double>(), bool ORDERinc = true);
    void t( int t_qbit, std::vector<double> arr= std::vector<double>(), bool ORDERinc = true);
    void tdg( int t_qbit, std::vector<double> arr= std::vector<double>(), bool ORDERinc = true);
    void s( int t_qbit, std::vector<double> arr= std::vector<double>(), bool ORDERinc = true);
    void sdg( int t_qbit, std::vector<double> arr= std::vector<double>(), bool ORDERinc = true);
    void ry( int t_qbit, std::vector<double> arr= std::vector<double>(), bool ORDERinc = true);

    void toffoli(int c1_qbit, int c2_qbit, int t_qbit);
    void ccx(int c1_qbit, int c2_qbit, int t_qbit);

    void cx(int cbit, int qbit);
    void ch(int cbit, int qbit);
    
    void print_creg() const;
    void measure (int tbit, int cbit, // cbit = classical // TODO: int->false
             int toPrint =0,
             int basis=0 ); // 0 basis  = along 0 or 1 , 1 bais = along + or - , 2 basis = along +i or -i 
    void setCbit1 (int cbit);
    void setCbit0 (int cbit);
    
    bool accessCreg(int cbit);

    experiment simulation();

};
    typedef Qcircuit::experiment Experiment;
