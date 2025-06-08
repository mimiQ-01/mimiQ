
#ifndef _mimiq
#define _mimiq

#include <vector>
#include <fstream>
#include <cstring>
#include <map>
#include <functional>
#include <array>
#include <cstdint>
#include <utility>
#include <cmath>
#include <iostream>


struct mimiqHandler 
{
    bool circuittDrawn, reportGenerated, qasmgen , canqasm ;
    std::ofstream wr; // to write Latex for pdf output
    std::string oqsm, // to write openQASM 2.0 code
    dir_path ; // where pdf should be stored

    mimiqHandler(std::string path = "");
    void clean();
    void writeInPdf(std::string msg);
    void generateReport();


};
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



// a qbit as from our understanding can be represented as a sphere will 3 angles which correspond to the rotations from reference axes x,y,z
// in this project, a qbit is represented in the form of a|0> + b|1> 
// the reason is for its simplicity and ease to manipulate it mathematically
// a is called the real part, and b is called the complex part (b is imaginary part where a,b together is called complex number, but bare with the naming sense lol)



static std::vector<std::vector<Coeff> > empty_matrix ; // for applyUtility
// our prime basis states as we told are {|0>, |1>}
// so the number of basis states from n qbits is 2^n

// the state_vector is the array of the Coefficients of all the basis states from a system of n_qbits - number of qbits
struct state_vector;

struct result 
{
    uint64_t n_cbits;
    std::map < int, int > m; // TODO rename, restructure 
    std::vector<bool> creg;
    struct state_vector* state;
    mimiqHandler* handler; 
    
    void generate_openqasm();
    void print_counts();
    std::pair<int,int> get_counts_of(int cbit);
    
};


struct Qcircuit
{
private:
    std::vector<std::vector<int> > *ORDER; // for drawing circuits
    std::vector<std::array<double, 3> > *euler_container; // for drawing circuits
    int n_qbits, n_cbits;
    bool pdfDone, describe = false; // TODO deal with describe
    uint64_t c_reg;
    struct state_vector* state;
    void printVector();
    void drawCircuit(int num);
    std::pair<std::string, int> drawCircuitUtility(int key, int eulerIndex );
    void applyUtility(const std::string &gate, int& G, int& esize,std::vector<double> arr , std::vector<std::vector<Coeff> >& matrix);
    void applyGate (const std::string &gate, int t_qbit, std::vector<double> arr= std::vector<double>(), bool ORDERinc = true); // overloading applyGate

public:

    std::string name; 
    mimiqHandler* handler;    
    struct experiment
    {
        uint64_t c_reg_value,n_cbits ;
        struct state_vector* final_state;
    };

    void initialize(mimiqHandler* h, int nQ, int nC, std::string name );
    Qcircuit ( mimiqHandler* h, int nQ, int nC, std::string );
    Qcircuit ( mimiqHandler* h, int nQ, std::string );
    ~Qcircuit()
    {
       // delete handler; TODO analyse
    }
    void applyControlledGate (const std::string &gate, const std::string c_qbit_str,int t_qbit, std::vector<double> arr = std::vector<double>() );


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
    
    void print_creg ();
    void measure (int tbit, int cbit, // cbit = classical // TODO: int->false
             int toPrint =0,
             int basis=0 ); // 0 basis  = along 0 or 1 , 1 bais = along + or - , 2 basis = along +i or -i 
    void setCbit1 (int cbit);
    void setCbit0 (int cbit);
    
    bool accessCreg(int cbit);

    experiment simulation();


};
    typedef Qcircuit::experiment Experiment;


struct result simulate(mimiqHandler* handler, std::function<Qcircuit::experiment(mimiqHandler* )> func = nullptr , int shots = 1 );

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
}


#endif