/**
 * @file mimiq.cpp
 * @brief implementation of mimiq.h functions
 * 
 *  
 * @author Rushikesh Muraharisetty
 * @date last updated: Mar 10 '25
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <ctime>
#include <functional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "mimiq.h"

#define root2 std::sqrt(2)
#define PI acos(0.0) * 2
#define R3(x) (std::round((x) * 1000.0) / 1000.0)



mimiqHandler::mimiqHandler(std::string path) {
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
    dir_path = path;
    wr.open(dir_path + "report.tex");
wr << "\\documentclass{article}\n";
wr << "\\usepackage[margin=1in]{geometry}\n";
wr << "\\usepackage{datetime}\n";
wr << "\\usepackage{qcircuit}\n";
wr << "\\usepackage{amsmath}\n";
wr << "\\usepackage{pgfplots}\n";
wr << "\\pgfplotsset{compat=1.18}\n";
wr << "\\usepackage[utf8]{inputenc}\n";
wr << "\\begin{document}\n";
wr << "\\begin{center}\n";
wr << "    {\\LARGE \\textbf{mimiQ++ Report}} \\\\\n";
wr << "    \\large \\today \\quad \\currenttime\n";
wr << "\\end{center}\n";
wr << "\\hrule\n";
wr << "\\vspace{1cm}\n";

    circuittDrawn = false;
    qasmgen = true;
    canqasm= true;
    oqsm += "OPENQASM 2.0;\ninclude \"qelib1.inc\";\n";

}
void mimiqHandler::clean() { circuittDrawn = false; qasmgen = true; canqasm = true; oqsm = "";}

void mimiqHandler::writeInPdf(std::string msg) {
    if (wr.is_open() && wr.good()) wr << msg << "\n\n";
}

void mimiqHandler::generateReport() {
    if (wr.is_open()) {
        reportGenerated = true;
        wr << "\\end{document}";
        wr.close();
    } else
        std::cerr << "Unable to generate report as either already generated "
                     "earlier (or) writer closed\n";
    std::string latexFilename = dir_path + "report.tex";
    // Compile LaTeX file into PDF
    std::string compileCommand = "pdflatex -output-directory=" + dir_path +
                                 " " + latexFilename + " > nul 2>&1";
    int exitCode = std::system(compileCommand.c_str());

    if (exitCode == 0)
        std::cout << "PDF generated: "
                  << latexFilename.substr(0, latexFilename.find_last_of('.'))
                  << ".pdf" << std::endl;
    else
        std::cerr << "Error occurred during latex compilation." << std::endl;

    std::string tmp2 = dir_path + "report.log", tmp3 = dir_path + "report.aux";
    std::remove(tmp2.c_str());
    std::remove(tmp3.c_str());
}

// a qbit as from our understanding can be represented as a sphere will 3 angles
// which correspond to the rotations from reference axes x,y,z in this project,
// a qbit is represented in the form of a|0> + b|1> the reason is for its
// simplicity and ease to manipulate it mathematically a is called the real
// part, and b is called the complex part (b is imaginary part where a,b
// together is called complex number, but bare with the naming sense lol)

// prime basis state: |0> , |1>
// for other basis sets: {|+>, |->}, {|i>, |-i>}, spontaneous computation is
// done as they are not often used

// default constructors

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
        other)  // returns a tensor product of 2 state vectors
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
    handler->wr << "\\text{The state vector for the last shot is as follows: }";
    handler->wr << "\\[\n\\begin{array}{@{}llll@{}}\n";
    std::cout << std::endl << "number of qbits = " << n_qbits << std::endl;
    int cnd = 1 << n_qbits;
    for (int i = 0; i < cnd; i++) {
        std::string binaryString;
        for (int j = n_qbits - 1; j >= 0; --j) {
            int bit = (i >> j) & 1;
            binaryString += (bit == 0) ? '0' : '1';
        }
        handler->wr << "\\text{" << binaryString << ":} & "
                    << coeffs[i].amp_sq() * 100 << "\\% & " << coeffs[i].real
                    << " |0\\rangle &  ";
        handler->wr << coeffs[i].complex << " |1\\rangle \\\\" << std::endl;
        std::cout << binaryString << " =>" << coeffs[i].amp_sq() * 100 << "% ( "
                  << coeffs[i].real << " |0> + " << coeffs[i].complex
                  << " i |1> )\n";
    }
    handler->wr << "\\end{array}\n\\]\n";
    std::cout << std::endl;
}

std::pair<Coeff, Coeff> state_vector::measureAlong(int bit) {
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
    handler->wr << "\\begin{align*}\n";
    int cnd = 1 << n_qbits;
    for (int bit = 0; bit < n_qbits; bit++) {
        auto tmask = 1 << (n_qbits - 1 - bit);
        double prob = 0.0;
        for (int i = 0; i < cnd; i++)
            if (i & tmask) prob += coeffs[i].amp_sq();
        std::cout << "for qbit " << bit << ": " << "Prob of |1> : " << prob
                  << ",  ";
        std::cout << "Prob of |0> : " << 1 - prob << std::endl;
        handler->wr << "\\text{For qbit " << bit
                    << "} \\quad & \\text{Probability of} |1\\rangle: " << prob
                    << ", \\text{Probability of} |0\\rangle: " << 1 - prob
                    << " \\\\\n";
    }
    std::cout << std::endl;
    handler->wr << "\\end{align*}\n";
}

std::string reverse(const std::string& str) {
    std::string reversed = str;  // Create a copy of the input string
    std::reverse(reversed.begin(), reversed.end());  // Reverse the string
    return reversed;
}
std::string generateLatexBarChart(
    const std::vector<std::pair<std::string, std::string>>& data) {
    std::ostringstream latex;
    latex << "\\begin{center}\n";
    latex << "    \\begin{tikzpicture}\n";
    latex << "        \\begin{axis}[\n";
    latex << "            ybar,\n";
    latex << "            symbolic x coords={";

    // Add x-coordinates (Quantum States)
    for (size_t i = 0; i < data.size(); ++i) {
        latex << data[i].first;
        if (i != data.size() - 1) latex << ", ";
    }
    latex << "},\n";

    latex << "            xtick=data,\n";
    latex << "            xlabel={Quantum States},\n";
    latex << "            ylabel={Counts},\n";
    latex << "            ymin=0,\n";
    latex << "            bar width=20pt,\n";
    latex << "            width=10cm,\n";
    latex << "            height=7cm,\n";
    latex << "            nodes near coords,\n";
    latex << "            nodes near coords align={vertical},\n";
    latex << "            enlarge x limits=0.3,\n";
    latex << "            title={Quantum State Counts}\n";
    latex << "        ]\n";

    // Add plot coordinates
    latex << "        \\addplot coordinates {";
    for (const auto& [state, count] : data) {
        latex << "(" << state << "," << count << ") ";
    }
    latex << "};\n";

    latex << "        \\end{axis}\n";
    latex << "    \\end{tikzpicture}\n";
    latex << "\\end{center}\n";

    return latex.str();
}

void result::print_counts() {
    std::vector<std::pair<std::string, std::string>> gr;
   /* handler->wr
        << "\n\n\n\n\n\n\\text{Classical register readings (left to right: "
           "cn,cn-1,..c2,c1,c0) for the simulation:} \n\n\n\n\n\n";*/
    std::cout << "classical register readings for the simulation: "
              << std::endl;
    for (auto i = m.begin(); i != m.end(); i++) {
        std::string binaryString;
        for (int j = n_cbits - 1; j >= 0; --j) {
            int bit = (i->first >> j) & 1;
            binaryString += (bit == 0) ? '0' : '1';
        }
        auto rev = reverse(binaryString);
        std::cout << rev << ": " << i->second << std::endl;
       // handler->wr << rev << ": " << i->second << "\n\n\n";
        gr.push_back({rev, std::to_string(i->second)});
    }
    handler->wr << generateLatexBarChart(gr);
    std::cout << std::endl;
    handler->wr << "\n\n";
}

void insertBackslashes(std::string& str) {
    for (size_t i = 0; i < str.size(); ++i) {
        if (str[i] == ';') {
            str.insert(i + 1, "\\\\");
            i += 2; // Move past the inserted "\\"
        }
    }
}

void result::generate_openqasm()
{
    std::cout<< handler->oqsm<<"\n";
    //insertBackslashes(handler->oqsm);
    handler->writeInPdf("the OpenQASM 2.0 code for the above qircuit is: \n");
    handler->wr<<"\\begin{verbatim}\n";
    handler->wr<< handler->oqsm << "\\end{verbatim}\n";
}

std::pair<int, int> result::get_counts_of(int cbit) {
    int cmask = n_cbits - 1 - cbit;
    int ones = 0, zeroes = 0;
    int cnd = 1 << state->n_qbits;

    for (int i = 0; i < cnd; i++) {
        if (i & cmask)  // x1x val = 1
            ones += m[i];
        else
            zeroes += m[i];
    }
    return {zeroes, ones};
}

struct wire {
    std::vector<std::pair<std::string, int>> line;
    int pos;
    bool nature;  // TODO
};

std::string trimZeroes(std::string str) {
    size_t decimalIndex = str.find('.');
    if (decimalIndex != std::string::npos) {
        size_t lastNonZeroIndex = std::string::npos;
        for (auto i = str.size() - 1; i > decimalIndex; --i) {
            if (str[i] != '0') {
                lastNonZeroIndex = i;
                break;
            }
        }
        if (lastNonZeroIndex != std::string::npos) {
            if (lastNonZeroIndex != decimalIndex) {
                return str.substr(0, lastNonZeroIndex + 1);
            } else {
                return str.substr(0, decimalIndex);
            }
        } else {
            return str.substr(0, decimalIndex);
        }
    } else {
        return str;
    }
}

void Qcircuit::printVector() {
    std::cout << "ORDER: \n";
    for (const auto& row : (*ORDER)) {
        for (int element : row) std::cout << element << " ";
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void Qcircuit::initialize(mimiqHandler* h, int nQ, int nC, std::string n) {
    name = n;
    pdfDone = false;
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
    
    state = new state_vector;
    (*state) = sv[nQ - 1];
    (*state).handler = handler;
    if(handler->qasmgen) h->oqsm += ("qreg q[" + std::to_string(nQ) + "];\n");
    if(handler->qasmgen) h->oqsm += ("creg c[" + std::to_string(nC) + "];\n");

    ORDER = new std::vector<std::vector<int>>;
    euler_container = new std::vector<std::array<double, 3>>;
}
Qcircuit::Qcircuit(mimiqHandler* h, int nQ, std::string name = "") {
    initialize(h, nQ, 0, name);
}
Qcircuit::Qcircuit(mimiqHandler* h, int nQ, int nC, std::string name = "") {
    initialize(h, nQ, nC, name);
}

void Qcircuit::print_creg() { std::cout << "creg: " << c_reg << std::endl; }

void Qcircuit::measure(int tbit,
                       int cbit,  // cbit = classical // TODO: int->false
                       int toPrint,
                       int basis)  // 0 basis  = along 0 or 1 , 1 bais = along +
                                   // or - , 2 basis = along +i or -i
{
    if(handler->qasmgen) handler->oqsm+="measure q["+ std::to_string(tbit) +"] -> c[" +std::to_string(cbit) +"];\n" ;
    std::pair<Coeff, Coeff> res1;
    if (basis == 1) {
        std::vector<double> placebo;
        applyGate("h", tbit, placebo, false);
        res1 = (*state).measureAlong(tbit);
        applyGate("h", tbit, placebo, false);
    } else
        res1 = (*state).measureAlong(tbit);

    if (describe || toPrint) {
        std::cout << "for 0: " << res1.first.real << " " << res1.first.complex
                  << "\n";
        std::cout << "for 1: " << res1.second.real << " " << res1.second.complex
                  << "\n";
    }
    // res = a | 0> + b | 1 > ;
    int cnd = 1 << n_qbits;                  // for 2 qbits , cnd = 4
    auto tmask = 1 << (n_qbits - 1 - tbit);  // qbit: 0 means tmask
    // 2-1-0 = 1 shit to left.

    double prob = res1.second.amp_sq();  // b*b
    double randomValue;
    for (int l = 1; l <= 3; l++)
        randomValue = static_cast<double>(rand()) / RAND_MAX;
    bool res = (randomValue <= prob);
    if (toPrint == 1 || describe) {
        std::cout << "Probability of measuring |1> "
                     "on qubit "
                  << tbit << ": " << prob << std::endl;
        std::cout << "Probability of measuring |0> "
                     "on qubit "
                  << tbit << ": " << 1 - prob << std::endl;
        std::cout << "rand/prob: " << randomValue << "/ " << prob
                  << " res: " << res << " " << (1 << (n_cbits - 1 - cbit))
                  << std::endl;
    }
    if (res)
        c_reg = c_reg | (1 << (n_cbits - 1 - cbit));
    else
        c_reg = c_reg & (~(1 << (n_cbits - 1 - cbit)));

    // now update state vector to cause collapse
    // TODO: the search is exhaustive, u only need to check limited. Find a
    // better way if possible
    Coeff normaliser = {0, 0};
    if (res)  // if res== 1 , take the zero value and add it to 1 value
    {
        if (describe) std::cout << "measured 1 so.. tmask " << tmask << " \n";
        for (int i = 0; i < cnd; ++i) {
            if (i & (tmask))  // we got x1x
            {
                if (describe)
                    std::cout << "for i =  " << i << "adding to "
                              << (i | (tmask)) << " from " << (i & (~tmask))
                              << std::endl;

                normaliser.real =
                    normaliser.real + (*state).coeffs[i | tmask].amp_sq();
                (*state).coeffs[i & (~tmask)] = {0, 0};
            }
        }
    } else  // res == 0
    {
        if (describe) std::cout << "measured 0 so.. tmask is " << tmask << "\n";
        for (int i = 0; i < cnd; ++i) {
            if (i & (tmask))  // we got x1x
            {
                if (describe)
                    std::cout << "for i =  " << i << " adding to "
                              << (i & (~tmask)) << " from " << (i | (tmask))
                              << std::endl;

                normaliser.real =
                    normaliser.real + (*state).coeffs[i & (~tmask)].amp_sq();
                (*state).coeffs[i | tmask] = {0, 0};
            }
        }
    }
    normaliser.real = std::sqrt(normaliser.real);
    if (res)  // if res== 1 , take the zero value and add it to 1 value
    {
        if (describe) std::cout << "measured 1 so.. tmask " << tmask << " \n";
        for (int i = 0; i < cnd; ++i) {
            if (i & (tmask))  // we got x1x
                (*state).coeffs[i | tmask] = (*state).coeffs[i | tmask] / normaliser;
        }
    } else  // res == 0
    {
        if (describe) std::cout << "measured 0 so.. tmask is " << tmask << "\n";
        for (int i = 0; i < cnd; ++i) {
            if (i & (tmask))  // we got x1x
                (*state).coeffs[i & (~tmask)] =
                    (*state).coeffs[i & (~tmask)] / normaliser;
        }
    }
    if (describe) std::cout << "updated creg: " << c_reg << std::endl;

    // std::cout<<"m "<< tbit << " "<<cbit << std::endl;
    (*ORDER).push_back({3, tbit, cbit, toPrint, basis});
}
void Qcircuit::setCbit1(int cbit) {
    c_reg = c_reg | (1 << (n_cbits - 1 - cbit));
}
void Qcircuit::setCbit0(int cbit) {
    c_reg = c_reg & (~(1 << (n_cbits - 1 - cbit)));
}

std::pair<std::string, int> Qcircuit::drawCircuitUtility(int key,
                                                         int eulerIndex) {
    int length = 1;
    // std::cout<<"received with key: "<< key <<" , index: "<< eulerIndex ;
    std::string res;
    auto& euler_array =
        (*euler_container)[eulerIndex];  // Explicitly reference the array
    switch (key) {
        case 1:
            res = "\\gate{H} & ";
            break;
        case 2:
            res = "\\gate{X} & ";
            break;
        case 3:
            res = "\\gate{Y} & ";
            break;
        case 4:
            res = "\\gate{Z} & ";
            break;
        case 5:
            res = "\\gate{S} & ";
            break;
        case 6:
            res = "\\gate{S^\\dag} & ";
            break;
        case 7:
            res = "\\gate{T} & ";
            break;
        case 8:
            res = "\\gate{T^\\dag} & ";
            break;
        case 101:
            res =
                "\\gate{U(" + trimZeroes(std::to_string(R3(euler_array[0]))) +
                "," +
                trimZeroes(std::to_string(R3((*euler_container)[eulerIndex][1]))) +
                "," +
                trimZeroes(std::to_string(R3((*euler_container)[eulerIndex][2]))) +
                ")} & ";
            length += 3;
            break;
        case 102:
            res =
                "\\gate{IU3(" +
                trimZeroes(std::to_string(R3((*euler_container)[eulerIndex][0]))) +
                "," +
                trimZeroes(std::to_string(R3((*euler_container)[eulerIndex][1]))) +
                "," +
                trimZeroes(std::to_string(R3((*euler_container)[eulerIndex][2]))) +
                ")} & ";
            length += 3;
            break;
        case 103:
            res =
                "\\gate{U1(" +
                trimZeroes(std::to_string(R3((*euler_container)[eulerIndex][0]))) +
                ")} & ";
            length++;
            break;
        case 104:
            res = ("\\gate{U2(" +
                   trimZeroes(
                       std::to_string(R3((*euler_container)[eulerIndex][0]))) +
                   "," +
                   trimZeroes(
                       std::to_string(R3((*euler_container)[eulerIndex][1]))) +
                   ")} & ");
            length += 2;
            break;

        case 105:
            res =
                "\\gate{Rx(" +
                trimZeroes(std::to_string(R3((*euler_container)[eulerIndex][0]))) +
                ")} & ";
            length++;
            break;
        case 106:
            res =
                "\\gate{Ry(" +
                trimZeroes(std::to_string(R3((*euler_container)[eulerIndex][0]))) +
                ")} & ";
            length++;
            break;
    }
    // std::cout<<" returns: "<<res<<std::endl;
    return {res, length};
}
void Qcircuit::drawCircuit(int num = 0) {
    /*
        this function is responsible for drawing the circuit using Qcircuit
       typesetting langauge (latex)

    */
    // printVector();
    int maxt = 1;
    struct wire qw[n_qbits], cw[n_cbits];
    // init the lines
    for (int i = 0; i < n_qbits; i++) {
        qw[i].line.push_back({("\\lstick{q" + std::to_string(i) + "} & "), 1});
        qw[i].pos = 1;
    }
    for (int i = 0; i < n_cbits; i++) {
        cw[i].line.push_back({("\\lstick{c" + std::to_string(i) + "} & "), 1});
        cw[i].pos = 1;
    }

    // int tmp1, tmp2; // utilities
    for (std::vector<std::vector<int>>::size_type i = 0; i < (*ORDER).size();
         i++) {
        // std::cout<<"current: ";
        // for(auto q : (*ORDER)[i]) std::cout<< q <<" ";
        // std::cout<<std::endl;
        if ((*ORDER)[i][0] == 1) {
            qw[(*ORDER)[i][1]].pos++;
            if (qw[(*ORDER)[i][1]].pos > maxt) {
                maxt = qw[(*ORDER)[i][1]].pos;
            }
            auto tmp = drawCircuitUtility((*ORDER)[i][2], (*ORDER)[i][3]);
            qw[(*ORDER)[i][1]].line.push_back({tmp.first, tmp.second});

        }
        // TODO OPTIMISE DRAWCIRCUIT
        else if ((*ORDER)[i][0] == 4) {
            // (*ORDER).push_back ({ 4, t_qbit, G, c_or_q1, control_bit1,
            // control_bit2, esize});
            int tbit = (*ORDER)[i][1];
            int gate = (*ORDER)[i][2];
            int cq = (*ORDER)[i][3];
            int c1 = (*ORDER)[i][4];
            int c2 = (*ORDER)[i][5];
            int eindex = (*ORDER)[i][6];
            /*
                cases:
                1. c c t
                2. c t c
                3. t c c

                top: min of all 3
                bottom: max of all 3
                mid: the middle one

                maxt++
                from top+1 to bottom-1 excluding mid, fill till maxt
                for top,mid,bottom, fill till maxt-1

                add cntl, target for them

            */
            int top = std::min({tbit, c1, c2});
            int bottom = std::max({tbit, c1, c2});
            int mid = tbit + c1 + c2 - top - bottom;

            for (int i = top + 1; i < bottom; i++) {
                if (i == mid) continue;
                while (qw[i].pos <= maxt) {
                    qw[i].line.push_back(
                        {"\\qw & ",
                         1});  // Append "\\qw & " with a length increment of 1
                    qw[i].pos++;
                }
            }
            while (qw[mid].pos < maxt) {
                qw[mid].line.push_back(
                    {"\\qw & ",
                     1});  // Append "\\qw & " with a length increment of 1
                qw[mid].pos++;
            }
            while (qw[top].pos < maxt) {
                qw[top].line.push_back(
                    {"\\qw & ",
                     1});  // Append "\\qw & " with a length increment of 1
                qw[top].pos++;
            }
            while (qw[bottom].pos < maxt) {
                qw[bottom].line.push_back(
                    {"\\qw & ",
                     1});  // Append "\\qw & " with a length increment of 1
                qw[bottom].pos++;
            }

            qw[c1].pos++;
            qw[c1].line.push_back(
                {"\\ctrl{" + std::to_string(tbit - c1) + "} & ",
                 1});  // Append control gate

            qw[c2].pos++;
            qw[c2].line.push_back(
                {"\\ctrl{" + std::to_string(tbit - c2) + "} & ",
                 1});  // Append control gate

            qw[tbit].pos++;
            auto tmp = drawCircuitUtility(gate, 0);
            qw[tbit].line.push_back(
                {tmp.first,
                 tmp.second});  // Append the result of drawCircuitUtility
            maxt++;

        } else if ((*ORDER)[i][0] == 2) {
            //  (*ORDER).push_back ({ 2, t_qbit, GATE, c_or_q, control_bit});
            if ((*ORDER)[i][3] == 1)  // type: q
            {
                // std::cout<<"control: G:  "<<(*ORDER)[i][2] << " from
                // "<<(*ORDER)[i][4]<< " to "<< (*ORDER)[i][1]  <<" \n";
                int low = std::min((*ORDER)[i][4], (*ORDER)[i][1]);
                int high = (low == (*ORDER)[i][4]) ? (*ORDER)[i][1] : (*ORDER)[i][4];
                // fill the middle ones completely
                // quantum control has cntrl direct single line
                // classical control has doublle line
                for (int x = low + 1; x < high; ++x) {
                    while (qw[x].pos <= maxt) {
                        qw[x].line.push_back(
                            {"\\qw & ", 1});  // Append "\\qw & " with a length
                                              // increment of 1
                        qw[x].pos++;
                    }
                }

                // Fill for the lower bound
                while (qw[low].pos < maxt) {
                    qw[low].line.push_back(
                        {"\\qw & ",
                         1});  // Append "\\qw & " with a length increment of 1
                    qw[low].pos++;
                }

                // Fill for the higher bound
                while (qw[high].pos < maxt) {
                    qw[high].line.push_back(
                        {"\\qw & ",
                         1});  // Append "\\qw & " with a length increment of 1
                    qw[high].pos++;
                }

                // Add control gate for (*ORDER)[i][4]
                qw[(*ORDER)[i][4]].line.push_back({
                    "\\ctrl{" + std::to_string((*ORDER)[i][1] - (*ORDER)[i][4]) +
                        "} & ",
                    1  // Length increment
                });
                qw[(*ORDER)[i][4]].pos++;

                // Add the gate for (*ORDER)[i][1]
                qw[(*ORDER)[i][1]].pos++;
                auto tmp = drawCircuitUtility((*ORDER)[i][2], 0);
                qw[(*ORDER)[i][1]].line.push_back(
                    {tmp.first,
                     tmp.second});  // Append the result of drawCircuitUtility
                maxt++;

            } else if ((*ORDER)[i][3] == 0) {
                // (*ORDER).push_back ({ 2, t_qbit, G, c_or_q, control_bit,
                // esize});

                // [control, target] - add wires till the max
                for (int x = (*ORDER)[i][1]; x < n_qbits; x++) {
                    while (qw[x].pos < maxt) {
                        qw[x].line.push_back(
                            {"\\qw & ", 1});  // Append "\\qw & " with a length
                                              // increment of 1
                        qw[x].pos++;
                    }
                }

                for (int x = 0; x < n_cbits; x++) {
                    while (cw[x].pos < maxt) {
                        cw[x].line.push_back(
                            {"\\cw & ", 1});  // Append "\\cw & " with a length
                                              // increment of 1
                        cw[x].pos++;
                    }
                }

                // The wires between get extra "\\qw \\cwx"
                for (int x = (*ORDER)[i][1] + 1; x < n_qbits; x++) {
                    qw[x].line.push_back(
                        {"\\qw \\cwx & ", 1});  // Append "\\qw \\cwx & " with a
                                                // length increment of 1
                    qw[x].pos++;
                }

                // The control and target get specific control/wire commands
                for (int x = 0; x < n_cbits; x++) {
                    if (x != (*ORDER)[i][4]) {
                        if (x < (*ORDER)[i][4]) {
                            cw[x].line.push_back(
                                {"\\cw \\cwx & ",
                                 1});  // Append "\\cw \\cwx & " with a length
                                       // increment of 1
                        } else {
                            cw[x].line.push_back(
                                {"\\cw & ", 1});  // Append "\\cw & " with a
                                                  // length increment of 1
                        }
                        cw[x].pos++;
                    } else {
                        cw[x].line.push_back(
                            {"\\control \\cw \\cwx & ",
                             1});  // Append "\\control \\cw \\cwx & " with a
                                   // length increment of 1
                        cw[x].pos++;
                    }
                }

                qw[(*ORDER)[i][1]].pos++;
                auto tmp = drawCircuitUtility((*ORDER)[i][2], (*ORDER)[i][5]);
                qw[(*ORDER)[i][1]].line.push_back(
                    {tmp.first,
                     tmp.second});  // Append the result of drawCircuitUtility
                maxt++;
            }
        } else if ((*ORDER)[i][0] == 3)  // measure
        {
            // For the main, put wires
            while (qw[(*ORDER)[i][1]].pos < maxt) {
                qw[(*ORDER)[i][1]].line.push_back(
                    {"\\qw & ",
                     1});  // Append "\\qw & " with a length increment of 1
                qw[(*ORDER)[i][1]].pos++;
            }
            qw[(*ORDER)[i][1]].line.push_back(
                {"\\meter & ",
                 1});  // Append "\\meter & " with a length increment of 1
            qw[(*ORDER)[i][1]].pos++;

            // Q: Putting wires in between
            for (int x = (*ORDER)[i][1] + 1; x < n_qbits; x++) {
                while (qw[x].pos < maxt) {
                    qw[x].line.push_back(
                        {"\\qw & ",
                         1});  // Append "\\qw & " with a length increment of 1
                    qw[x].pos++;
                }
            }

            // Q: Edgy wires for the between
            for (int x = (*ORDER)[i][1] + 1; x < n_qbits; x++) {
                qw[x].line.push_back(
                    {"\\qw \\cwx & ", 1});  // Append "\\qw \\cwx & " with a
                                            // length increment of 1
                qw[x].pos++;
            }

            // C: Putting cwires
            for (int x = 0; x <= (*ORDER)[i][2]; ++x) {
                while (cw[x].pos < maxt) {
                    cw[x].line.push_back(
                        {"\\cw & ",
                         1});  // Append "\\cw & " with a length increment of 1
                    cw[x].pos++;
                }
            }

            // C: Putting edgy wires
            for (int x = 0; x <= (*ORDER)[i][2]; ++x) {
                cw[x].line.push_back(
                    {"\\cw \\cwx & ", 1});  // Append "\\cw \\cwx & " with a
                                            // length increment of 1
                cw[x].pos++;
            }

            maxt++;
        }
    }
    for (int i = 0; i < n_qbits; ++i) {
        while (qw[i].pos <= maxt) {
            qw[i].line.push_back(
                {"\\qw & ",
                 1});  // Append "\\qw & " with a length increment of 1
            qw[i].pos++;
        }
        qw[i].line.push_back({"\\\\ ", 0});  // Add a newline in LaTeX format
    }

    for (int i = 0; i < n_cbits; ++i) {
        while (cw[i].pos <= maxt) {
            cw[i].line.push_back(
                {"\\cw & ",
                 1});  // Append "\\cw & " with a length increment of 1
            cw[i].pos++;
        }
        cw[i].line.push_back({"\\\\ ", 0});  // Add a newline in LaTeX format
    }

handler->wr << "\\clearpage\n\\section*{Quantum Circuit}\n";
handler->wr << "\\begin{figure}[htbp]\n";
handler->wr << "    \\centering\n";
handler->wr << "    \\[\n";
handler->wr << "    \\Qcircuit @C=1em @R=.7em {\n";

    // int max_length = 0;

/*

Breaking Down Large Quantum Circuits into Manageable Units

When rendering large quantum circuits, maintaining clarity is essential. Instead of displaying an overwhelming full-length circuit, we break it into structured segments. The goal is to align each segment within a single figure, ensuring that every column of gates aligns with the widest (fat-est) block in that section. This method optimizes readability while respecting the constraints of Qcircuit (LaTeX).
How It Works

    Finding the Shortest Path
        We determine the shortest path that reaches the end of the page among the qubits q0,q1,...,qnq0​,q1​,...,qn​.
        The “length” of a path is measured in blocks (circuit elements like gates, controls, or wires).
        Most blocks have a unit length, but larger gates (e.g., UU gates) count as multiple blocks since they take extra space.
        The qubit with the shortest path to the end of the page dictates where we segment the circuit.

    Segmenting the Circuit
        Suppose the shortest path in the first iteration is 10 blocks.
        We reset the count and start a new segment, repeating Step 1 from that point onward.
        The shortest path may shift to a different qubit in the next iteration, depending on circuit structure.

    Aligning Based on the Widest Segment
        Each iteration produces a segment, but all iterations must align to the widest segment encountered. I mean shortest length score but has fat-ass elements
        This ensures that the shortest fatest segment never exceeds the page width, even if it introduces slight space inefficiencies.
        This tradeoff is necessary due to Qcircuit’s limitations in LaTeX formatting.

By structuring circuits this way, we maintain readability and logical flow while ensuring that large quantum circuits fit neatly within a figure
*/

    std::vector<int> inds;
    int csum = 0;
    inds.push_back(0);
    int starter = 0;

    for (int z = starter; z < qw[0].line.size(); z++) {
        int max_block = -1;
        for (int i = 0; i < n_qbits; i++) {
            if (max_block < qw[i].line[z].second)
                max_block = qw[i].line[z].second;
        }
        csum += max_block;
        if (csum > 15) {
            inds.push_back(z);
            csum = 0;
        }
    }

    //std::cout << "finally \n";
    // 0 5 10 17
    uint32_t min_diff = 0xffffffff;
    for (int i = 0; i < inds.size() - 1; i++) {
        min_diff = (min_diff > (inds[i + 1] - inds[i]))
                       ? ((inds[i + 1] - inds[i]))
                       : (min_diff);
    }
    int tmp = inds.size();
    inds.clear();
    for (int i = 0; i <= qw[0].line.size(); i += min_diff) {
        inds.push_back(i);
        std::cout << i << " ";
    }
    std::cout << std::endl;

    for (int i = 0; i < inds.size() - 1; i++)  // i about iterations
    {
        for (int j = 0; j < n_qbits; j++)  // j about q/c number
        {
            std::string combinedLine = qw[j].line[0].first;
            for (int k = inds[i]; k < inds[i + 1];
                 k++)  // k about line indicies
                if (k) combinedLine += qw[j].line[k].first;

            handler->wr << combinedLine << " \\qw & \\\\" << std::endl;
        }
        for (int j = 0; j < n_cbits; j++) {
            std::string combinedLine = cw[j].line[0].first;
            for (int k = inds[i]; k < inds[i + 1]; k++)
                if (k) combinedLine += cw[j].line[k].first;

            handler->wr << combinedLine << " \\cw & \\\\" << std::endl;
        }
        handler->wr << "\\\\ \n\\\\ \n";
    }
    for (int i = 0; i < n_qbits; i++) {
        std::string combinedLine = qw[i].line[0].first;
        for (int j = inds[inds.size() - 1]; j < qw[0].line.size(); j++)
            if (j) combinedLine += qw[i].line[j].first;

        handler->wr << combinedLine << std::endl;
    }
    for (int i = 0; i < n_cbits; i++) {
        std::string combinedLine = cw[i].line[0].first;
        for (int j = inds[inds.size() - 1]; j < qw[0].line.size(); j++)
            if (j) combinedLine += cw[i].line[j].first;

        handler->wr << combinedLine << std::endl;
    }
    handler->wr << "\\\\ \n\\\\ \n";

    // std::cout << "}\n\\end{displaymath}\n";
    if (handler->wr.good()) {
        handler->wr << "}\n\\]\n";
        if (name != "") handler->wr << "\\caption{" << name << "}\n";
        handler->wr << "\\end{figure}\n";
        handler->wr.flush();
    } else
        std::cerr << "writing issue";

    // std::cout<<"CIRCUIT maxt: "<<maxt << " len: "<< max_length << std::endl;
    // handler->writeInPdf("CIRCUIT maxt: "+ std::to_string(maxt)+ " len: " +
    // std::to_string(max_length));
}

bool Qcircuit::accessCreg(int cbit) {
    return (c_reg & (1 << (n_cbits - 1 - cbit)));
}

// std::vector<struct experiment> simulation()
Experiment Qcircuit::simulation() {
    if (handler->circuittDrawn == false) {
        drawCircuit();
        handler->circuittDrawn = true;
    }
    handler->qasmgen=false;
    Experiment exp1;
    exp1.c_reg_value = c_reg;
    exp1.final_state = state;
    exp1.n_cbits = static_cast<uint64_t>(n_cbits);
    if (!handler->wr.is_open())
        std::cerr << "already closed before returning shot";
    // handler->wr = NULL;
    return exp1;
}

struct result simulate(mimiqHandler* handler,
                       std::function<Qcircuit::experiment(mimiqHandler*)> func,
                       int shots) {
    if (!func) {
        std::cerr << "NULL Experiment error \n";
        exit(0);
    }
    struct result res;
    Experiment shot_result;

    while (shots--) {
        shot_result = func(handler);
        if (!handler->wr.is_open()) std::cerr << "end of shot closed";
        res.m[shot_result.c_reg_value]++;
    }

    // last shot results
    res.n_cbits = shot_result.n_cbits;

    res.creg.resize(res.n_cbits);
    for (uint64_t i = 0; i < res.n_cbits; ++i)
        res.creg[shot_result.n_cbits - 1 - i] =
            (shot_result.c_reg_value >> i) & 1;
    res.state = shot_result.final_state;  // TODO final state doesnt make sense
    // res.(*state).print();
    res.handler = handler;
    return res;
}

// all basisc gates like unitary gate, hadamard gate, pauli's x,y,z are all
// derieved from this unitary gate which takes euler angles as input the
// mathematical formula for this matrix is based on qiskit's implementation
std::vector<std::vector<Coeff>> GATES::unitaryGate(double theta, double psi,
                                                   double lam) {
    return {
        {{cos(theta / 2), 0},
         {-cos(lam) * sin(theta / 2), -sin(lam) * sin(theta / 2)}},
        {{cos(psi) * sin(theta / 2), sin(psi) * sin(theta / 2)},
         {cos(theta / 2) * cos(lam + psi), cos(theta / 2) * sin(lam + psi)}}};
}

// while applying a nomral independant gate can be generalized, controlled gates
// are not. Is that even possible to do os ? controlled pauli X gate - when
// control qbit is 1 on {|0>, |1>} basis, then target qbit is flipped.
void GATES::controlled_pauli_X(struct state_vector& sv, int cbit, int tbit) {
    // the cnd is just the number of states possible from n_qbits, which is
    // 2^n_qbits the mask is created from the specified positions of control bit
    // and target bit logical operation with the mask will result in 1 or 0
    int cnd = 1 << sv.n_qbits, cmask = 1 << (sv.n_qbits - 1 - cbit),
        tmask = 1 << (sv.n_qbits - 1 - tbit), other;
    std::vector<bool> visit(cnd, false);
    // the reason why I have a visit array is because I want to swap only once
    // for a given set 1xxt where t(target bit) is 0 or 1 and the most
    // signioficant 1 is assumed to be the control bit so if I encounter either
    // of 1xx0 or 1xx1, I will swap the Coeffs and then I will set flag of both
    // states as true to prevent re-swapping which will be pointless
    double tmp1, tmp2;
    for (int i = 0; i < cnd; i++) {
        if ((i & cmask) &&
            (visit[i] ==
             false)) {      // when unvisited and mask condition satisfies --
                            // meaning the position is right
            if (i & tmask)  // if tbit 1 -> other = xxxx & 11011 => 0 at tbit
                other = i & ~(tmask);  // retrieving the basis state of the
                                       // complement of the current 1xxt
            else
                other = i | tmask;  // else, other = xxxx | 00100 => 1 at tbit

            // swapping
            tmp1 = sv.coeffs[i].real;
            tmp2 = sv.coeffs[i].complex;
            sv.coeffs[i].real = sv.coeffs[other].real;
            sv.coeffs[i].complex = sv.coeffs[other].complex;
            sv.coeffs[other].real = tmp1;
            sv.coeffs[other].complex = tmp2;
            visit[i] = visit[other] = true;
        }
    }
}
void GATES::toffoli(struct state_vector& sv, int cbit1, int cbit2, int tbit) {
    // the cnd is just the number of states possible from n_qbits, which is
    // 2^n_qbits the mask is created from the specified positions of control bit
    // and target bit logical operation with the mask will result in 1 or 0
    int cnd = 1 << sv.n_qbits, cmask1 = 1 << (sv.n_qbits - 1 - cbit1),
        cmask2 = 1 << (sv.n_qbits - 1 - cbit2),
        tmask = 1 << (sv.n_qbits - 1 - tbit), other;
    std::vector<bool> visit(cnd, false);
    // the reason why I have a visit array is because I want to swap only once
    // for a given set 1xxt where t(target bit) is 0 or 1 and the most
    // signioficant 1 is assumed to be the control bit so if I encounter either
    // of 1xx0 or 1xx1, I will swap the Coeffs and then I will set flag of both
    // states as true to prevent re-swapping which will be pointless
    double tmp1, tmp2;
    for (int i = 0; i < cnd; i++) {
        if (((i & cmask1) && (i & cmask2)) &&
            (visit[i] ==
             false)) {      // when unvisited and mask condition satisfies --
                            // meaning the position is right
            if (i & tmask)  // if tbit 1 -> other = xxxx & 11011 => 0 at tbit
                other = i & ~(tmask);  // retrieving the basis state of the
                                       // complement of the current 1xxt
            else
                other = i | tmask;  // else, other = xxxx | 00100 => 1 at tbit

            // swapping
            tmp1 = sv.coeffs[i].real;
            tmp2 = sv.coeffs[i].complex;
            sv.coeffs[i].real = sv.coeffs[other].real;
            sv.coeffs[i].complex = sv.coeffs[other].complex;
            sv.coeffs[other].real = tmp1;
            sv.coeffs[other].complex = tmp2;
            visit[i] = visit[other] = true;
        }
    }
}

// controlled hadamard-  qbit in control and target also qbit, along with state
// vector of course
void GATES::controlled_hadamard(struct state_vector& sv, int cbit, int tbit) {
    // masking as explaine in controlled pauli x
    int cnd = 1 << sv.n_qbits, cmask = 1 << (sv.n_qbits - 1 - cbit),
        tmask = 1 << (sv.n_qbits - 1 - tbit), other;
    // I am not directly updating the passed state vector because I need to use
    // all the states passed before updating any, so I just create a new vector
    // and assign it in the end
    struct state_vector tmp(sv.n_qbits);
    for (int i = 0; i < cnd; i++) {
        if (i & cmask) {
            // std::cout << "in "<< i << std::endl ;
            //  in high level, we apply hadamard gate to a qbit. A single qbit,
            //  hence whatever is the Coeff value, eg. C|1> or C|0> it will be
            //  1/root2 * C|0> +/- 1/root2* C|1>
            // so in the new vector, we add the 2 parts in their respective
            // places
            struct Coeff coeff1(1 / root2, 0);  // 1/root2
            if (i & tmask)  // if target bit is 1 here, x|1>x --> 1/root2 *[
                            // x|0>x - x|1>x ]
            {
                other = i & ~(tmask);
                tmp.coeffs[i] = tmp.coeffs[i] - (coeff1)*sv.coeffs[i];
            } else  // if target bit is 0 here, x|0>x --> 1/root2 *[ x|0>x +
                    // x|1>x ]
            {
                other = i | tmask;
                tmp.coeffs[i] = tmp.coeffs[i] + (coeff1)*sv.coeffs[i];
            }
            tmp.coeffs[other] = tmp.coeffs[other] + (coeff1)*sv.coeffs[i];
        }
        // if not involved in control bit, just assign whats existing as there
        // wont be any change required
        else
            tmp.coeffs[i] = sv.coeffs[i];
    }
    // assign to make changes
    sv = tmp;
    return;
}

// simple 2x2 matrix multiplication
void GATES::matrixmultiply2x2(const std::vector<std::vector<Coeff>>& g1,
                              const std::vector<std::vector<Coeff>>& g2,
                              std::vector<std::vector<Coeff>>& g3) {
    g3[0][0] = (g1[0][0] * g2[0][0]) + (g1[0][1] * g2[1][0]);
    g3[0][1] = (g1[0][0] * g2[0][1]) + (g1[0][1] * g2[1][1]);
    g3[1][0] = (g1[1][0] * g2[0][0]) + (g1[1][1] * g2[1][0]);
    g3[1][1] = (g1[1][0] * g2[0][1]) + (g1[1][1] * g2[1][1]);
}

// when control bit is 1, then target bit coeff should be multiplied with (-1 +
// 0i)
void GATES::controlled_pauli_Z(struct state_vector& sv, int cbit, int tbit) {
    int cnd = 1 << sv.n_qbits, cmask = 1 << (sv.n_qbits - 1 - cbit),
        tmask = 1 << (sv.n_qbits - 1 - tbit);
    struct Coeff tmp(-1, 0);
    for (int i = 0; i < cnd; i++)
        if ((i & cmask) && (i & tmask))  // if cbit and tbit 1 flip sign.
            sv.coeffs[i] = sv.coeffs[i] * tmp;
}

// to print matrix
void GATES::printMatrix(const std::vector<std::vector<Coeff>>& matrix) {
    if (matrix.size() == 0) {
        std::cout << "empty " << std::endl;
        return;
    }
    for (const auto& row : matrix) {
        for (const auto& element : row) {
            std::cout << "{" << element.real << ", " << element.complex << "} ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
// kroneckerSUB is the utility of KnroneckerBUFF, it will create a cross product
// (or kronecker product ) of 2 gates and assign it to m3
void GATES::kroneckerSUB(const std::vector<std::vector<Coeff>>& m1,
                         const std::vector<std::vector<Coeff>>& m2,
                         std::vector<std::vector<Coeff>>& m3) {
    // Handle edge case
    if (m1.empty()) {
        m3 = m2;
        return;
    }

    // Precompute dimensions
    size_t m1_rows = m1.size();
    size_t m1_cols = m1[0].size();
    size_t m2_rows = m2.size();
    size_t m2_cols = m2[0].size();

    size_t m3_rows = m1_rows * m2_rows;
    size_t m3_cols = m1_cols * m2_cols;

    // Allocate memory for m3
    m3.resize(m3_rows);
    for (auto& row : m3) {
        row.resize(m3_cols);
    }

    // Compute Kronecker product
    for (size_t i = 0; i < m1_rows; ++i) {
        for (size_t j = 0; j < m1_cols; ++j) {
            for (size_t k = 0; k < m2_rows; ++k) {
                for (size_t l = 0; l < m2_cols; ++l) {
                    m3[i * m2_rows + k][j * m2_cols + l] = m1[i][j] * m2[k][l];
                }
            }
        }
    }
}

// when a gate has to be applied to a qbit.. but we have a state vector of
// multiple states of multiple qbits, we need to build a matrix which when
// multiplied with the state vector, the areas of state_vector which correspond
// to the atrget qbit only has to be changed. the kroneckerbuff will create the
// matrix. for a statevector which has 5 qbits like xxxxx, if we want to apply a
// unitary gate U to 2nd qbit, then our matrix would like I x U x I x I x I
// which is multiplying identity gates using cross product in appropriate order
void GATES::kroneckerBUFF(const std::vector<std::vector<Coeff>>& gate,
                          std::vector<std::vector<Coeff>>& res, int n_qbits,
                          int tbit) {
    if (n_qbits == 1) {
        // if only 1 qbit, then no need of buffing so just directly return gate
        res = gate;
        return;
    }
    std::vector<std::vector<Coeff>>
        tmp;  // tmp will be buffed from empty, and will be assigned to res on
              // completion
    for (int i = 0; 1;) {
        // loop for n_qbit times, cross product with Identity every iteration
        // except for the position of target qbit
        if (i == tbit)
            kroneckerSUB(tmp, gate, res);
        else
            kroneckerSUB(tmp, GATES::unitaryGate(0, 0, 0), res);
        if (++i < n_qbits)
            tmp = res;
        else
            break;
    }
    return;
}

// given a 2x2 mateix of complex numbers (or type Coeff), this function will
// return the inverse of it
std::vector<std::vector<Coeff>> GATES::inverse2x2(
    const std::vector<std::vector<Coeff>>& matrix) {
    std::vector<std::vector<Coeff>> res;
    res.resize(2, std::vector<Coeff>(2));
    Coeff num,
        det = ((matrix[0][0] * matrix[1][1]) - (matrix[0][1] * matrix[1][0]));
    num.real = det.real;
    num.complex = -det.complex;
    auto denom = (det.real * det.real) + (det.complex * det.complex);
    num.real /= denom;
    num.complex /= denom;
    // std::cout<<"1/det: "<<num.real<<std::endl;
    // 1 / a + i b
    // a - i b /  a^2 - (ib)^2
    res[0][0] = matrix[1][1];
    res[0][1].real = -matrix[0][1].real;
    res[0][1].complex = -matrix[0][1].complex;
    res[1][0].real = -matrix[1][0].real;
    res[1][0].complex = -matrix[1][0].complex;
    res[1][1].real = matrix[0][0].real;
    res[1][1].complex = (matrix[0][0].complex);
    res[0][0] = num * res[0][0];
    res[0][1] = num * res[0][1];
    res[1][0] = num * res[1][0];
    res[1][1] = num * res[1][1];
    return res;
}

void Qcircuit::applyUtility(
    const std::string& gate, int& G, int& esize,
    std::vector<double> arr = std::vector<double>(),
    std::vector<std::vector<Coeff>>& matrix = empty_matrix) {
    if (gate == "H" || gate == "h") {
        G = 1;
        matrix = GATES::unitaryGate(PI / 2, 0, PI);
    } else if (gate == "x" || gate == "X") {
        G = 2;
        matrix = GATES::unitaryGate(PI, 0, PI);
    } else if (gate == "y" || gate == "Y") {
        G = 3;
        matrix = GATES::unitaryGate(PI, PI / 2, PI / 2);
    } else if (gate == "z" || gate == "Z") {
        G = 4;
        matrix = GATES::unitaryGate(0, 0, PI);
    } else if (gate == "s" || gate == "S") {
        G = 5;
        matrix = GATES::unitaryGate(0, 0, PI / 2);
    } else if (gate == "sd" || gate == "Sd") {
        G = 6;
        matrix = GATES::unitaryGate(0, 0, -PI / 2);
    } else if (gate == "T" || gate == "t") {
        G = 7;
        matrix = GATES::unitaryGate(0, 0, PI / 4);
    } else if (gate == "Td" || gate == "td") {
        G = 8;
        matrix = GATES::unitaryGate(0, 0, -PI / 4);
    } else if (gate == "rx" || gate == "Rx") {
        G = 105;
        esize = (*euler_container).size();
        (*euler_container).push_back({arr[0], 0, 0});

        matrix = GATES::unitaryGate((*euler_container)[esize][0], -PI / 2, PI / 2);
    } else if (gate == "ry" || gate == "Ry") {
        G = 106;
        esize = (*euler_container).size();
        (*euler_container).push_back({arr[0], 0, 0});

        matrix = GATES::unitaryGate((*euler_container)[esize][0], 0, 0);
    } else if (gate == "U3" || gate == "u3" || gate == "u" || gate == "U") {
        esize = (*euler_container).size();
        (*euler_container).push_back({arr[0], arr[1], arr[2]});
        G = 101;
        matrix = GATES::unitaryGate((*euler_container)[esize][0],
                                    (*euler_container)[esize][1],
                                    (*euler_container)[esize][2]);
    } else if (gate == "IU3" || gate == "iu3" || gate == "IU" ||
               gate == "iu")  // todo: generalise u1 u2 u3 inverse
    {
        esize = (*euler_container).size();
        (*euler_container).push_back({arr[0], arr[1], arr[2]});
        G = 102;
        matrix = GATES::inverse2x2(GATES::unitaryGate(
            (*euler_container)[esize][0], (*euler_container)[esize][1],
            (*euler_container)[esize][2]));
    } else if (gate == "U1" || gate == "u1") {
        esize = (*euler_container).size();
        (*euler_container).push_back({arr[0], 0, 0});
        G = 103;
        matrix = GATES::unitaryGate(0, 0, (*euler_container)[esize][0]);
    } else if (gate == "U2" || gate == "u2") {
        esize = (*euler_container).size();
        (*euler_container).push_back({arr[0], arr[1], 0});
        G = 104;
        matrix = GATES::unitaryGate(PI / 2, (*euler_container)[esize][0],
                                    (*euler_container)[esize][1]);
    }
}

void Qcircuit::applyGate(const std::string& gate, int t_qbit,
                         std::vector<double> arr, bool ORDERinc) {
    std::vector<std::vector<Coeff>> matrix, buffedmatrix;
    int esize = -1, G;

    applyUtility(gate, G, esize, arr, matrix);

    // std::cout<<"p: "<<1 <<" "<<t_qbit << " " << G << " " << esize <<
    // std::endl;
    if (ORDERinc) (*ORDER).push_back({1, t_qbit, G, esize});  // send index

    GATES::kroneckerBUFF(matrix, buffedmatrix, n_qbits, t_qbit);
    int n = 1 << n_qbits;
    std::vector<Coeff> result(n, {0, 0});
    for (int ii = 0; ii < n; ++ii) {
        for (int j = 0; j < n; ++j)
            result[ii] = result[ii] + (buffedmatrix[ii][j] * (*state).coeffs[j]);
    }
    (*state).coeffs = result;
}
void Qcircuit::applyControlledGate(const std::string& gate,
                                   const std::string c_qbit_str, int t_qbit,
                                   std::vector<double> arr) {
    std::istringstream iss(c_qbit_str);
    char cq1, cq2, d = '0';
    int control_bit1, control_bit2, c_or_q1, c_or_q2, G;
    iss >> cq1 >> control_bit1 >> d >> cq2 >> control_bit2;

    if (d == '0') { // WHEN DELIMETER NOT RECEIVED, MEANING SINGLE CONTROL, not toffoli
        c_or_q1 = (cq1 == 'q') ? 1 : 0;

        if (c_or_q1 == 0) {
            
            if(handler->qasmgen) handler->oqsm +=             
            "if ( c == " + std::to_string(1 << control_bit1) + ") " + gate +  " q[" + std::to_string(t_qbit) + "];\n";
            // todo gate validate the user written gate 


            if (c_reg & (1 << (n_cbits - 1 - control_bit1)))
                applyGate(gate, t_qbit, arr, false);
        }

        int esize;
        applyUtility(gate, G, esize, arr);

        (*ORDER).push_back({2, t_qbit, G, c_or_q1, control_bit1, esize});

        if (c_or_q1 == 1) {  // CONTROL ON ZERO WITH Q needs to be done. TODO
            if (G == 1)      // gate
                GATES::controlled_hadamard((*state), control_bit1,
                                           t_qbit);  // qbit, tqbit

            else if (G == 2)
                GATES::controlled_pauli_X((*state), control_bit1, t_qbit);

            else if (G == 4)  // 4 - pauli z , 3 - pauli y
                GATES::controlled_pauli_Z((*state), control_bit1, t_qbit);
        }
    } else if (d == ',') {
        c_or_q1 = (cq1 == 'q') ? 1 : 0;
        c_or_q2 = (cq2 == 'q') ? 1 : 0;
        if ((c_or_q1 == c_or_q2) && (c_or_q1 == 0)) {
            if ((c_reg & (1 << (n_cbits - 1 - control_bit1))) &&
                (c_reg & (1 << (n_cbits - 1 - control_bit2))))
                applyGate(gate, t_qbit, arr, false);
        }
        int esize;
        applyUtility(gate, G, esize, arr);
        // circuit
        (*ORDER).push_back(
            {4, t_qbit, G, c_or_q1, control_bit1, control_bit2, esize});

        if ((c_or_q1 == c_or_q2) &&
            (c_or_q1 == 1)) {  // CONTROL ON ZERO WITH Q needs to be done. TODO
            if (G == 2)
                GATES::toffoli((*state), control_bit1, control_bit2, t_qbit);
        }
    }
}



    void Qcircuit::u( int t_qbit, std::vector<double> arr, bool ORDERinc)
    {
        applyGate("u",t_qbit,arr,ORDERinc);
        if(handler->qasmgen) handler->oqsm += "u("+ std::to_string(arr[0]) + ","  + std::to_string(arr[1]) + ","  + std::to_string(arr[2])  + ") q[" + std::to_string(t_qbit) + "];\n";
    }
    void Qcircuit::x( int t_qbit, std::vector<double> arr, bool ORDERinc)
    {
        applyGate("x",t_qbit,arr,ORDERinc);
        if(handler->qasmgen) handler->oqsm += "x q[" + std::to_string(t_qbit) + "];\n";
    }
    void Qcircuit::y( int t_qbit, std::vector<double> arr, bool ORDERinc)
    {
        applyGate("y",t_qbit,arr,ORDERinc);
        if(handler->qasmgen) handler->oqsm += "y q[" + std::to_string(t_qbit) + "];\n";
    }
    void Qcircuit::z( int t_qbit, std::vector<double> arr, bool ORDERinc)
    {
        applyGate("z",t_qbit,arr,ORDERinc);
        if(handler->qasmgen) handler->oqsm += "z q[" + std::to_string(t_qbit) + "];\n";
    }
    void Qcircuit::h( int t_qbit, std::vector<double> arr, bool ORDERinc)
    {
        applyGate("h",t_qbit,arr,ORDERinc);
        if(handler->qasmgen) handler->oqsm += "h q[" + std::to_string(t_qbit) + "];\n";
    }
    void Qcircuit::u2( int t_qbit, std::vector<double> arr, bool ORDERinc)
    {
        applyGate("u2",t_qbit,arr,ORDERinc);
        if(handler->qasmgen) handler->oqsm += "u(pi/2,"+ std::to_string(arr[0]) + ","  + std::to_string(arr[1]) + ") q[" + std::to_string(t_qbit) + "];\n";
    }
    void Qcircuit::u1( int t_qbit, std::vector<double> arr, bool ORDERinc)
    {
        applyGate("u",t_qbit,arr,ORDERinc);
        if(handler->qasmgen) handler->oqsm += "u(0,0,"+ std::to_string(arr[0])  + ") q[" + std::to_string(t_qbit) + "];\n";
    }
    void Qcircuit::t( int t_qbit, std::vector<double> arr, bool ORDERinc)
    {
        applyGate("t",t_qbit,arr,ORDERinc);
        if(handler->qasmgen) handler->oqsm += "t q[" + std::to_string(t_qbit) + "];\n";
    }
    void Qcircuit::tdg( int t_qbit, std::vector<double> arr, bool ORDERinc)
    {
        applyGate("td",t_qbit,arr,ORDERinc);
        if(handler->qasmgen) handler->oqsm += "tdg q[" + std::to_string(t_qbit) + "];\n";
    }
    void Qcircuit::s( int t_qbit, std::vector<double> arr, bool ORDERinc)
    {
        applyGate("s",t_qbit,arr,ORDERinc);
        if(handler->qasmgen) handler->oqsm += "s q[" + std::to_string(t_qbit) + "];\n";
    }
    void Qcircuit::sdg( int t_qbit, std::vector<double> arr, bool ORDERinc )
    {
        applyGate("sd",t_qbit,arr,ORDERinc);
        if(handler->qasmgen) handler->oqsm += "sdg q[" + std::to_string(t_qbit) + "];\n";
    }
    void Qcircuit::ry( int t_qbit, std::vector<double> arr, bool ORDERinc )
    {
        applyGate("ry",t_qbit,arr,ORDERinc);
        if(handler->qasmgen) handler->oqsm += "ry("+ std::to_string(arr[0])  + ") q[" + std::to_string(t_qbit) + "];\n";
    }
    void Qcircuit::toffoli(int c1_qbit, int c2_qbit, int t_qbit)
    {
        applyControlledGate("x",("q" + std::to_string(c1_qbit) + ",q" + std::to_string(c2_qbit)),t_qbit );
        if(handler->qasmgen) handler->oqsm += "ccx q[" + std::to_string(c1_qbit) + "], q[" + std::to_string(c2_qbit) + "], q["  + std::to_string(t_qbit) +  "];\n";
    }
    void Qcircuit::ccx(int c1_qbit, int c2_qbit, int t_qbit)
    {
        applyControlledGate("x",("q" + std::to_string(c1_qbit) + ",q" + std::to_string(c2_qbit)),t_qbit );
        if(handler->qasmgen) handler->oqsm += "ccx q[" + std::to_string(c1_qbit) + "], q[" + std::to_string(c2_qbit) + "], q["  + std::to_string(t_qbit) +  "];\n";
    }
    void Qcircuit::cx(int cbit, int qbit)
    {
        applyControlledGate("x","q"+std::to_string(cbit), qbit);
        if(handler->qasmgen) handler->oqsm += "cx q[" + std::to_string(cbit) + "], q[" + std::to_string(qbit) + "];\n";
    }
    void Qcircuit::ch(int cbit, int qbit)
    {
        applyControlledGate("h","q"+std::to_string(cbit), qbit);
        if(handler->qasmgen) handler->oqsm += "ch q[" + std::to_string(cbit) + "], q[" + std::to_string(qbit) + "];\n";
    }
