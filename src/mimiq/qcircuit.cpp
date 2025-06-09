#include "qcircuit.hpp"
#include "gates.hpp"
#include <iostream>
#include <cmath>
#include <array>
#include <algorithm>
#include <sstream>

#define PI acos(0.0) * 2
#define R3(x) (std::round((x) * 1000.0) / 1000.0)
static std::vector<std::vector<Coeff> > empty_matrix ; // for applyUtility
struct wire {
    std::vector<std::pair<std::string, int>> line;
    int pos;
    bool nature;  // TODO
};

void Qcircuit::printVector() const{
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

void Qcircuit::print_creg() const 
{ std::cout << "creg: " << c_reg << std::endl; }

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


std::pair<std::string, int> Qcircuit::drawCircuitUtility(int key,
                                                         int eulerIndex) {
    int length = 1;
    // std::cout<<"received with key: "<< key <<" , index: "<< eulerIndex ;
    std::string res;
    auto& euler_array = 
    
    (*euler_container)[eulerIndex];  // Explicitly reference the array
    switch (key) {
        case GATES::gateNumber::h :
            res = "\\gate{H} & ";
            break;
        case GATES::gateNumber::x:
            res = "\\gate{X} & ";
            break;
        case GATES::gateNumber::y:
            res = "\\gate{Y} & ";
            break;
        case GATES::gateNumber::z:
            res = "\\gate{Z} & ";
            break;
        case GATES::gateNumber::s:
            res = "\\gate{S} & ";
            break;
        case GATES::gateNumber::sd:
            res = "\\gate{S^\\dag} & ";
            break;
        case GATES::gateNumber::t:
            res = "\\gate{T} & ";
            break;
        case GATES::gateNumber::td:
            res = "\\gate{T^\\dag} & ";
            break;
        case GATES::gateNumber::u3:
            res =
                "\\gate{U(" + trimZeroes(std::to_string(
                        R3(
                            (*euler_container)[eulerIndex][0]
                        )
                    )) +
                "," +
                trimZeroes(std::to_string(R3((*euler_container)[eulerIndex][1]))) +
                "," +
                trimZeroes(std::to_string(R3((*euler_container)[eulerIndex][2]))) +
                ")} & ";
            length += 3;
            break;
        case GATES::gateNumber::iu3:
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
        case GATES::gateNumber::u1:
            res =
                "\\gate{U1(" +
                trimZeroes(std::to_string(R3((*euler_container)[eulerIndex][0]))) +
                ")} & ";
            length++;
            break;
        case GATES::gateNumber::u2:
            res = ("\\gate{U2(" +
                   trimZeroes(
                       std::to_string(R3((*euler_container)[eulerIndex][0]))) +
                   "," +
                   trimZeroes(
                       std::to_string(R3((*euler_container)[eulerIndex][1]))) +
                   ")} & ");
            length += 2;
            break;

        case GATES::gateNumber::rx:
            res =
                "\\gate{Rx(" +
                trimZeroes(std::to_string(R3((*euler_container)[eulerIndex][0]))) +
                ")} & ";
            length++;
            break;
        case GATES::gateNumber::ry:
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



void Qcircuit::applyUtility(
    const std::string& gate, int& G, int& esize,
    std::vector<double> arr = std::vector<double>(),
    std::vector<std::vector<Coeff>>& matrix = empty_matrix) {
    if (gate == "H" || gate == "h") {
        G =  GATES::gateNumber::h ;
        matrix = GATES::unitaryGate(PI / 2, 0, PI);
    } else if (gate == "x" || gate == "X") {
        G = GATES::gateNumber::x;
        matrix = GATES::unitaryGate(PI, 0, PI);
    } else if (gate == "y" || gate == "Y") {
        G = GATES::gateNumber::y;
        matrix = GATES::unitaryGate(PI, PI / 2, PI / 2);
    } else if (gate == "z" || gate == "Z") {
        G = GATES::gateNumber::z;
        matrix = GATES::unitaryGate(0, 0, PI);
    } else if (gate == "s" || gate == "S") {
        G = GATES::gateNumber::s;
        matrix = GATES::unitaryGate(0, 0, PI / 2);
    } else if (gate == "sd" || gate == "Sd") {
        G = GATES::gateNumber::sd;
        matrix = GATES::unitaryGate(0, 0, -PI / 2);
    } else if (gate == "T" || gate == "t") {
        G = GATES::gateNumber::t;
        matrix = GATES::unitaryGate(0, 0, PI / 4);
    } else if (gate == "Td" || gate == "td") {
        G = GATES::gateNumber::td;
        matrix = GATES::unitaryGate(0, 0, -PI / 4);
    } else if (gate == "rx" || gate == "Rx") {
        G = GATES::gateNumber::rx;
        esize = (*euler_container).size();
        (*euler_container).push_back({arr[0], 0, 0});

        matrix = GATES::unitaryGate((*euler_container)[esize][0], -PI / 2, PI / 2);
    } else if (gate == "ry" || gate == "Ry") {
        G = GATES::gateNumber::ry;
        esize = (*euler_container).size();
        (*euler_container).push_back({arr[0], 0, 0});

        matrix = GATES::unitaryGate((*euler_container)[esize][0], 0, 0);
    } else if (gate == "U3" || gate == "u3" || gate == "u" || gate == "U") {
        esize = (*euler_container).size();
        (*euler_container).push_back({arr[0], arr[1], arr[2]});
        G = GATES::gateNumber::u3;
        matrix = GATES::unitaryGate((*euler_container)[esize][0],
                                    (*euler_container)[esize][1],
                                    (*euler_container)[esize][2]);
    } else if (gate == "IU3" || gate == "iu3" || gate == "IU" ||
               gate == "iu")  // todo: generalise u1 u2 u3 inverse
    {
        esize = (*euler_container).size();
        (*euler_container).push_back({arr[0], arr[1], arr[2]});
        G = GATES::gateNumber::iu3;
        matrix = GATES::inverse2x2(GATES::unitaryGate(
            (*euler_container)[esize][0], (*euler_container)[esize][1],
            (*euler_container)[esize][2]));
    } else if (gate == "U1" || gate == "u1") {
        esize = (*euler_container).size();
        (*euler_container).push_back({arr[0], 0, 0});
        G = GATES::gateNumber::u1;
        matrix = GATES::unitaryGate(0, 0, (*euler_container)[esize][0]);
    } else if (gate == "U2" || gate == "u2") {
        esize = (*euler_container).size();
        (*euler_container).push_back({arr[0], arr[1], 0});
        G = GATES::gateNumber::u2;
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
