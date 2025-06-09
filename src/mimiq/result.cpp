#include "result.hpp"
#include "state_vector.hpp"
#include "qcircuit.hpp"
#include "result.hpp"
#include <string>
#include <algorithm>
#include <vector>
#include <sstream>
#include <functional>
#include <iostream>
#include <map>
#include <cstdint>

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

void result::print_counts() const {
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

std::pair<int, int> result::get_counts_of(int cbit){
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
