#pragma once 

#include <fstream>

#define R3(x) (std::round((x) * 1000.0) / 1000.0)
struct MimiqHandler 
{
    bool circuit_drawn, report_generated, qasm_gen , can_qasm ;
    std::ofstream latex_writer, // to write Latex for pdf output
    debug_writer;
    std::string oqsm, // to write openQASM 2.0 code
    dir_path ; // where pdf should be stored

    MimiqHandler(std::string path = "");
    void clean();
    void write_in_pdf(std::string msg);
    void generate_report();
};

