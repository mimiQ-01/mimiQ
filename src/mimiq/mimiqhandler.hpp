#pragma once 

#include <fstream>

#define R3(x) (std::round((x) * 1000.0) / 1000.0)
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

