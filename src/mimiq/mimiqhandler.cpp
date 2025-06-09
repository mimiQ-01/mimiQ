#include "mimiqhandler.hpp"
#include <iostream>
#include <ctime>
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
