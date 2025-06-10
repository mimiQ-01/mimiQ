#include "mimiqhandler.hpp"
#include <iostream>
#include <ctime>
MimiqHandler::MimiqHandler(std::string path) {
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
    dir_path = path;
    latex_writer.open(dir_path + "report.tex");
latex_writer << "\\documentclass{article}\n";
latex_writer << "\\usepackage[margin=1in]{geometry}\n";
latex_writer << "\\usepackage{datetime}\n";
latex_writer << "\\usepackage{qcircuit}\n";
latex_writer << "\\usepackage{amsmath}\n";
latex_writer << "\\usepackage{pgfplots}\n";
latex_writer << "\\pgfplotsset{compat=1.18}\n";
latex_writer << "\\usepackage[utf8]{inputenc}\n";
latex_writer << "\\begin{document}\n";
latex_writer << "\\begin{center}\n";
latex_writer << "    {\\LARGE \\textbf{mimiQ++ Report}} \\\\\n";
latex_writer << "    \\large \\today \\quad \\currenttime\n";
latex_writer << "\\end{center}\n";
latex_writer << "\\hrule\n";
latex_writer << "\\vspace{1cm}\n";

    circuit_drawn = false;
    qasm_gen = true;
    can_qasm= true;
    oqsm += "OPENQASM 2.0;\ninclude \"qelib1.inc\";\n";

}
void MimiqHandler::clean() { circuit_drawn = false; qasm_gen = true; can_qasm = true; oqsm = "";}

void MimiqHandler::write_in_pdf(std::string msg) {
    if (latex_writer.is_open() && latex_writer.good()) latex_writer << msg << "\n\n";
}

void MimiqHandler::generate_report() {
    if (latex_writer.is_open()) {
        report_generated = true;
        latex_writer << "\\end{document}";
        latex_writer.close();
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
