#include "../qcircuit.hpp"
#include "../gates.hpp"
#include <array>
#include <functional>
#include <cmath>
#include <iostream>

#define R3(x) (std::round((x) * 1000.0) / 1000.0)
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


std::pair<std::string, int> Qcircuit::draw_circuit_utility_(int key,
                                                         int eulerIndex) {
    int length = 1;
    // std::cout<<"received with key: "<< key <<" , index: "<< eulerIndex ;
    std::string res;
    auto& euler_array = (*euler_container_)[eulerIndex];  // Explicitly reference the array
    switch (key) {
        case gates::gateNumber::h :
            res = "\\gate{H} & ";
            break;
        case gates::gateNumber::x:
            res = "\\gate{X} & ";
            break;
        case gates::gateNumber::y:
            res = "\\gate{Y} & ";
            break;
        case gates::gateNumber::z:
            res = "\\gate{Z} & ";
            break;
        case gates::gateNumber::s:
            res = "\\gate{S} & ";
            break;
        case gates::gateNumber::sd:
            res = "\\gate{S^\\dag} & ";
            break;
        case gates::gateNumber::t:
            res = "\\gate{T} & ";
            break;
        case gates::gateNumber::td:
            res = "\\gate{T^\\dag} & ";
            break;
        case gates::gateNumber::u3:
            res =
                "\\gate{U(" + trimZeroes(std::to_string(
                        R3(euler_array[0])
                    )) +
                "," +
                trimZeroes(std::to_string(R3((*euler_container_)[eulerIndex][1]))) +
                "," +
                trimZeroes(std::to_string(R3((*euler_container_)[eulerIndex][2]))) +
                ")} & ";
            length += 3;
            break;
        case gates::gateNumber::iu3:
            res =
                "\\gate{IU3(" +
                trimZeroes(std::to_string(R3(euler_array[0]))) +
                "," +
                trimZeroes(std::to_string(R3(euler_array[1]))) +
                "," +
                trimZeroes(std::to_string(R3(euler_array[2]))) +
                ")} & ";
            length += 3;
            break;
        case gates::gateNumber::u1:
            res =
                "\\gate{U1(" +
                trimZeroes(std::to_string(R3(euler_array[0]))) +
                ")} & ";
            length++;
            break;
        case gates::gateNumber::u2:
            res = ("\\gate{U2(" +
                   trimZeroes(
                       std::to_string(R3(euler_array[0]))) +
                   "," +
                   trimZeroes(
                       std::to_string(R3(euler_array[1]))) +
                   ")} & ");
            length += 2;
            break;

        case gates::gateNumber::rx:
            res =
                "\\gate{Rx(" +
                trimZeroes(std::to_string(R3(euler_array[0]))) +
                ")} & ";
            length++;
            break;
        case gates::gateNumber::ry:
            res =
                "\\gate{Ry(" +
                trimZeroes(std::to_string(R3(euler_array[0]))) +
                ")} & ";
            length++;
            break;
    }
    // std::cout<<" returns: "<<res<<std::endl;
    return {res, length};
}
void Qcircuit::draw_circuit_(int num ) {
    /*
        this function is responsible for drawing the circuit using Qcircuit
       typesetting langauge (latex)

    */
    // print_vector_();
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
    for (std::vector<std::vector<int>>::size_type i = 0; i < (*order_).size();
         i++) {
        // std::cout<<"current: ";
        // for(auto q : (*order_)[i]) std::cout<< q <<" ";
        // std::cout<<std::endl;
        if ((*order_)[i][0] == Operation::singleGateOPeration ) {
            qw[(*order_)[i][1]].pos++;
            if (qw[(*order_)[i][1]].pos > maxt) {
                maxt = qw[(*order_)[i][1]].pos;
            }
            auto tmp = draw_circuit_utility_((*order_)[i][2], (*order_)[i][3]);
            qw[(*order_)[i][1]].line.push_back({tmp.first, tmp.second});

        }
        // TODO OPTIMISE DRAWCIRCUIT
        else if ((*order_)[i][0] == Operation::twoControlledGateOperation) {
            // (*order_).push_back ({ 4, t_qbit, G, c_or_q1, control_bit1,
            // control_bit2, esize});
            int tbit = (*order_)[i][1];
            int gate = (*order_)[i][2];
            int cq = (*order_)[i][3]; //will use it in future

            int c1 = (*order_)[i][4];
            int c2 = (*order_)[i][5];
             int eindex = (*order_)[i][6]; //will use in future
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
            auto tmp = draw_circuit_utility_(gate, 0);
            qw[tbit].line.push_back(
                {tmp.first,
                 tmp.second});  // Append the result of draw_circuit_utility_
            maxt++;

        } else if ((*order_)[i][0] == Operation::controlledGateOperation ) {
            //  (*order_).push_back ({ 2, t_qbit, GATE, c_or_q, control_bit});
            if ((*order_)[i][3] == 1)  // type: q
            {
                // std::cout<<"control: G:  "<<(*order_)[i][2] << " from
                // "<<(*order_)[i][4]<< " to "<< (*order_)[i][1]  <<" \n";
                int low = std::min((*order_)[i][4], (*order_)[i][1]);
                int high = (low == (*order_)[i][4]) ? (*order_)[i][1] : (*order_)[i][4];
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

                // Add control gate for (*order_)[i][4]
                qw[(*order_)[i][4]].line.push_back({
                    "\\ctrl{" + std::to_string((*order_)[i][1] - (*order_)[i][4]) +
                        "} & ",
                    1  // Length increment
                });
                qw[(*order_)[i][4]].pos++;

                // Add the gate for (*order_)[i][1]
                qw[(*order_)[i][1]].pos++;
                auto tmp = draw_circuit_utility_((*order_)[i][2], 0);
                qw[(*order_)[i][1]].line.push_back(
                    {tmp.first,
                     tmp.second});  // Append the result of draw_circuit_utility_
                maxt++;

            } else if ((*order_)[i][3] == 0) {
                // (*order_).push_back ({ 2, t_qbit, G, c_or_q, control_bit,
                // esize});

                // [control, target] - add wires till the max
                for (int x = (*order_)[i][1]; x < n_qbits; x++) {
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
                for (int x = (*order_)[i][1] + 1; x < n_qbits; x++) {
                    qw[x].line.push_back(
                        {"\\qw \\cwx & ", 1});  // Append "\\qw \\cwx & " with a
                                                // length increment of 1
                    qw[x].pos++;
                }

                // The control and target get specific control/wire commands
                for (int x = 0; x < n_cbits; x++) {
                    if (x != (*order_)[i][4]) {
                        if (x < (*order_)[i][4]) {
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

                qw[(*order_)[i][1]].pos++;
                auto tmp = draw_circuit_utility_((*order_)[i][2], (*order_)[i][5]);
                qw[(*order_)[i][1]].line.push_back(
                    {tmp.first,
                     tmp.second});  // Append the result of draw_circuit_utility_
                maxt++;
            }
        } else if ((*order_)[i][0] == Operation::measureOperation)  // measure
        {
            // For the main, put wires
            while (qw[(*order_)[i][1]].pos < maxt) {
                qw[(*order_)[i][1]].line.push_back(
                    {"\\qw & ",
                     1});  // Append "\\qw & " with a length increment of 1
                qw[(*order_)[i][1]].pos++;
            }
            qw[(*order_)[i][1]].line.push_back(
                {"\\meter & ",
                 1});  // Append "\\meter & " with a length increment of 1
            qw[(*order_)[i][1]].pos++;

            // Q: Putting wires in between
            for (int x = (*order_)[i][1] + 1; x < n_qbits; x++) {
                while (qw[x].pos < maxt) {
                    qw[x].line.push_back(
                        {"\\qw & ",
                         1});  // Append "\\qw & " with a length increment of 1
                    qw[x].pos++;
                }
            }

            // Q: Edgy wires for the between
            for (int x = (*order_)[i][1] + 1; x < n_qbits; x++) {
                qw[x].line.push_back(
                    {"\\qw \\cwx & ", 1});  // Append "\\qw \\cwx & " with a
                                            // length increment of 1
                qw[x].pos++;
            }

            // C: Putting cwires
            for (int x = 0; x <= (*order_)[i][2]; ++x) {
                while (cw[x].pos < maxt) {
                    cw[x].line.push_back(
                        {"\\cw & ",
                         1});  // Append "\\cw & " with a length increment of 1
                    cw[x].pos++;
                }
            }

            // C: Putting edgy wires
            for (int x = 0; x <= (*order_)[i][2]; ++x) {
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

handler->latex_writer << "\\clearpage\n\\section*{Quantum Circuit}\n";
handler->latex_writer << "\\begin{figure}[htbp]\n";
handler->latex_writer << "    \\centering\n";
handler->latex_writer << "    \\[\n";
handler->latex_writer << "    \\Qcircuit @C=1em @R=.7em {\n";

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

    for (size_t z = starter; z < qw[0].line.size(); z++) {
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
    int min_diff = 0xffffffff;
    for (size_t i = 0; i < inds.size() - 1; i++) {
        min_diff = (min_diff > (inds[i + 1] - inds[i]))
                       ? ((inds[i + 1] - inds[i]))
                       : (min_diff);
    }
    //int tmp = inds.size();
    inds.clear();
    for (size_t i = 0; i <= qw[0].line.size(); i += min_diff) {
        inds.push_back(i);
        std::cout << i << " ";
    }
    std::cout << std::endl;

    for (size_t i = 0; i < inds.size() - 1; i++)  // i about iterations
    {
        for (int j = 0; j < n_qbits; j++)  // j about q/c number
        {
            std::string combinedLine = qw[j].line[0].first;
            for (int k = inds[i]; k < inds[i + 1];
                 k++)  // k about line indicies
                if (k) combinedLine += qw[j].line[k].first;

            handler->latex_writer << combinedLine << " \\qw & \\\\" << std::endl;
        }
        for (int j = 0; j < n_cbits; j++) {
            std::string combinedLine = cw[j].line[0].first;
            for (int k = inds[i]; k < inds[i + 1]; k++)
                if (k) combinedLine += cw[j].line[k].first;

            handler->latex_writer << combinedLine << " \\cw & \\\\" << std::endl;
        }
        handler->latex_writer << "\\\\ \n\\\\ \n";
    }
    for (int i = 0; i < n_qbits; i++) {
        std::string combinedLine = qw[i].line[0].first;
        for (size_t j = inds[inds.size() - 1]; j < qw[0].line.size(); j++)
            if (j) combinedLine += qw[i].line[j].first;

        handler->latex_writer << combinedLine << std::endl;
    }
    for (int i = 0; i < n_cbits; i++) {
        std::string combinedLine = cw[i].line[0].first;
        for (size_t j = inds[inds.size() - 1]; j < qw[0].line.size(); j++)
            if (j) combinedLine += cw[i].line[j].first;

        handler->latex_writer << combinedLine << std::endl;
    }
    handler->latex_writer << "\\\\ \n\\\\ \n";

    // std::cout << "}\n\\end{displaymath}\n";
    if (handler->latex_writer.good()) {
        handler->latex_writer << "}\n\\]\n";
        if (name != "") handler->latex_writer << "\\caption{" << name << "}\n";
        handler->latex_writer << "\\end{figure}\n";
        handler->latex_writer.flush();
    } else
        std::cerr << "writing issue";

    // std::cout<<"CIRCUIT maxt: "<<maxt << " len: "<< max_length << std::endl;
    // handler->write_in_pdf("CIRCUIT maxt: "+ std::to_string(maxt)+ " len: " +
    // std::to_string(max_length));
}
