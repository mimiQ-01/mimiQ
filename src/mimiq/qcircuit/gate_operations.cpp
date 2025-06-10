#include "../qcircuit.hpp"
#include "../gates.hpp"
#include <cmath>
#include <vector>
#include <array>
#include <sstream>

#define PI acos(0.0) * 2
static std::vector<std::vector<Coeff> > empty_matrix ; // for applyUtility

void Qcircuit::apply_utility_(
    const std::string& gate, int& G, int& esize,
    std::vector<double> arr = std::vector<double>(),
    std::vector<std::vector<Coeff>>& matrix = empty_matrix) {
    if (gate == "H" || gate == "h") {
        G =  gates::gateNumber::h ;
        matrix = gates::unitaryGate(PI / 2, 0, PI);
    } else if (gate == "x" || gate == "X") {
        G = gates::gateNumber::x;
        matrix = gates::unitaryGate(PI, 0, PI);
    } else if (gate == "y" || gate == "Y") {
        G = gates::gateNumber::y;
        matrix = gates::unitaryGate(PI, PI / 2, PI / 2);
    } else if (gate == "z" || gate == "Z") {
        G = gates::gateNumber::z;
        matrix = gates::unitaryGate(0, 0, PI);
    } else if (gate == "s" || gate == "S") {
        G = gates::gateNumber::s;
        matrix = gates::unitaryGate(0, 0, PI / 2);
    } else if (gate == "sd" || gate == "Sd") {
        G = gates::gateNumber::sd;
        matrix = gates::unitaryGate(0, 0, -PI / 2);
    } else if (gate == "T" || gate == "t") {
        G = gates::gateNumber::t;
        matrix = gates::unitaryGate(0, 0, PI / 4);
    } else if (gate == "Td" || gate == "td") {
        G = gates::gateNumber::td;
        matrix = gates::unitaryGate(0, 0, -PI / 4);
    } else if (gate == "rx" || gate == "Rx") {
        G = gates::gateNumber::rx;
        esize = (*euler_container_).size();
        (*euler_container_).push_back({arr[0], 0, 0});

        matrix = gates::unitaryGate((*euler_container_)[esize][0], -PI / 2, PI / 2);
    } else if (gate == "ry" || gate == "Ry") {
        G = gates::gateNumber::ry;
        esize = (*euler_container_).size();
        (*euler_container_).push_back({arr[0], 0, 0});

        matrix = gates::unitaryGate((*euler_container_)[esize][0], 0, 0);
    } else if (gate == "U3" || gate == "u3" || gate == "u" || gate == "U") {
        esize = (*euler_container_).size();
        (*euler_container_).push_back({arr[0], arr[1], arr[2]});
        G = gates::gateNumber::u3;
        matrix = gates::unitaryGate((*euler_container_)[esize][0],
                                    (*euler_container_)[esize][1],
                                    (*euler_container_)[esize][2]);
    } else if (gate == "IU3" || gate == "iu3" || gate == "IU" ||
               gate == "iu")  // todo: generalise u1 u2 u3 inverse
    {
        esize = (*euler_container_).size();
        (*euler_container_).push_back({arr[0], arr[1], arr[2]});
        G = gates::gateNumber::iu3;
        matrix = gates::inverse2x2(gates::unitaryGate(
            (*euler_container_)[esize][0], (*euler_container_)[esize][1],
            (*euler_container_)[esize][2]));
    } else if (gate == "U1" || gate == "u1") {
        esize = (*euler_container_).size();
        (*euler_container_).push_back({arr[0], 0, 0});
        G = gates::gateNumber::u1;
        matrix = gates::unitaryGate(0, 0, (*euler_container_)[esize][0]);
    } else if (gate == "U2" || gate == "u2") {
        esize = (*euler_container_).size();
        (*euler_container_).push_back({arr[0], arr[1], 0});
        G = gates::gateNumber::u2;
        matrix = gates::unitaryGate(PI / 2, (*euler_container_)[esize][0],
                                    (*euler_container_)[esize][1]);
    }
}

void Qcircuit::apply_gate_(const std::string& gate, int t_qbit,
                         std::vector<double> arr, bool ORDERinc) {
    std::vector<std::vector<Coeff>> matrix, buffedmatrix;
    int esize = -1, G;

    apply_utility_(gate, G, esize, arr, matrix);

    // std::cout<<"p: "<<1 <<" "<<t_qbit << " " << G << " " << esize <<
    // std::endl;
    if (ORDERinc) (*order_).push_back({ Operation::singleGateOPeration , t_qbit, G, esize});  // send index

    gates::kronecker_buff(matrix, buffedmatrix, n_qbits, t_qbit);
    int n = 1 << n_qbits;
    std::vector<Coeff> result(n, {0, 0});
    for (int ii = 0; ii < n; ++ii) {
        for (int j = 0; j < n; ++j)
            result[ii] = result[ii] + (buffedmatrix[ii][j] * (*state_).coeffs[j]);
    }
    (*state_).coeffs = result;
}
void Qcircuit::apply_controlled_gate(const std::string& gate,
                                   const std::string c_qbit_str, int t_qbit,
                                   std::vector<double> arr) {
    std::istringstream iss(c_qbit_str);
    char cq1, cq2, d = '0';
    int control_bit1, control_bit2, c_or_q1, c_or_q2, G;
    iss >> cq1 >> control_bit1 >> d >> cq2 >> control_bit2;

    if (d == '0') { // WHEN DELIMETER NOT RECEIVED, MEANING SINGLE CONTROL, not toffoli
        c_or_q1 = (cq1 == 'q') ? 1 : 0;

        if (c_or_q1 == 0) {
            
            if(handler->qasm_gen) handler->oqsm +=             
            "if ( c == " + std::to_string(1 << control_bit1) + ") " + gate +  " q[" + std::to_string(t_qbit) + "];\n";
            // todo gate validate the user written gate 


            if (c_reg & (1 << (n_cbits - 1 - control_bit1)))
                apply_gate_(gate, t_qbit, arr, false);
        }

        int esize;
        apply_utility_(gate, G, esize, arr);

        (*order_).push_back({ Operation::controlledGateOperation , t_qbit, G, c_or_q1, control_bit1, esize});

        if (c_or_q1 == 1) {  // CONTROL ON ZERO WITH Q needs to be done. TODO
            if (G == 1)      // gate
                gates::controlled_hadamard((*state_), control_bit1,
                                           t_qbit);  // qbit, tqbit

            else if (G == 2)
                gates::controlled_pauli_X((*state_), control_bit1, t_qbit);

            else if (G == 4)  // 4 - pauli z , 3 - pauli y
                gates::controlled_pauli_Z((*state_), control_bit1, t_qbit);
        }
    } else if (d == ',') {
        c_or_q1 = (cq1 == 'q') ? 1 : 0;
        c_or_q2 = (cq2 == 'q') ? 1 : 0;
        if ((c_or_q1 == c_or_q2) && (c_or_q1 == 0)) {
            if ((c_reg & (1 << (n_cbits - 1 - control_bit1))) &&
                (c_reg & (1 << (n_cbits - 1 - control_bit2))))
                apply_gate_(gate, t_qbit, arr, false);
        }
        int esize;
        apply_utility_(gate, G, esize, arr);
        // circuit
        (*order_).push_back(
            { Operation::twoControlledGateOperation , t_qbit, G, c_or_q1, control_bit1, control_bit2, esize});

        if ((c_or_q1 == c_or_q2) &&
            (c_or_q1 == 1)) {  // CONTROL ON ZERO WITH Q needs to be done. TODO
            if (G == 2)
                gates::toffoli((*state_), control_bit1, control_bit2, t_qbit);
        }
    }
}



    void Qcircuit::u( int t_qbit, std::vector<double> arr, bool ORDERinc)
    {
        apply_gate_("u",t_qbit,arr,ORDERinc);
        if(handler->qasm_gen) handler->oqsm += "u("+ std::to_string(arr[0]) + ","  + std::to_string(arr[1]) + ","  + std::to_string(arr[2])  + ") q[" + std::to_string(t_qbit) + "];\n";
    }
    void Qcircuit::x( int t_qbit, std::vector<double> arr, bool ORDERinc)
    {
        apply_gate_("x",t_qbit,arr,ORDERinc);
        if(handler->qasm_gen) handler->oqsm += "x q[" + std::to_string(t_qbit) + "];\n";
    }
    void Qcircuit::y( int t_qbit, std::vector<double> arr, bool ORDERinc)
    {
        apply_gate_("y",t_qbit,arr,ORDERinc);
        if(handler->qasm_gen) handler->oqsm += "y q[" + std::to_string(t_qbit) + "];\n";
    }
    void Qcircuit::z( int t_qbit, std::vector<double> arr, bool ORDERinc)
    {
        apply_gate_("z",t_qbit,arr,ORDERinc);
        if(handler->qasm_gen) handler->oqsm += "z q[" + std::to_string(t_qbit) + "];\n";
    }
    void Qcircuit::h( int t_qbit, std::vector<double> arr, bool ORDERinc)
    {
        apply_gate_("h",t_qbit,arr,ORDERinc);
        if(handler->qasm_gen) handler->oqsm += "h q[" + std::to_string(t_qbit) + "];\n";
    }
    void Qcircuit::u2( int t_qbit, std::vector<double> arr, bool ORDERinc)
    {
        apply_gate_("u2",t_qbit,arr,ORDERinc);
        if(handler->qasm_gen) handler->oqsm += "u(pi/2,"+ std::to_string(arr[0]) + ","  + std::to_string(arr[1]) + ") q[" + std::to_string(t_qbit) + "];\n";
    }
    void Qcircuit::u1( int t_qbit, std::vector<double> arr, bool ORDERinc)
    {
        apply_gate_("u",t_qbit,arr,ORDERinc);
        if(handler->qasm_gen) handler->oqsm += "u(0,0,"+ std::to_string(arr[0])  + ") q[" + std::to_string(t_qbit) + "];\n";
    }
    void Qcircuit::t( int t_qbit, std::vector<double> arr, bool ORDERinc)
    {
        apply_gate_("t",t_qbit,arr,ORDERinc);
        if(handler->qasm_gen) handler->oqsm += "t q[" + std::to_string(t_qbit) + "];\n";
    }
    void Qcircuit::tdg( int t_qbit, std::vector<double> arr, bool ORDERinc)
    {
        apply_gate_("td",t_qbit,arr,ORDERinc);
        if(handler->qasm_gen) handler->oqsm += "tdg q[" + std::to_string(t_qbit) + "];\n";
    }
    void Qcircuit::s( int t_qbit, std::vector<double> arr, bool ORDERinc)
    {
        apply_gate_("s",t_qbit,arr,ORDERinc);
        if(handler->qasm_gen) handler->oqsm += "s q[" + std::to_string(t_qbit) + "];\n";
    }
    void Qcircuit::sdg( int t_qbit, std::vector<double> arr, bool ORDERinc )
    {
        apply_gate_("sd",t_qbit,arr,ORDERinc);
        if(handler->qasm_gen) handler->oqsm += "sdg q[" + std::to_string(t_qbit) + "];\n";
    }
    void Qcircuit::ry( int t_qbit, std::vector<double> arr, bool ORDERinc )
    {
        apply_gate_("ry",t_qbit,arr,ORDERinc);
        if(handler->qasm_gen) handler->oqsm += "ry("+ std::to_string(arr[0])  + ") q[" + std::to_string(t_qbit) + "];\n";
    }
    void Qcircuit::toffoli(int c1_qbit, int c2_qbit, int t_qbit)
    {
        apply_controlled_gate("x",("q" + std::to_string(c1_qbit) + ",q" + std::to_string(c2_qbit)),t_qbit );
        if(handler->qasm_gen) handler->oqsm += "ccx q[" + std::to_string(c1_qbit) + "], q[" + std::to_string(c2_qbit) + "], q["  + std::to_string(t_qbit) +  "];\n";
    }
    void Qcircuit::ccx(int c1_qbit, int c2_qbit, int t_qbit)
    {
        apply_controlled_gate("x",("q" + std::to_string(c1_qbit) + ",q" + std::to_string(c2_qbit)),t_qbit );
        if(handler->qasm_gen) handler->oqsm += "ccx q[" + std::to_string(c1_qbit) + "], q[" + std::to_string(c2_qbit) + "], q["  + std::to_string(t_qbit) +  "];\n";
    }
    void Qcircuit::cx(int cbit, int qbit)
    {
        apply_controlled_gate("x","q"+std::to_string(cbit), qbit);
        if(handler->qasm_gen) handler->oqsm += "cx q[" + std::to_string(cbit) + "], q[" + std::to_string(qbit) + "];\n";
    }
    void Qcircuit::ch(int cbit, int qbit)
    {
        apply_controlled_gate("h","q"+std::to_string(cbit), qbit);
        if(handler->qasm_gen) handler->oqsm += "ch q[" + std::to_string(cbit) + "], q[" + std::to_string(qbit) + "];\n";
    }
