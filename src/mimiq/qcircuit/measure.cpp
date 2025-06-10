#include "../qcircuit.hpp"
#include "../gates.hpp"
#include <cmath>
#include <iostream>

void Qcircuit::measure(int tbit,
                       int cbit,  // cbit = classical // TODO: int->false
                       int toPrint,
                       int basis)  // 0 basis  = along 0 or 1 , 1 bais = along +
                                   // or - , 2 basis = along +i or -i
{
    if(handler->qasm_gen) handler->oqsm+="measure q["+ std::to_string(tbit) +"] -> c[" +std::to_string(cbit) +"];\n" ;
    std::pair<Coeff, Coeff> res1;
    if (basis == 1) {
        std::vector<double> placebo;
        apply_gate_("h", tbit, placebo, false);
        res1 = (*state_).measure_along(tbit);
        apply_gate_("h", tbit, placebo, false);
    } else
        res1 = (*state_).measure_along(tbit);

    if (describe_ || toPrint) {
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
    if (toPrint == 1 || describe_) {
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

    // now update state_ vector to cause collapse
    // TODO: the search is exhaustive, u only need to check limited. Find a
    // better way if possible
    Coeff normaliser = {0, 0};
    if (res)  // if res== 1 , take the zero value and add it to 1 value
    {
        if (describe_) std::cout << "measured 1 so.. tmask " << tmask << " \n";
        for (int i = 0; i < cnd; ++i) {
            if (i & (tmask))  // we got x1x
            {
                if (describe_)
                    std::cout << "for i =  " << i << "adding to "
                              << (i | (tmask)) << " from " << (i & (~tmask))
                              << std::endl;

                normaliser.real =
                    normaliser.real + (*state_).coeffs[i | tmask].amp_sq();
                (*state_).coeffs[i & (~tmask)] = {0, 0};
            }
        }
    } else  // res == 0
    {
        if (describe_) std::cout << "measured 0 so.. tmask is " << tmask << "\n";
        for (int i = 0; i < cnd; ++i) {
            if (i & (tmask))  // we got x1x
            {
                if (describe_)
                    std::cout << "for i =  " << i << " adding to "
                              << (i & (~tmask)) << " from " << (i | (tmask))
                              << std::endl;

                normaliser.real =
                    normaliser.real + (*state_).coeffs[i & (~tmask)].amp_sq();
                (*state_).coeffs[i | tmask] = {0, 0};
            }
        }
    }
    normaliser.real = std::sqrt(normaliser.real);
    if (res)  // if res== 1 , take the zero value and add it to 1 value
    {
        if (describe_) std::cout << "measured 1 so.. tmask " << tmask << " \n";
        for (int i = 0; i < cnd; ++i) {
            if (i & (tmask))  // we got x1x
                (*state_).coeffs[i | tmask] = (*state_).coeffs[i | tmask] / normaliser;
        }
    } else  // res == 0
    {
        if (describe_) std::cout << "measured 0 so.. tmask is " << tmask << "\n";
        for (int i = 0; i < cnd; ++i) {
            if (i & (tmask))  // we got x1x
                (*state_).coeffs[i & (~tmask)] =
                    (*state_).coeffs[i & (~tmask)] / normaliser;
        }
    }
    if (describe_) std::cout << "updated creg: " << c_reg << std::endl;

    // std::cout<<"m "<< tbit << " "<<cbit << std::endl;
    (*order_).push_back({ Operation::measureOperation , tbit, cbit, toPrint, basis});
}