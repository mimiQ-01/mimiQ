/**
 * @file main.cpp
 * @brief main driver file
 * 
 *  
 * @author Rushikesh Muraharisetty
 * @date last updated: Mar 5 '25
 * 
 */

typedef void Lab; // semantic
#include "mimiq/mimiq.h"
#include "Experiments/myLab.h"
#include <iostream>

Lab myLab(MimiqHandler* handler)
{
    result res; // a result is needed to hold a simulationn result
    /*
        Quantum teleportation - variation 1 - qiskit example 
    */
    res = simulate(handler,quantum_teleportation_prtcl_var1, 1024 );
    res.print_counts();
    res.generate_openqasm();
    handler->clean();
    /*
        Quantum teleportation - variation 2 - IBM quantum learning example
    */
    res = simulate(handler,quantum_teleportation_prtcl_var2, 1024 );
    res.print_counts();
    res.generate_openqasm();
    handler->clean();
    /*
        Quantum superdense coding - IBM quantum learning 1st example 
    */
    res = simulate(handler,quantum_superdense, 1024 );
    res.print_counts();
    res.generate_openqasm();
    handler->clean();
    /*
        Quantum superdense coding - IBM quantum learning 2nd example 
    */
    res = simulate(handler,quantum_superdense_random, 1024 );
    res.print_counts();
    res.generate_openqasm();
    handler->clean();

    /*
        Quantum key distribution - BB98 protocol - own implementation
    */
    int common_basis = 0 ,leaked_bits = 0, eve_attacks = 0 , checker=0;
    std::string key = "";
    int keySize = 50 ;

    while(keySize)
    {
        res = simulate(handler, BB84_QKD);
        //handler->clean();
        if(res.creg[1] == res.creg[2])
        {
            common_basis++;

            if(res.creg[0] != res.creg[3]) 
                leaked_bits++;

            else 
            {
                key += (res.creg[0] ? "1":"0" );
                keySize--;
            }

            if(res.creg[4])
            {
                checker++;
                eve_attacks++;
            }
            // ASSERT(checker == leaked_bits)
        }
    }
    std::cout<< "common bases: "<< common_basis << " eveattacks: "<< eve_attacks << " detected eve attcks: " << leaked_bits << std::endl;
    std::cout<<"\nFINAL KEY: "<< key <<" length: "<< key.length()<< std::endl;
    handler->write_in_pdf("Resultant key: " + key );
    handler->clean();
    
    /*
        Grovers search (quantum-inspire example)
    */
    res = simulate(handler,Grover, 1000 );
    res.print_counts();
    res.generate_openqasm();
    handler->clean();
    /*
        quantum full adder (quantum-inspire example)
    */
    res = simulate(handler,Quantum_full_adder, 1000 );
    res.print_counts();
    res.generate_openqasm();
    handler->clean();
    /*
        quantum classification (quantum-inspire example)
    */
    res = simulate(handler,quantum_classification, 1000 );
    res.print_counts();
    res.generate_openqasm();
    handler->clean();
    
}
Lab lab2(MimiqHandler* handler)
{
    result res;
    res = simulate(handler,trial, 10 );
    res.print_counts();
    res.generate_openqasm();
    //res.generate_openqasm2();
    
}
int main()
{  
    MimiqHandler* myhandler= new MimiqHandler("/home/zenyuzien/Downloads/"); // path to store report, default is current directory
    lab2(myhandler);
    myhandler->generate_report();
    delete myhandler;
    return 0 ;
}

// tasks for future:
/*
1. fix state_.print()
5. parallelism
6. cqasm code gen
7. openqasm code gen
8. add more features from net
*/
