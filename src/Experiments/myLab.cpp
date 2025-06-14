
#include "myLab.h"
#include <array>


Experiment quantum_teleportation_prtcl_var1(MimiqHandler* handler) 
{
    Qcircuit qc(handler, 3,3, "quantum teleportation 1"); 

    qc.u(0, {0.3,0.2,0.1});
    qc.h(1);
    qc.cx(1,2);
    qc.cx(0,1);
    qc.h(0);  
    qc.measure(0,0);
    qc.measure(1,1);
    qc.apply_controlled_gate("x","c1",2);
    qc.apply_controlled_gate("z","c0",2);
    qc.measure(2,2);
    
    return qc.simulation();
}

Experiment quantum_teleportation_prtcl_var2(MimiqHandler* handler)
{
    Qcircuit qc(handler, 3,3, "quantum teleportation ibm"); 
    
    qc.u(0, {6.28,0.408,2.22});
    qc.h(1);
    qc.cx(1,2);
    qc.cx(0,1);
    qc.h(0);  
    qc.measure(0,0);
    qc.measure(1,1);
    qc.apply_controlled_gate("x","c1",2);
    qc.apply_controlled_gate("z","c0",2);
    qc.u(2, {-6.28,-2.22,-0.408});
    qc.measure(2,2);
    
    return qc.simulation();   
}

Experiment quantum_superdense_random(MimiqHandler* handler)
{
    Qcircuit qc(handler, 3,4, "superdense coding - random");
    qc.h(2);
    qc.measure(2,0);
    qc.apply_controlled_gate("x","c0",2);
    qc.h(2);
    qc.measure(2,1);
    
    qc.h(0);
    qc.cx(0,1);
    qc.apply_controlled_gate("z","c1",0); 
    qc.apply_controlled_gate("x","c0",0);
    qc.cx(0,1);
    qc.h(0);
    
    qc.measure(0,3);
    qc.measure(1,2);
    return qc.simulation();
}

Experiment quantum_superdense(MimiqHandler* handler)
{
    Qcircuit qc(handler, 2,4, "superdense coding"); 

    qc.x(0);
    qc.measure(0,1);
    qc.x(0);
    qc.h(0);
    qc.cx(0,1);
    qc.apply_controlled_gate("z","c1",0); 
    qc.apply_controlled_gate("x","c0",0);
    qc.cx(0,1);
    qc.h(0);
    qc.measure(0,3);
    qc.measure(1,2);
    
    return qc.simulation();
}

Experiment Grover(MimiqHandler* handler)
{
    Qcircuit qc(handler,3,3,"grover");
    
    qc.h(2);
    qc.h(1);
    qc.h(0);

    qc.x(0);
    qc.h(2);
    qc.toffoli(0,1,2);
    qc.h(2);
    qc.x(0);

    qc.h(0);
    qc.h(1);
    qc.h(2);
    qc.x(0);
    qc.x(1);
    qc.x(2);
    qc.h(2);
    qc.toffoli(0,1,2);
    qc.h(2);
    qc.x(1);
    qc.x(0);
    qc.x(2);
    qc.h(0);
    qc.h(1);
    qc.h(2);

    qc.x(0);
    qc.h(2);
    qc.toffoli(0,1,2);
    qc.h(2);
    qc.x(0);

    qc.h(0);
    qc.h(1);
    qc.h(2);
    qc.x(0);
    qc.x(1);
    qc.x(2);
    qc.h(2);
    qc.toffoli(0,1,2);
    qc.h(2);
    qc.x(1);
    qc.x(0);
    qc.x(2);
    qc.h(0);
    qc.h(1);
    qc.h(2);

    qc.measure(0,0);
    qc.measure(1,1);
    qc.measure(2,2);

    return qc.simulation();
}

Experiment Grover2(MimiqHandler* handler)
{
    Qcircuit qc(handler,3,3,"grover by circuit");
    
    qc.h(0);
    qc.h(1);
    qc.h(2); 
    qc.x(0);

    qc.h(2); 
    qc.toffoli(0,1,2);
    
    qc.x(0);
    qc.h(2); 
   
    qc.h(0);
    qc.h(1);
    qc.h(2); 
    qc.x(1);
    qc.x(0);
    qc.x(2);  

    qc.h(2); 
    qc.toffoli(0,1,2);

    qc.h(2); 

    qc.x(1);
    qc.x(0);
    qc.x(2);     
    qc.h(0);
    qc.h(1);
    qc.h(2);

    qc.x(0);

    qc.h(2); 
    qc.toffoli(0,1,2);

    qc.x(0);
    qc.h(2); 
 
    qc.h(0);
    qc.h(1);
    qc.h(2);
    qc.x(1);
    qc.x(0);
    qc.x(2);    
    qc.h(2); 
    qc.toffoli(0,1,2);
    qc.h(2);
    qc.x(1);
    qc.x(0);
    qc.x(2);     
    qc.h(0);
    qc.h(1);
    qc.h(2); 

    qc.measure(0,0);
    qc.measure(1,1);
    qc.measure(2,2);

    return qc.simulation();
}

Experiment BB84_QKD(MimiqHandler* handler)
{
    Qcircuit qc(handler,2,6, "BB84-ptcl QKD");
    // alice

    qc.h(1);                                                // init randomiser
    qc.measure(1,0);                                         // c0 has a random -> alice bit 
    qc.apply_controlled_gate("x","c0",1);                                 // reset randomiser
    qc.apply_controlled_gate("x","c0",0);                                 // alice bit applied
    qc.h(1);                                                // init randomiser
    qc.measure(1,1);                                         // c1 has a random -> alice basis
    qc.apply_controlled_gate("x","c1",1);                                 // reset randomiser
    qc.apply_controlled_gate("h","c1",0);                                 // alice basis applied
    // eve 
    qc.h(1);                                                // init randomiser
    qc.measure(1,4);                                         // c4 has a random value -> eve presence
    qc.apply_controlled_gate("x","c4",1);                                 // reset randomiser
    if( qc.accessCreg(4) )
    {
        qc.h(1);                                            // init randomiser
        qc.measure(1,5);                                     // c5 has random value -> eve bit
        qc.apply_controlled_gate("x","c5",1);                             // reset randomiser
        qc.apply_controlled_gate("x","c5",0);                             // eve bit applied
        qc.h(1);                                            // init randomiser
        qc.measure(1,5);                                     // c5 has random value -> eve basis
        qc.apply_controlled_gate("x","c5",1);                             // reset randomiser
        qc.apply_controlled_gate("h","c5",0);                             // eve basis applied
    }
    // bob 
    qc.h(1);                                                // init randomiser
    qc.measure(1,2);                                         // c2 has randomvalue -> bob basis
    qc.measure(0,3,0,qc.accessCreg(2)); // c3 has bob's measurement -> bob bit
    return qc.simulation();
} // qc.ifcreg(0).hx

Experiment Quantum_full_adder(MimiqHandler* handler)
{
    Qcircuit qc(handler, 4,4, "full adder");
    /*
        q0 -> A 
        q1 -> B
        q2 -> carryin_sumout
        q3 -> carryout
    */    

    // initialize inputs to some values
    // initialize inputs A=1, B=0 and carry_in=1
    qc.x(0);
    qc.x(2);

    // perform addition
    // decomposition of toffoli q[0], q[1], q[2]
    qc.toffoli(0,1,2);
    /*
    qc.h(2);
    qc.cx(1,2);
    qc.apply_gate_("td",2);
    qc.cx(0,2);
    qc.apply_gate_("t",2);
    qc.cx(1,2);
    qc.apply_gate_("td",2);
    qc.cx(0,2);
    qc.apply_gate_("t",1);
    qc.apply_gate_("t",2);
    qc.cx(0,1);
    qc.h(2);
    qc.apply_gate_("t",0);
    qc.apply_gate_("td",1);
    qc.cx(0,1);
    */

    qc.cx(0,1);

    // decomposition of toffoli q[1], q[2], q[3]
    qc.toffoli(1,2,3);
    /*
    qc.h(3);
    qc.apply_controlled_gate("x","q2",3);
    qc.apply_gate_("td",3);
    qc.cx(1,3);
    qc.apply_gate_("t",3);
    qc.apply_controlled_gate("x","q2",3);
    qc.apply_gate_("td",3);
    qc.cx(1,3);
    qc.apply_gate_("t",2);
    qc.apply_gate_("t",3);
    qc.cx(1,2);
    qc.h(3);
    qc.apply_gate_("t",1);
    qc.apply_gate_("td",2);
    qc.cx(1,2); */

    qc.cx(1,2);
    qc.cx(0,1);

    qc.measure(0,0);
    qc.measure(1,1);
    qc.measure(2,2);
    qc.measure(3,3);

    return qc.simulation();
}

Experiment trial(MimiqHandler* handler)
{
    Qcircuit qc(handler,3,3,"trial");

    qc.h(0);
    qc.x(1);
    qc.toffoli(0,2,1);
    qc.measure(1,2);

    return qc.simulation();
}

Experiment quantum_classification(MimiqHandler* handler)
{
    Qcircuit qc(handler,4,4,"quantum classification");
    qc.h(0);
    qc.h(1);

    /* Encode new data point */
    // encode the new data point, this implements a cRy(omega)
    qc.cx(1,2);
    qc.ry(2,{0.1105});// (Ry q[2], -omega/2)
    qc.cx(1,2);
    qc.ry(2,{-0.1105});// (Ry q[2], omega/2)
    qc.x(1);


    /* Encode first training point */
    // encode the first data point, this implements a ccRy(theta)

    // decomposition of toffoli q[0], q[1], q[2]
    qc.toffoli(0,1,2);

    qc.cx(0,2);
    qc.ry(2,{0});// (Ry q[2], theta/4)
    qc.cx(0,2);
    qc.ry(2,{0});// (Ry q[2], -theta/4)

    // decomposition of toffoli q[0], q[1], q[2]
    qc.toffoli(0,1,2);
    qc.cx(0,2);
    qc.ry(2,{0});// (Ry q[2], theta/4)
    qc.cx(0,2);
    qc.ry(2,{0});// (Ry q[2], -theta/4)
    qc.x(0);

    /* .Encode second training point */
    // encode the second data point, this implements a ccRy(phi)

    // decomposition of toffoli q[0], q[1], q[2]
    qc.toffoli(0,1,2);
    qc.cx(0,2);
    qc.ry(2,{-1.511125});// (Ry q[2], phi/4)
    qc.cx(0,2);
    qc.ry(2,{1.511125});// (Ry q[2], -phi/4)

    // decomposition of toffoli q[0], q[1], q[2]
    qc.toffoli(0,1,2);

    qc.cx(0,2);
    qc.ry(2,{1.511125});// (Ry q[2], -phi/4)
    qc.cx(0,2);
    qc.ry(2,{-1.511125});// (Ry q[2], phi/4)

    /* Labels */
    // encode the labels
    qc.cx(0,3);

    /* Algorithm */
    // The actual algorithm
    qc.h(1);

    qc.measure(1,1);
    qc.measure(3,3);
    return qc.simulation();
}
