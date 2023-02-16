void defineRxnKinetics(std::vector<double> B, std::vector<double> C, std::vector<double> &subsR, std::vector<double> &bioR, plint mask )
{
    // B = biomass density, C = solute concentration, subsR = substrate reaction rate, bioR = microbial reaction rate, mask = local mask number

    // This function is used only when solver_type "kinetics" is used in CompLaB.xml file.
    // It is up to users to define the reaction rate (R) and stoichiometry vectors (S) in this function.
    // The length of B and C vectors should be equal to num_of_microbes and total_num_of_subs defined in CompLaB.xml file.
    // Indexing starts from 0 following C/C++ convention.

    // Example 1 (Lines 10-24)
    // The reaction rate expressions
    // subsR0 = dC1/dt = -kc1*[B1]*[C1]*[C2]
    // subsR1 = dC2/dt = subsR0
    // subsR2 = dC3/dt = -subsR0-subsR3
    // subsR3 = dC4/dt = kc2*[B2]*[C3]/([C2] + K2)
    // bioR0 = dB1/dt = Y1*subsR0-kd1*[B1] // microbial growth yield (Y1) and death rate constant (kd1)
    // bioR1 = dB2/dt = Y2*subsR3-kd2*[B2] // microbial growth yield (Y2) and death rate constant (kd2)

        
   
       //subsR[0] = 1.722e-4 * C[0]; // 
        







    //R[0] = -0.8135 * B[0] * C[0] * C[1]; // where kc1 = 0.8135
    // R[1] = R[0];
    // R[3] = 1.333 * B[1] * C[2] / (C[1] + 3.5); // where kc2 = 1.333, K2 = 3.5
    // R[2] = -R[0]-R[3];
    // R[4] = 2.1*R[0]-0.004*B[0]; // where Y1 = 2.1, kd1 = 0.0004
    // R[5] = 1.05*R[3]-0.055*B[1]; // where Y2 = 1.05, kd2 = 0.055
   
    double kads = 2.1e-3; // 1/s
    double kdes = 1.2e-5; // 1/s
    double mu1 = 1.74e-6; // 1/s
    double k_nh4 = 1.722e-4; // mM/s (Burdige, 1989, Biogeochemistry, The effects of sediment slurrying on microbial processes, and the role of amino acids as substrates for sulfate reduction in anoxic marine sediments)
    // double k_nh4 = 5.903e-6; // mM/s (Klump and Martens, 1989, Limnol. Oceanogr.)
    /* if (mask == 3) {
        bioR[0] = kads*B[1]-kdes*B[0];
        bioR[1] = -kads*B[1]+kdes*B[0]-mu1*B[1];
        subsR[1] = k_nh4;
    }
    else if (mask == 0) {
        subsR[1] = k_nh4;
        bioR[1] = -mu1*B[1];
    }
    else {
        bioR[1] = -mu1*B[1];
    }
    */

    subsR[1] = k_nh4;
    bioR[1] = -mu1 * B[1];

}
