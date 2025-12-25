// File       : main.cpp
// Description: CahnHilliard application code
// Copyright 2023 Harvard University. All Rights Reserved.
#include "CahnHilliard.h"
#include <mpi.h>

int main(int argc, char *argv[])
{
    // TODO: make sure you properly support multiple threads
    MPI_Init(&argc, &argv);

    const int N = 350; 
    CahnHilliard *sim = new CahnHilliard(N);
    sim->integrate();
    sim->dump("final.bin");
    
    delete sim;

    MPI_Finalize();
    return 0;
}
