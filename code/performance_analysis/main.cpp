// File       : main.cpp
// Description: CahnHilliard application code
// Copyright 2023 Harvard University. All Rights Reserved.
#include "papi.h"
#include "CahnHilliard.h"
#include <mpi.h>
#include <cstdlib>

int main(int argc, char *argv[])
{   
    // check correct input of grid size from command line
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <N>" << std::endl;
        return EXIT_FAILURE;
    }
    
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    if (provided < MPI_THREAD_FUNNELED) {
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    }

    const int N = std::atoi(argv[1]);

    // int event_set = PAPI_NULL;
    // int events[4] = {PAPI_TOT_CYC, PAPI_TOT_INS, PAPI_LST_INS};
    // long long int counters[3];
    // PAPI_library_init(PAPI_VER_CURRENT);
    // PAPI_create_eventset(&event_set);
    // PAPI_add_events(event_set, events, 3);

   
    CahnHilliard *sim = new CahnHilliard(N);

    // const long long int t0 = PAPI_get_real_nsec();
    // PAPI_start(event_set);
    sim->integrate();
    // const long long int t1 = PAPI_get_real_nsec();
    // PAPI_stop(event_set, counters);

    sim->dump("steps/step_35000.bin", 35000);
    delete sim;
    
    // const long long total_cycles = counters[0];       // cpu cycles
    // const long long total_instructions = counters[1]; // any
    // const long long total_load_stores = counters[2];  // number of such instructions

    MPI_Finalize();

    // const size_t flops = 35 * N * N;
    // const size_t mem_ops = 18 * N * N;
    // const double twall = (static_cast<double>(t1) - t0) * 1.0e-9; // seconds
    // const double IPC = static_cast<double>(total_instructions) / total_cycles;
    // const double OI =
    //     static_cast<double>(flops) / (total_load_stores * sizeof(double));
    // const double OI_theory =
    //     static_cast<double>(flops) / (mem_ops * sizeof(double));
    // const double float_perf = flops / twall * 1.0e-9; // Gflop/s

    // std::cout << "Total cycles:                 " << total_cycles << '\n';
    // std::cout << "Total instructions:           " << total_instructions << '\n';
    // std::cout << "Instructions per cycle (IPC): " << IPC << '\n';
    // std::cout << "Total load/store:             " << total_load_stores
    //           << " (expected: " << mem_ops << ")\n";
    // std::cout << "Operational intensity:        " << std::scientific << OI
    //           << " (expected: " << OI_theory << ")\n";
    // std::cout << "Performance [Gflop/s]:        " << float_perf << '\n';
    // std::cout << "Wall-time   [micro-seconds]:  " << twall * 1.0e6 << '\n';
    return 0;
}
