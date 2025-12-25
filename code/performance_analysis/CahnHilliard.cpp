// File       : CahnHilliard.cpp
// Description: CahnHilliard model implementation
// Copyright 2023 Harvard University. All Rights Reserved.
#include "CahnHilliard.h"
#include <cassert>
#include <iostream>
#include <mpi.h>
#include <string>

// MPI_Cart_shift2 Customized function implementation
int MyMPI_Cart_shift2(MPI_Comm cart_comm, int disp0, int disp1,
                      int *rank_source, int *rank_dest) {
    int rank, coords[2];
    MPI_Comm_rank(cart_comm, &rank);
    MPI_Cart_coords(cart_comm, rank, 2, coords);
    
    coords[0] += disp0;
    coords[1] += disp1;
    MPI_Cart_rank(cart_comm, coords, rank_dest);

    coords[0] -= 2*disp0;
    coords[1] -= 2*disp1;
    MPI_Cart_rank(cart_comm, coords, rank_source);

    return 0;
}

// Constructor
CahnHilliard::CahnHilliard(const int N)
    : nodes_(N)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // if (0 == rank) {
    //         std::cout<<"running with "<<size<<" ranks"<<"\n" ;
    //     }
    
    p_ = static_cast<int>(std::sqrt(size));
    if (size != p_ * p_) {
        if (0 == rank) {
            std::cerr << "ERROR: number of ranks must be a square number\n";
        }
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_TOPOLOGY);
    }

    if (0 == rank) {
        std::cout<<"Job size N = "<<N<<"\n";
    }

    const int dims[2] = {p_, p_};
    const int period[2] = {true, true};
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, period, true, &cart_comm_);

    MPI_Cart_shift(cart_comm_, 0, 1, &neighbor_[Top], &neighbor_[Bottom]);
    MPI_Cart_shift(cart_comm_, 1, 1, &neighbor_[Left], &neighbor_[Right]);

    int cart_rank;
    MPI_Comm_rank(cart_comm_, &cart_rank);
    MPI_Cart_coords(cart_comm_, cart_rank, 2, mpi_coords_);

    MyMPI_Cart_shift2(cart_comm_, 1, 1, &neighbor_[LT], &neighbor_[RB]);
    MyMPI_Cart_shift2(cart_comm_, -1, 1, &neighbor_[LB], &neighbor_[RT]);


    // allocate memory for computational grid
    stride_ = nodes_ + 4; // number of nodes plus 4 ghost cells
                          // per dimension (per rank)
    u_ = new double[stride_ * stride_];
    utmp_ = new double[stride_ * stride_]; // for time integration
    derivative_mu_ = new double[stride_ * stride_]; 
    derivative_mu_tmp_ = new double[stride_ * stride_]; // for time integration

    top_  = new double[nodes_*2];
    bottom_ = new double[nodes_*2];
    left_   = new double[nodes_*2];
    right_  = new double[nodes_*2];

    toprcv_  = new double[nodes_*2];
    bottomrcv_ = new double[nodes_*2];
    leftrcv_   = new double[nodes_*2];
    rightrcv_  = new double[nodes_*2];

    // MPI type used in file view and local data
    MPI_Type_vector(nodes_, nodes_, p_ * nodes_, MPI_DOUBLE, &MyFileType_);
    MPI_Type_commit(&MyFileType_);

    MPI_Type_vector(nodes_, nodes_, stride_, MPI_DOUBLE, &MyGridType_);
    MPI_Type_commit(&MyGridType_);

    // initialize grid
    initialize_();
}

// Destructor
CahnHilliard::~CahnHilliard()
{
    delete[] u_;                // free memory
    delete[] utmp_;             // free memory
    delete[] derivative_mu_ ;    // free memory
    delete[] derivative_mu_tmp_ ; // free memory
    delete[] top_;             // free memory
    delete[] bottom_;             // free memory
    delete[] left_;             // free memory
    delete[] right_;             // free memory
    delete[] toprcv_;             // free memory
    delete[] bottomrcv_;             // free memory
    delete[] leftrcv_;             // free memory
    delete[] rightrcv_;             // free memory
    MPI_Comm_free(&cart_comm_); // free Cartesian communicator
    MPI_Type_free(&MyFileType_);
    MPI_Type_free(&MyGridType_);
}

// Public methods
void CahnHilliard::dump(const std::string fname, int it_num)
{
    if (it_num % 100 == 0) {
        MPI_File fh;
        MPI_Status st;

        if (MPI_File_open(cart_comm_,
                        fname.c_str(),
                        MPI_MODE_CREATE | MPI_MODE_WRONLY,
                        MPI_INFO_NULL,
                        &fh)) {
            std::cerr << "Can not create MPI file '" << fname << "'\n";
            MPI_Abort(cart_comm_, MPI_ERR_BAD_FILE);
        }
        
        int cart_rank;
        MPI_Comm_rank(cart_comm_, &cart_rank);
        const int dim = p_ * nodes_;
        if (0 == cart_rank) {
            // write file dimension
            MPI_File_write_at(fh, 0, &dim, 1, MPI_INT, &st);
        }
        MPI_Offset offset = sizeof(dim) + (mpi_coords_[0] * dim + mpi_coords_[1]) *
                                            nodes_ * sizeof(double);
        MPI_File_set_view(
            fh, offset, MPI_DOUBLE, MyFileType_, "native", MPI_INFO_NULL);
        MPI_File_write_all(fh, &u(0, 0), 1, MyGridType_, &st);
        MPI_File_close(&fh);
    }
}
