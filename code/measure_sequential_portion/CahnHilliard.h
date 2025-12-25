// File       : CahnHilliard.h
// Description: CahnHilliard model API
// Copyright 2023 Harvard University. All Rights Reserved.
#include <cassert>
#include <cmath>
#include <iostream>
#include <mpi.h>
#include <string>
#include <vector>
#include <omp.h>
#include <numeric>
#include <random>
#include <sys/time.h>




const double INITIAL_MEAN = 0.5;
const double INITIAL_VARIANCE = 0.01;
const double DIFFUSION_COEFFICIENT = 1.0;
const double MU_SCALE = 1.0;
const double C_MINUS = 0.0;
const double C_PLUS = 1.0;
const double PENALTY_SCALE = 0.2;
const double SPATIAL_SCALE = 1.0;

class CahnHilliard
{
public:
    CahnHilliard(const int N = 350);

    virtual ~CahnHilliard();

    // Public methods
    double &u(const int i, const int j)
    {
        assert(i >= -2 && i < nodes_ + 2);
        assert(j >= -2 && j < nodes_ + 2);
        return u_[(i + 2) * stride_ + j + 2];
    }

    double &derivative_mu(const int i, const int j)
    {
        assert(i >= -2 && i < nodes_ + 2);
        assert(j >= -2 && j < nodes_ + 2);
        return derivative_mu_[(i + 2) * stride_ + j + 2];
    }

    void dump(const std::string fname);

    void integrate()
    {
        int cart_rank;
        MPI_Comm_rank(cart_comm_, &cart_rank);
        const bool isroot = (0 == cart_rank);

        //double dexp = -5;
        double elapsed = 0;
        double duration = 35;
        int it_num = 0;
        double dt = 0.001;
        
        while (elapsed < duration) {
            if (isroot && it_num % 100 == 0) {
                std::cout << "step = " << it_num << '\t' << "time = " << elapsed
                          << '\n';
            }

            //dt = std::min(100.0, std::exp(dexp));
            elapsed += dt;
            //dexp += 0.01;
            it_num++;
            // update solution in time
            update_async_(dt);

           
            // swap data arrays and advance time
            std::swap(u_, utmp_);
            std::swap(derivative_mu_, derivative_mu_tmp_);
        }
        std::cout<<"parallel time is"<<t1_ - t0_<<"\n";
    }

private:
    // Data
    const int nodes_;
    int p_, stride_;
    double *u_, *utmp_, *derivative_mu_, *derivative_mu_tmp_; // species u (including ghosts, per rank)

    double *top_, *bottom_, *left_,*right_;
    double *toprcv_, *bottomrcv_, *leftrcv_, *rightrcv_;

    double t0_ = 0.0;
    double t1_= 0.0;
   

    // MPI
    MPI_Comm cart_comm_; // Cartesian communicator
    int mpi_coords_[2];  // rank coordinates in Cartesian topology
    MPI_Datatype MyFileType_; // MPI type used in file view
    MPI_Datatype MyGridType_; // MPI type used for local data
    // enum Neighbor { Left, Right, Bottom, Top, NNeighbors };
    enum Neighbor { Left, Right, Bottom, Top, LB, LT, RB, RT, NNeighbors };
    int neighbor_[NNeighbors];

    // Private methods
    double get_wtime(void)
    {
    struct timeval t;
    gettimeofday(&t, NULL);
    return (double)t.tv_sec + (double)t.tv_usec * 1.0e-6; // seconds
    }


    void initialize_()
    {
    t0_ += get_wtime();
    #pragma omp parallel  
    { 
        std::default_random_engine generator(omp_get_thread_num()+(mpi_coords_[0]*32+mpi_coords_[1])*32);
        std::normal_distribution<double> distribution(INITIAL_MEAN, std::sqrt(INITIAL_VARIANCE)); 
    #pragma omp for nowait schedule(static)
        for (int i = 0; i < nodes_; ++i) {
            for (int j = 0; j < nodes_; ++j) {
                u(i, j) = distribution(generator);
                derivative_mu(i,j) = derive_mu(u(i,j));
            }
        }
    }
    t1_ += get_wtime();
    }

    double derive_mu(const double c)
    {
        return MU_SCALE * (c - C_MINUS) * (c - C_PLUS) * (2 * c - C_MINUS - C_PLUS); 
    }

    // Spatio-temporal finite difference schemes for species u (will be
    // inlined)
    void update_u_(const int i, const int j, const double dt)
    {

        utmp_[(i + 2) * stride_ + j + 2] =
            u(i, j) +
            // suboptimal for now, mu computed 4 times
            dt * DIFFUSION_COEFFICIENT * 
            (derivative_mu(i+1,j) + derivative_mu(i-1,j) + derivative_mu(i,j+1) + derivative_mu(i,j-1) - 4 * derivative_mu(i,j) 
            - PENALTY_SCALE*(
                       u(i+2, j) + u(i-2, j) + u(i, j+2) + u(i, j-2) 
                + 2 * (u(i+1, j+1) + u(i+1, j-1) + u(i-1, j+1) + u(i-1, j-1))
                - 8 * (u(i+1, j) + u(i-1, j) + u(i, j+1) + u(i, j-1))
                + 20 * u(i,j))
                /(SPATIAL_SCALE*SPATIAL_SCALE))
            /(SPATIAL_SCALE*SPATIAL_SCALE) ;
        
        derivative_mu_tmp_[(i + 2) * stride_ + j + 2] = derive_mu( utmp_[(i + 2) * stride_ + j + 2] );
    }



    void update_async_(const double dt)
    {
        // Algorithm:
        // 1.) Pack ghosts (see pack_() below)
        pack_();
        // 2.) Initiate asynchronous communication
        MPI_Request reqs[16];
        MPI_Irecv(&leftrcv_[0], nodes_*2, MPI_DOUBLE, neighbor_[Left], 0, cart_comm_, &reqs[0]);
        MPI_Irecv(&rightrcv_[0], nodes_*2, MPI_DOUBLE, neighbor_[Right], 1, cart_comm_, &reqs[1]);
        MPI_Irecv(&bottomrcv_[0], nodes_*2, MPI_DOUBLE, neighbor_[Bottom], 2, cart_comm_, &reqs[2]);
        MPI_Irecv(&toprcv_[0], nodes_*2, MPI_DOUBLE, neighbor_[Top], 3, cart_comm_, &reqs[3]);
        
        MPI_Isend(&left_[0], nodes_*2, MPI_DOUBLE, neighbor_[Left], 1, cart_comm_, &reqs[4]);
        MPI_Isend(&right_[0], nodes_*2, MPI_DOUBLE, neighbor_[Right], 0, cart_comm_, &reqs[5]);
        MPI_Isend(&bottom_[0], nodes_*2, MPI_DOUBLE, neighbor_[Bottom], 3, cart_comm_, &reqs[6]);
        MPI_Isend(&top_[0], nodes_*2, MPI_DOUBLE, neighbor_[Top], 2, cart_comm_, &reqs[7]);

////////////////////////////
        // Send and receive diagonal elements
        MPI_Irecv(&u(-1,-1), 1, MPI_DOUBLE, neighbor_[LT], 4, cart_comm_, &reqs[8]);
        MPI_Irecv(&u(nodes_,-1), 1, MPI_DOUBLE, neighbor_[LB], 5, cart_comm_, &reqs[9]);
        MPI_Irecv(&u(nodes_,nodes_), 1, MPI_DOUBLE, neighbor_[RB], 6, cart_comm_, &reqs[10]);
        MPI_Irecv(&u(-1,nodes_), 1, MPI_DOUBLE, neighbor_[RT], 7, cart_comm_, &reqs[11]);
        
        MPI_Isend(&u(0,0), 1, MPI_DOUBLE, neighbor_[LT], 6, cart_comm_, &reqs[12]);
        MPI_Isend(&u(nodes_-1,0), 1, MPI_DOUBLE, neighbor_[LB], 7, cart_comm_, &reqs[13]);
        MPI_Isend(&u(nodes_-1,nodes_-1), 1, MPI_DOUBLE, neighbor_[RB], 4, cart_comm_, &reqs[14]);
        MPI_Isend(&u(0,nodes_-1), 1, MPI_DOUBLE, neighbor_[RT], 5, cart_comm_, &reqs[15]);
////////////////////////////

        // 3.) Update interior domain
    t0_ += get_wtime();
    #pragma omp parallel for schedule(static)    
        for (int i = 2; i < nodes_ - 2; ++i) {
            for (int j = 2; j < nodes_ - 2; ++j) {
                update_u_(i, j, dt);
            }
        }
    t1_ += get_wtime();

        // 4.) Wait for completion of MPI operations
        MPI_Waitall(16, reqs, MPI_STATUSES_IGNORE);

        // 5.) Unpack ghosts (see unpack_() below)
        unpack_();

        // 6.) Update boundary domain
    t0_ += get_wtime();
    #pragma omp parallel for schedule(static)        
        for (int i = 0; i < nodes_; ++i) {
            for (int j: {0, 1, nodes_-2, nodes_-1}) {
                update_u_(i, j, dt);
                update_u_(j, i, dt);
            }
        }
    t1_ += get_wtime();

    }

    void pack_()
    {
        // Pack domain data (ghosts) into send buffers
    t0_ += get_wtime();
    #pragma omp parallel for schedule(static)        
        for (int i = 0; i < nodes_; ++i) {
            bottom_[i] = u(nodes_ - 1, i);
            top_[i] = u(0, i);
            right_[i] = u(i, nodes_ - 1);
            left_[i] = u(i, 0);

            bottom_[nodes_ + i] = u(nodes_ - 2, i);
            top_[nodes_ + i] = u(1, i);
            right_[nodes_ + i] = u(i, nodes_ - 2);
            left_[nodes_ + i] = u(i, 1);
        }
    t1_ += get_wtime();
        


    }

    void unpack_()
    {
        // Unpack receive buffers (ghosts) into domain data
    t0_ += get_wtime();
    #pragma omp parallel for schedule(static)        
        for (int i = 0; i < nodes_; ++i) {
            u(nodes_, i) = bottomrcv_[i];
            derivative_mu(nodes_, i) = derive_mu(u(nodes_, i));
            u(-1, i) = toprcv_[i];
            derivative_mu(-1, i) = derive_mu(u(-1, i));
            u(i, nodes_) = rightrcv_[i];
            derivative_mu(i, nodes_) = derive_mu(u(i, nodes_));
            u(i, -1) = leftrcv_[i];
            derivative_mu(i, -1) = derive_mu(u(i, -1));

            u(nodes_ + 1, i) = bottomrcv_[nodes_ + i];
            derivative_mu(nodes_ + 1, i) = derive_mu(u(nodes_ + 1, i));
            u(-2, i) = toprcv_[nodes_ + i];
            derivative_mu(-2, i) = derive_mu(u(-2, i));
            u(i, nodes_ + 1) = rightrcv_[nodes_ + i];
            derivative_mu(i, nodes_ + 1) = derive_mu(u(i, nodes_ + 1));
            u(i, -2) = leftrcv_[nodes_ + i];
            derivative_mu(i, -2) = derive_mu(u(i, -2));


        }
    t1_ += get_wtime();
    }
};
