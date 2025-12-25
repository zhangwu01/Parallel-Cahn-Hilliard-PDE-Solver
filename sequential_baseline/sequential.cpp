#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <vector>
#include <algorithm>

// define constants
const double INITIAL_MEAN = 0.5;
const double INITIAL_VARIANCE = 0.01;
const int NX = 100;
const int NY = 100;
const double DX = 1.0;
const double DY = 1.0;
const double DIFFUSION_COEFFICIENT = 1.0;
const double MU_SCALE = 1.0;
const double C_MINUS = 0.0;
const double C_PLUS = 1.0;
const double PENALTY_SCALE = 0.2;
const int nn = 1;

using namespace std;

vector<vector<double>> expand_c(vector<vector<double>> c) 
{
    int NX = c.size();
    int NY = c[0].size();

    vector<vector<double>> c_expand(nn + NX + nn, vector<double>(nn + NY + nn, 0)); // (nn+NX+nn)*(nn+NY+nn)*sizeof(double) bytes

    for (int i = 1; i <= NX; i++)
    {
        for (int j = 1; j <= NY; j++)
        {
            c_expand[i + nn - 1][j + nn - 1] = c[i - 1][j - 1]; //1 memory access, NX*NY times
        }
    }

    for (int j = nn; j <= nn + NY - 1; j++)
    {
        c_expand[0][j] = c[NX - 1][j - nn]; //1 memory access, NY times
        c_expand[nn + NX][j] = c[0][j - nn]; //1 memory access, NY times
    }

    for (int i = nn; i <= nn + NX - 1; i++)
    {
        c_expand[i][0] = c[i - nn][NY - 1]; //1 memory access, NX times
        c_expand[i][nn + NY] = c[i - nn][0]; //1 memory access, NX times
    }

    return c_expand;
}

vector<vector<double>> laplacian(vector<vector<double>> c) //
{
    int NX = c.size();
    int NY = c[0].size();
    vector<vector<double>> c_expand = expand_c(c); // (nn+NX+nn)*(nn+NY+nn)*sizeof(double) bytes
    vector<vector<double>> claplace(NX, vector<double>(NY, 0));

    for (int i = 0; i < NX; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            double cxx = (c_expand[i + nn][j + nn + 1] + c_expand[i + nn][j + nn - 1] - 2 * c_expand[i + nn][j + nn]) / (DX * DX); //5 FLOPS, 4 memory accesses
            double cyy = (c_expand[i + nn + 1][j + nn] + c_expand[i + nn - 1][j + nn] - 2 * c_expand[i + nn][j + nn]) / (DY * DY); //5 FLOPS, 4 memory accesses
            claplace[i][j] = cxx + cyy; //1 FLOP, 1 memory access
        }
    }

    return claplace;
}

std::vector<std::vector<double>> derivative_mu(std::vector<std::vector<double>> c, double cm = 0, double cp = 1, double a_squared = 1) 
{
    /*
    The derivative of chemical potential with respect to concentration
    input: c -- 2d vector, concentration
           cm -- concentration corresponding to the left well
           cp -- concentration corresponding to the right well
           a_squared -- energy scale
    */
    int nrow = c.size();
    int ncol = c[0].size();
    std::vector<std::vector<double>> partial_mu(nrow, std::vector<double>(ncol, 0));
    for (int i = 0; i < nrow; i++)
    {
        for (int j = 0; j < ncol; j++)
        {
            partial_mu[i][j] = a_squared * (c[i][j] - cm) * (c[i][j] - cp) * (2 * c[i][j] - cm - cp); 
        }
    }
    return partial_mu;
}

void ch_kernal(std::vector<std::vector<double>> &c, double dt, double D, double a_squared, double cm, double cp, double epsilon_squared, int NX, int NY) // 46*NX*NY FLOPS, 33*NX*NY memory accesses
{
    // Allocate memory for temporary concentration field
    std::vector<std::vector<double>> c_new(NX, std::vector<double>(NY)); // NX*NY*sizeof(double) bytes

    // Compute derivative of chemical potential
    std::vector<std::vector<double>> der_mu = derivative_mu(c, cm, cp, a_squared); //8*NX*NY FLOPS, 4*NX*NY memory accesses

    // Compute Laplacian of concentration
    std::vector<std::vector<double>> lap_c = laplacian(c); //11*NX*NY FLOPS, 8 memory accesses

    // Compute Laplacian of Laplacian of concentration
    std::vector<std::vector<double>> lap2_c = laplacian(lap_c); //11*NX*NY FLOPS, 8 memory accesses

    // Compute laplacian of der_mu
    std::vector<std::vector<double>> lap_der_mu = laplacian(der_mu); //11*NX*NY FLOPS, 8 memory accesses, NX*NY times.

    // Update concentration field
    for (int i = 0; i < NX; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            c_new[i][j] = c[i][j] + D * dt * (lap_der_mu[i][j] - epsilon_squared * lap2_c[i][j]); // 5 FLOPS, 4 memory accesses
        }
    }

    // Copy updated concentration field to original concentration field
    c = c_new; // 1 memory access
}

// main function for computation

int main()
{

    // set up initial conditions
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(INITIAL_MEAN, std::sqrt(INITIAL_VARIANCE));
    std::vector<std::vector<double>> Initial_Conditions(NX, std::vector<double>(NY));
    for (int i = 0; i < NX; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            Initial_Conditions[i][j] = distribution(generator);
        }
    }

    // set up variables
    double dexp = -5;
    double elapsed = 0;
    double duration = 35;
    int it_num = 0;

    // iterate
    while (elapsed < duration)
    {
        // save to file
        if (it_num % 10 == 0)
        {
            std::ofstream outfile;
            outfile.open("test_data/raw/timepoint_" + std::to_string((int)elapsed) + ".csv");

            // add column names to the csv file
            for (int j = 0; j < NY - 1; j++)
            {
                outfile << "Column" << j + 1 << ",";
            }
            outfile << "Column" << NY;
            outfile << "\n";

            for (int i = 0; i < NX; i++)
            {
                for (int j = 0; j < NY - 1; j++)
                {
                    outfile << Initial_Conditions[i][j] << ",";
                }
                outfile << Initial_Conditions[i][NY - 1];
                outfile << "\n";
            }
            outfile.close();
        }

        // update variables
        double dt = std::min(100.0, std::exp(dexp));
        elapsed += 0.1;
        dexp += 0.01;
        it_num++;

        ch_kernal(Initial_Conditions, dt, DIFFUSION_COEFFICIENT, MU_SCALE, C_MINUS, C_PLUS, PENALTY_SCALE, NX, NY); 

    }

    std::cout << it_num << " iterations were computed\n";
    std::cout << "Elapsed time: " << elapsed << " seconds\n";

    return 0;
}
