
#include "simulatedBifurcationOptimizer.hpp"
#include "numberpartition.hpp"
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <random>
#define dcal_t float

/*
template <class T>
void checkSBMCoeff(T** J, int matrix_size, float c0)
{
    std::cout << "Coefficient check\n";
    float c0_debug = 0;
    for (int i = 0; i < matrix_size; i++)
    {
        for (int j = 0; j < matrix_size; j++)
        {
            c0_debug += J[i][j] * J[i][j];
        }
    }
    c0_debug = 0.5 * sqrt(matrix_size / c0_debug);
    std::cout << "c0_theoretical: " << c0_debug << "\n";
    std::cout << "c0_in_use     : " << c0 << "\n";
    std::cout << "End of coefficient check\n\n";
}
*/

int main(int argc, char *argv[])
{
    auto timer = std::time(nullptr);
    auto start_time = std::localtime(&timer);
    std::cerr << std::put_time(start_time, "%Y-%m-%d %T");
    std::cerr << " Test starts.\n";

    // --Read problem information--
    std::string problemFilePath = "num_par_256_ising_0.txt";
    if (argc == 2) problemFilePath = argv[1];

    float *Q_flatten;
    float **Q;
    NumberPartitionProblem problem = NumberPartitionProblem(problemFilePath, Q_flatten, Q);

    auto timer1 = std::time(nullptr);
    auto read_time = std::localtime(&timer1);
    std::cerr << std::put_time(read_time, "%Y-%m-%d %T");
    std::cerr << " Read problem information done.\n";

    // TODO: Measure problem read time

    /*
    ** Simulated bifurcation algorithm setup
    */
    int steps = 200;
    // float dt = (float)200 / (float)steps;
    dcal_t dt = 0.01;
    // dcal_t a0 = 0.5;
    dcal_t c0 = 0.001;
    // dcal_t Kerr = 1.;
    uint problem_size = problem.getProblemSize();
    uint matrix_size = problem.getMatrixSize();
    dcal_t x[matrix_size] = {0};
    dcal_t y[matrix_size] = {0};
    for (int i = 0; i < matrix_size; ++i) {
        y[i] = 0.05;
    }
    dcal_t best_energy = MAXFLOAT;
    int best_step = 0;
    uint best_spin[matrix_size / 32] = {0};
    float myC0 = 4.2479e-06;
    float myDT = 0.01;
    float delta_a[steps]{0};
    for (int i = 0; i < steps; ++i)
    {
        delta_a[i] = float(steps - i) / steps;
    }

    float* x_history = new float[matrix_size * steps];

    auto t1 = clock();
    sbOptimizerTop(Q_flatten, x, y, delta_a, myC0, myDT, matrix_size, steps, best_spin, x_history);
    auto t2 = clock();
    
    auto timer2 = std::time(nullptr);
    auto execute_time = std::localtime(&timer2);
    std::cerr << std::put_time(execute_time, "%Y-%m-%d %T") << " Simulated bifurcation done.\n";
    std::cout << "Execution time = " << ((float)(t2 - t1))/CLOCKS_PER_SEC << " seconds.\n\n";

    // checkSBMCoeff<float>(Q, matrix_size, myC0);

    //
    // Direct SBM history (x) to file
    //
    std::cout << "---Dumping x history---\n";
    std::ofstream x_outfile;
    std::ofstream energy_outfile;
    // outfile.open(argv[2], std::ios::out);
    // TODO
    // Wrap the utilities in functions
    x_outfile.open("x_history.log", std::ios::out);

    for (int i = 0; i < steps; ++i)
    {
        for (int j = 0; j < matrix_size; ++j)
        {
            x_outfile << x_history[i * matrix_size + j] << " ";
        }
        x_outfile << "\n";
    }

    std::cout << "---Dumping x history done---\n\n";

    //
    // SBM summary
    //

    // TODO
    // Save the report as file
    std::cout << "---SBM summary---\n";
    // TODO: Print overflow counts
    std::cout << "Best solution at step " << best_step << " of " << steps << "\n";
    std::cout << "Best energy: " << best_energy << "\n";
    std::cout << "Best spin: ";
    for (int i = 0; i < (matrix_size / 32); ++i)
    {
        ap_uint<32> u = best_spin[i];
        for (int j = 0; j < 32; ++j) {
            std::cout << u[j] << " ";
        }
    }
    std::cout << "\n";
    std::cout << "\n";

    std::cout << "---Number partition result---\n";
    // Members of the first group
    std::cout << "0: ";
    for (int i = 0; i < (problem_size / 32); ++i)
    {
        ap_uint<32> u = best_spin[i];
        for (int j = 0; j < 32; ++j) {
            if (!u[j]) {
                std::cout << problem.getNumber(j) << " ";
            }
        }
    }
    std::cout << "\n";
    // Members of the second group
    std::cout << "1: ";
    for (int i = 0; i < (problem_size / 32); ++i)
    {
        ap_uint<32> u = best_spin[i];
        for (int j = 0; j < 32; ++j) {
            if (u[j]) {
                std::cout << problem.getNumber(j) << " ";
            }
        }
    }
    std::cout << "\n";
    float diff = 0;
    for (int i = 0; i < (problem_size / 32); ++i)
    {
        ap_uint<32> u = best_spin[i];
        for (int j = 0; j < 32; ++j) {
            diff += problem.getNumber(32*i + j) * (2 * u[j] - 1);
        }
    }
    for (int i = 0; i < problem_size; ++i)
    {
        ap_uint<32> u = best_spin[0];
    }
    std::cout << "Difference: " << diff << "\n";

    // Clean up problem matrix if it is dynamically allocated
    // delete[] numbers;
    for (int i = 0; i < problem_size; i++)
    {
        delete[] Q[i];
    }
    delete[] Q;
    delete[] Q_flatten;
}
