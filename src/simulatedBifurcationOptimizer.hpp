
#ifndef __SB_OPTIMIZER_HPP__
#define __SB_OPTIMIZER_HPP__

// TODO: decide include
#include "ap_fixed.h"

#define PRAGMA_SUB(x) _Pragma (#x)
#define DO_PRAGMA(x) PRAGMA_SUB(x)
#define MAX_QUBITS 4096
#define CACHE_FACTOR 64


class SimulatedBifurcationOptimizer
{
private:
    //
    // Parameters
    //
    // TODO: a0 placeholder
    float           _dt = 0.01;
    float           _c0 = 4.2479e-06;
    float           _c1 = 2 * _c0 * _dt;
    float           _c2 = 0;
    unsigned int    _qubits = 1024; // TODO: To be decided
    unsigned short  _steps = 100; // TODO: To be decided

    //
    // Result
    //
    unsigned short  _best_step = 0;
    unsigned short  _errorCode = 0; // TODO: To be designed
    
    //
    // Data buffer
    //
    bool            _signOfX[MAX_QUBITS] = {};
    // TODO: decide float cache size
    float           _jCache[CACHE_FACTOR][CACHE_FACTOR] = {};
    float           _xCache[CACHE_FACTOR] = {};
    float           _yCache[CACHE_FACTOR] = {};
    float           _productCache[CACHE_FACTOR] = {};
    bool            _best_solution[MAX_QUBITS] = {};
public:
    SimulatedBifurcationOptimizer();
    ~SimulatedBifurcationOptimizer();

    void discretize(float* x, bool* signOfX);
    void matrixVectorProduct(int blockNum, unsigned short blockSize);
    void bound(float& x, float& y);
    void simulatedBifurcationStep(float* isingJ, float* x, float* y);
    // TODO: Allocate memory through DDR for runSimulation function?
    // void runSimulation(float* isingJ, float* x, float* y, unsigned short numOfSteps, float* preserved1, bool* preserved2);
    void runSimulation(float* isingJ, float* x, float* y, float* delta_a, unsigned int* solution, float* x_history);
    void runSimulation(float* isingJ, unsigned short numOfSteps); // TODO: Randomly generate x and y
    void setNumOfQubits(unsigned int numOfQubits);
    void setNumOfSteps(unsigned short numOfSteps);
    void setParams(int paramCode, float& param);
    void getParams(int paramCode, float& param);
    void packSolution(unsigned int* solution_uint, bool* solution_bool); // TODO: Pack bool into unsigned int and write
    // TODO: Modify parameters
};

extern "C" void sbOptimizerTop(float* isingJ,
                               float* x,
                               float* y,
                               float* delta_a,
                               float& c0,
                               float& dt,
                               unsigned int numOfQubits,
                               unsigned short numOfSteps,
                               unsigned int* solution,
                               float* x_history);

#endif