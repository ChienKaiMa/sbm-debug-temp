

#include "simulatedBifurcationOptimizer.hpp"

typedef union {
    unsigned int u;
    float f;
} ap32;

void SimulatedBifurcationOptimizer::bound(float& x, float& y) {
    if (x > 1) {
        x = 1;
        y = 0;
    } else if (x < -1) {
        x = -1;
        y = 0;
    }
}

void SimulatedBifurcationOptimizer::discretize(float* x, bool* signOfX) {
    for (unsigned int i = 0; i < _qubits; ++i) {
DO_PRAGMA(HLS loop_tripcount min=64 max=MAX_QUBITS)
        signOfX[i] = (x[i] > 0);
    }
}

void SimulatedBifurcationOptimizer::matrixVectorProduct(int blockNum, unsigned short blockSize) {
    // TODO: blockSize tells the program how many _signOfX should we get
Multiply_matrix:
    for (unsigned int i = 0; i < CACHE_FACTOR; ++i) {
    Multiply_rows:
        for (unsigned int j = 0; j < CACHE_FACTOR; ++j) {
            ap32 temp = {0};
            temp.f = _jCache[i][j];
            ap_uint<32> tempu = temp.u;
            // Flip the sign bit when sign of x is -1 (0)
            tempu[31] = tempu[31] ^ !(_signOfX[blockNum * CACHE_FACTOR + j]);
            temp.u = tempu;
            if (i == 0)
                _productCache[j] = temp.f;
            else
                _productCache[j] += temp.f;
        }
    }
}

void SimulatedBifurcationOptimizer::simulatedBifurcationStep(float* isingJ, float* x, float* y) {
Discretize_all:
    discretize(x, _signOfX);
    int counter = 0;
Update_blocks:
    for (unsigned int i = 0; i < _qubits / float(CACHE_FACTOR); ++i) {
DO_PRAGMA(HLS loop_tripcount min=1 max=4)
// TODO 補 a0
// constexpr, MACRO...?
        // TODO: Prefetch x and y to xCache and yCache
    Prepare_data_xy:
        for (unsigned int k = 0; k < CACHE_FACTOR; ++k) {
            _yCache[k] = y[counter + k];
            _xCache[k] = x[counter + k];
        }
    Update_y_stage1:
        for (unsigned int k = 0; k < CACHE_FACTOR; ++k) {
            _yCache[k] -= _xCache[k] * _c2;
        }
    Update_y_stage2:
        for (unsigned int j = 0; j < _qubits / CACHE_FACTOR; ++j) {
    #pragma HLS loop_tripcount min=1 max=16
            // TODO: Prefetch J to jCache
        Prepare_data_j:
            for (unsigned int k = 0; k < CACHE_FACTOR; ++k) {
                for (unsigned int l = 0; l < CACHE_FACTOR; ++l) {
                    _jCache[k][l] = isingJ[(i * CACHE_FACTOR + k) * CACHE_FACTOR + j * CACHE_FACTOR + l];
                }
            }
            matrixVectorProduct(j, CACHE_FACTOR); // TODO: The last block may have blockSize less than CACHE_FACTOR
        }
    Update_and_write_xy:
        for (unsigned int k = 0; k < CACHE_FACTOR; ++k) {
            _yCache[k] -= _productCache[k] * _c1;
            _xCache[k] += _yCache[k] * _dt;
            bound(_xCache[k], _yCache[k]);
            // matrixVectorProduct initializes the variable
            // _productCache[k] = 0;
            x[counter + k] = _xCache[k];
            y[counter + k] = _yCache[k];
        }
        counter += CACHE_FACTOR;
    }

}

void SimulatedBifurcationOptimizer::runSimulation(float* isingJ, float* x, float* y, float* delta_a, unsigned int* solution, float* x_history) {
Simulated_bifurcation_steps:
    for (unsigned int i = 0; i < _steps; ++i) {
#pragma HLS loop_tripcount min=100 max=2000
        _c2 = delta_a[i] * _dt;
        simulatedBifurcationStep(isingJ, x, y);
        for (unsigned int j = 0; j < _qubits; ++j)
        {
            x_history[i * _qubits + j] = x[j];
        }
    }
    discretize(x, _best_solution);
    packSolution(solution, _best_solution);
}

void SimulatedBifurcationOptimizer::runSimulation(float* isingJ, unsigned short numOfSteps) {
    // If numOfQubits is large, where should we store x and y?
    // float* x;
    // float* y;
    // runSimulation(isingJ, x, y, numOfSteps);
}

void SimulatedBifurcationOptimizer::setParams(int paramCode, float& param) {
    switch (paramCode) {
    case 0:
        _dt = param;
        break;
    case 1:
        _c0 = param;
        break;
    }
    _c1 = 2 * _c0 * _dt;
}

void SimulatedBifurcationOptimizer::setNumOfQubits(unsigned int numOfQubits) {
    // TODO Error-handling
    // If numOfQubits > MAX_QUBITS, set errorCode to non-zero and prohibit any calculation
    _qubits = numOfQubits;
}

void SimulatedBifurcationOptimizer::setNumOfSteps(unsigned short numOfSteps) {
    _steps = numOfSteps;
}

void SimulatedBifurcationOptimizer::getParams(int paramCode, float& param) {
    switch (paramCode) {
    case 0:
        param = _dt;
        break;
    case 1:
        param = _c0;
        break;
    }
}

void SimulatedBifurcationOptimizer::packSolution(unsigned int * solution_uint, bool * solution_bool) {
Copy_solution_blocks:
    for (int i = 0; i < _qubits / 32; ++i) {
#pragma HLS loop_tripcount min=2 max=32
        ap_uint<32> u = {0};
    Copy_spins:
        for (int j = 0; j < 32; ++j) {
            u[j] = solution_bool[32 * i + j];
        }
#ifndef __SYNTHESIS__
        std::cout << "u = " << u << "\n";
#endif
        solution_uint[i] = u;
    }
}


SimulatedBifurcationOptimizer::SimulatedBifurcationOptimizer()
{

}


SimulatedBifurcationOptimizer::~SimulatedBifurcationOptimizer()
{
}

extern "C" void sbOptimizerTop(float* isingJ,
                               float* x,
                               float* y,
                               float* delta_a,
                               float& c0,
                               float& dt,
                               unsigned int numOfQubits,
                               unsigned short numOfSteps,
                               unsigned int* solution,
                               float* x_history)
{
    static SimulatedBifurcationOptimizer kernel;

#pragma HLS DATAFLOW disable_start_propagation

    kernel.setParams(0, dt);
    kernel.setParams(1, c0);
    kernel.setNumOfQubits(numOfQubits);
    kernel.setNumOfSteps(numOfSteps);

    kernel.runSimulation(isingJ,
                         x,
                         y,
                         delta_a,
                         solution,
                         x_history);

}
