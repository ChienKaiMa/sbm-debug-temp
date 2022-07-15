
#include <iostream>
#include <string>

class NumberPartitionProblem
{
protected:
    /* data */
    bool _is_from_file = false;
    std::string _filename;
    uint _problem_size = 0;
    uint _padding = 0;
    float* _numbers;
public:
    NumberPartitionProblem(/* args */); // TODO
    NumberPartitionProblem(std::string filename, float *& Q_flatten, float **& Q);
    ~NumberPartitionProblem();
    bool fromFile(std::string filename, float *& Q_flatten, float **& Q);
    uint getProblemSize() { return _problem_size; }
    uint getMatrixSize() { return _problem_size + _padding; }
    float getNumber(int index) { return _numbers[index]; }
    void printNumbers();
};
