#include "writer.hpp"
#include <Eigen/Core>
#include <Eigen/LU>
#include <vector>
#include <iostream>
#include <chrono>
#include "coronaoutbreak.hpp"
#include "dirksolver.hpp"

//----------------mainBegin----------------

int main(int argc, char **argv)
{

    double T = 101;
    CoronaOutbreak outbreak(0, 0, 0, 0, 0, 0.03, 0.02, 0, 0);
    std::vector<double> u0(5);
    u0[0] = 500;
    u0[1] = 0;
    u0[2] = 0;
    u0[3] = 0;
    u0[4] = 0;

    // Compute the exact solution for the parameters above
    std::vector<double> exact = outbreak.computeExactNoCorona(T, u0[0]);

    // Initialize solver object for the parameters above
    DIRKSolver dirkSolver(outbreak);

    int minExp = 0;
    int maxExp = 8;
    int countExponents = maxExp - minExp + 1;
    std::vector<double> numbers(countExponents);
    std::vector<double> walltimes(countExponents);
    std::vector<double> errors(countExponents);

    // (write your solution here)

    for (int i = minExp; i <= maxExp; i++)
    {

        int N = 200 * std::pow(2, i);
        std::vector<double> time(N + 1);
        std::vector<std::vector<double>> u(5, std::vector<double>(N + 1));
        u.at(0).at(0) = u0.at(0);
        //measure the time
        auto start = std::chrono::high_resolution_clock::now();
        dirkSolver.solve(u, time, T, N);
        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed;
        elapsed = finish - start;
        walltimes.at(i) = elapsed.count();
        numbers.at(i) = N;
        errors.at(i) = 0;
        for (int j = 0; j < 5; j++)
        {
            errors.at(i) += std::abs(u.at(j).at(N) - exact.at(j));
        }
    }

    writeToFile("numbers.txt", numbers);
    writeToFile("errors.txt", errors);
    writeToFile("walltimes.txt", walltimes);
}
//----------------mainEnd----------------
