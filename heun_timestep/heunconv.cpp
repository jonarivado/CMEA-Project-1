#include <Eigen/Core>
#include <vector>
#include <iostream>
#include "writer.hpp"
//I added cmath for the exponent
#include <cmath>

/// Uses the Heun's method to compute u from time 0 to time T
/// for the ODE $u'=e^{-2t}-2u$
///
/// @param[out] u at the end of the call (solution at all time steps)
/// @param[out] time contains the time levels
/// @param[in] u0 the initial data
/// @param[in] dt the step size
/// @param[in] T the final time up to which to compute the solution.
///

double TrueSol(double time)
{
    return time * exp(-2 * time);
}

void Heun(std::vector<double> &u, std::vector<double> &time,
          const double &u0, double dt, double T)
{
    const unsigned int nsteps = ceil(T / dt);
    // (write your solution here)
    u.resize(nsteps + 1);
    time.resize(nsteps + 1);
    double y1;
    double y2;
    time.at(0) = 0;
    u.at(0) = u0;
    double Time = 0;
    for (int i = 0; i < nsteps; i++)
    {
        y1 = u.at(i);
        y2 = u.at(i) + dt * (exp(-2 * Time) - 2 * u.at(i));
        u.at(i + 1) = u.at(i) + dt / 2 * (exp(-2 * Time) - 2 * y1) + dt / 2 * (exp(-2 * Time + dt) - 2 * y2);
        Time += dt;
        time.at(i + 1) = Time;
    }
}

int main(int argc, char **argv)
{

    double T = 10.0;

    double dt = 0.2;

    // To make some plotting easier, we take the dt parameter in as an optional
    // parameter.
    if (argc == 2)
    {
        dt = atof(argv[1]);
    }

    const double u0 = 0.;
    std::vector<double> time, u, truesol(8), error(8), timestep(8);

    for (int i = 1; i <= 8; i++)
    {
        Heun(u, time, u0, pow(2, -i), T);
        truesol.at(i - 1) = TrueSol(T);
        error.at(i - 1) = std::abs(u.at(u.size() - 1) - truesol.at(i - 1));
        timestep.at(i - 1) = pow(2, -i);
    }
    writeToFile("solution.txt", u);
    writeToFile("time.txt", timestep);
    writeToFile("truesol.txt", truesol);
    writeToFile("error.txt", error);

    return 0;
}
