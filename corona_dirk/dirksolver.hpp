#pragma once
#include <Eigen/Core>
#include <vector>
#include "coronaoutbreak.hpp"

class DIRKSolver
{
public:
	DIRKSolver(const CoronaOutbreak &coronaOutbreak_)
		: coronaOutbreak(coronaOutbreak_)
	{
	}

	///
	/// Evaluates function G1(y) from task b)
	/// @param[out] G evaluation of function
	/// @param[in] y input to G
	/// @param[in] tn current time
	/// @param[in] Un computed value of U = [S,E,I,R] at time tn
	/// @param[in] dt timestep
	///
	//----------------G1Start----------------
	void computeG1(Eigen::VectorXd &G, Eigen::VectorXd y, double tn,
				   Eigen::VectorXd Un, double dt)
	{
		// (write your solution here)
		Eigen::VectorXd F;
		F.resize(4);
		coronaOutbreak.computeF(F, tn + mu * dt, y);
		G = Un + dt * mu * F - y;
	}
	//----------------G1End----------------

	///
	/// Evaluates function G2(y1, y) from task b)
	/// @param[out] G evaluation of function
	/// @param[in] y input to G
	/// @param[in] tn current time
	/// @param[in] Un computed value of U = [S,E,I,R] at time tn
	/// @param[in] dt timestep
	/// @param[in] y1 computed value for first intermediate RK value
	///
	void computeG2(Eigen::VectorXd &G, Eigen::VectorXd y, double tn,
				   Eigen::VectorXd Un, double dt, Eigen::VectorXd y1)
	{
		// (write your solution here)
		Eigen::VectorXd F1;
		Eigen::VectorXd F2;
		F1.resize(4);
		F2.resize(4);
		coronaOutbreak.computeF(F1, tn + mu * dt, y1);
		coronaOutbreak.computeF(F2, tn + (mu - nu) * dt, y);
		G = Un + dt * (-nu * F1 + mu * F2) - y;
	}

	///
	/// Find a solution to JG1*x = -G1 with Newton's method
	/// @param[out] x solution to the system
	/// @param[in] Un computed value of U = [S,E,I,R] at time tn
	/// @param[in] dt timestep
	/// @param[in] tn current time
	/// @param[in] tolerance if Newton increment smaller, successfully converged
	/// @param[in] maxIterations  max Newton iterations to try before failing
	///
	void newtonSolveY1(Eigen::VectorXd &u, Eigen::VectorXd Un, double dt,
					   double tn, double tolerance, int maxIterations)
	{
		int dim = Un.size();
		Eigen::VectorXd RHSG1(dim), Fu(dim), x(dim);
		Eigen::MatrixXd JG1(dim, dim), JFu(dim, dim);
		u = Un;

		// Write your loop for the Newton solver here.
		// You will need to use coronaOutbreak.computeF(...)
		// and coronaOutbreak.computeJF(...)
		// (write your solution here)

		for (int i = 0; (i < maxIterations); ++i)
		{

			coronaOutbreak.computeJF(JFu, tn + mu * dt, u);
			JG1 = dt * mu * JFu - Eigen::MatrixXd::Identity(4, 4);
			computeG1(RHSG1, u, tn, Un, dt);
			x = JG1.lu().solve(-RHSG1); //JG1*x = -G1

			if (x.norm() <= tolerance)
				return;
			u += x;
		}
		//u = x;

		// If we reach this point, something wrong happened.
		throw std::runtime_error("Did not reach tolerance in Newton iteration in Y1");
	}

	///
	/// Find a solution to JG2*x = -G2 with Newton's method
	/// @param[out] x solution to the system
	/// @param[in] Un computed value of U = [S, E, I, R] at time tn
	/// @param[in] y1 previous intermediate value for RK method
	/// @param[in] dt timestep
	/// @param[in] tn current time
	/// @param[in] tolerance if Newton increment smaller, successfully converged
	/// @param[in] maxIterations  max Newton iterations to try before failing
	///
	//----------------NewtonG2Start----------------
	void newtonSolveY2(Eigen::VectorXd &v, Eigen::VectorXd Un,
					   Eigen::VectorXd y1, double dt, double tn, double tolerance, int maxIterations)
	{

		// Use newtonSolveY1 as a model for this
		// (write your solution here)
		int dim = Un.size();
		Eigen::VectorXd RHSG2(dim), Fu(dim), x(dim);
		Eigen::MatrixXd JG2(dim, dim), JFu1(dim, dim), JFu2(dim, dim);
		v = Un;
		for (int i = 0; (i < maxIterations); ++i)
		{
			//coronaOutbreak.computeJF(JFu1, tn + mu * dt, y1);
			coronaOutbreak.computeJF(JFu2, tn + (mu - nu) * dt, v);
			JG2 = dt * mu * JFu2 - Eigen::MatrixXd::Identity(4, 4);
			computeG2(RHSG2, v, tn, Un, dt, y1);
			x = JG2.lu().solve(-RHSG2); //JG1*x = -G1

			if (x.norm() <= tolerance)
				return;
			v += x;
		}
	}
	//----------------NewtonG2End----------------

	///
	/// Compute N timesteps of DIRK(2,3)
	/// @param[in/out] u should be a vector of size 5, where each
	///                component is a vector of size N+1. u[i][0]
	///                should have the initial value stored before
	///                calling the funtion
	///
	/// @param[out] time should be of length N+1
	///
	/// @param[in] T the final time
	/// @param[in] N the number of timesteps
	///
	//----------------DirkStart----------------
	void solve(std::vector<std::vector<double>> &u, std::vector<double> &time,
			   double T, int N)
	{

		const double dt = T / N;

		// Your main loop goes here. At iteration n,
		// 1) Find Y_1 with newtonSolveY1 (resp. Y2)
		// 2) Compute U^{n+1} with F(Y1), F(Y2)
		// 3) Write the values at u[...][n]
		// 4) Compute D and write time[n]

		// (write your solution here)
		double tolerance = 1e-10;
		int iter = 100000;
		time.at(0) = 0;
		Eigen::VectorXd Y1, Y2, Un, F1, F2;
		Y1.resize(4);
		Y2.resize(4);
		Un.resize(4);
		F1.resize(4);
		F2.resize(4);
		Un(0) = u.at(0).at(0);
		Un(1) = u.at(1).at(0);
		Un(2) = u.at(2).at(0);
		Un(3) = u.at(3).at(0);
		//Un(4) = u.at(4).at(0);
		//double t = 0;
		for (int i = 1; i <= N; i++)
		{
			time.at(i) = i * dt;
			newtonSolveY1(Y1, Un, dt, time.at(i - 1), tolerance, iter);
			newtonSolveY2(Y2, Un, Y1, dt, time.at(i - 1), tolerance, iter);
			coronaOutbreak.computeF(F1, time.at(i - 1) + mu * dt, Y1);
			coronaOutbreak.computeF(F2, time.at(i - 1) + (mu - nu) * dt, Y2);
			Un += dt * ((mu - 0.5 * nu) * F1 + (mu - 0.5 * nu) * F2);
			u.at(0).at(i) = Un(0);
			u.at(1).at(i) = Un(1);
			u.at(2).at(i) = Un(2);
			u.at(3).at(i) = Un(3);
			//work in progress
		}
		coronaOutbreak.computeD(u, dt, N);
	}
	//----------------DirkEnd----------------

private:
	const double mu = 0.5 + 0.5 / sqrt(3);
	const double nu = 1. / sqrt(3);

	CoronaOutbreak coronaOutbreak;

}; // end class DIRKSolver
