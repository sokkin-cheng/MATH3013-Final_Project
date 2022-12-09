#include "solvers.h"
#include "stuff.h"
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

double PI = atan(1) * 4;

void alpha_warning(const double alpha, const double requirement) 
//checking the alpha in order to find a reliable results
{
    if (alpha > requirement)
    {
        cout << "alpha = (dt / h^2) = " << alpha << endl;
        cout << "(dt / h^2) must be < " << requirement <<" for reliable results." << endl;
        exit(1);
    }
}

void output(const std::string folder, const std::vector<double>& u, const double time_final)
//printing the output to create a .dat file 
//and use gnuplot to plot the data
{
    ofstream results(folder + "_results.dat");
    for (int i = 0; i < u.size(); i++)
    {
        results << setprecision(8) << (double)i / (u.size() - 1.) << setw(15) << u[i] << endl;
    }

    results.close();
}


void DiffusionEquation::expl(const int mpoints, const double time_final, const int time_steps, const std::string folder)
{
    //1D diffusion equation can be solves with a simple loop, as it is shown here

    const double h = 1. / (double)mpoints;          //  space-step
    const double dt = time_final / (double)time_steps;
    const double alpha = dt / (h * h);
    const double beta = 1. - 2. * alpha;  
    vector<double> u(mpoints + 1);                   //  solution vector

    alpha_warning(alpha, 0.005);  //  we require alpha > 0.005
    initial_conditions(u);

    for (double t = 0.; t <= time_final; t += dt)
    {
        for (int i = 1; i < mpoints; i++)
        {
            u[i] = alpha * (u[i + 1] + u[i - 1]) + beta * u[i];
        }
    }

    //Printing the output
    output(folder, u, time_final);
}

void DiffusionEquation::imp(const int mpoints, const double time_final, const int time_steps, const std::string folder)
{

    //Different from explicit scheme,suppose we have a linear system, we let a squared (n+1) 
    //tridiagonal matrix A with constant diagonals a, b and c.
    //Create two vectors u and y.
    //We solve A * u = u_old, u_old being u at a previous time-step with our solver

    //Remark: Implicit Scheme does not need to check the alpha

    const double h = 1. / (double)mpoints;          //  space-step
    const double dt = time_final / (double)time_steps;
    const double alpha = dt / (h * h);
    vector<double> u(mpoints + 1);                   //  solution vector
    vector<double> u_old(mpoints + 1);
    vector<double> b(mpoints + 1);                   //  vector of the main diagonal

    initial_conditions(u, u_old);

    for (double t = 0.; t < time_final; t += dt)
    {
        // Gaussian elimination, we enter the three diagonals in the following 
        tridia(mpoints, -alpha, 1. + 2 * alpha, -alpha, u, b, u_old);
        u_old = u;  //  iteration for the next time-step
    }

    //Printing outputs
    output(folder, u, time_final);
}

void DiffusionEquation::cranknicolson(const int mpoints, const double time_final, const int time_steps, const std::string folder)
{
    //Suppose we come up with a linear algebra system.
    //We first perform a matrix-vector multiplication.
    //Then we perform a matrix inversion with the new vector.
    //Finally we use our tridiagonal solver to solve it again.

    const double h = 1. / (double)mpoints;          //  space-step
    const double dt = time_final / (double)time_steps;
    const double alpha = dt / (h * h);
    const double beta = 2. - 2. * alpha;
    const double gamma = 2. + 2. * alpha;
    vector<double> u(mpoints + 1);                   //  solution vector
    vector<double> u_old(mpoints + 1);
    vector<double> b(mpoints + 1);                   //  vector of the main diagonal

    alpha_warning(alpha, 0.005);
    initial_conditions(u, u_old);

    for (double time = 0.; time < time_final; time += dt)
    {
        //  first we initialize y to use it after in the tridiag solver
        //  we compute (2I - alpha*B)*u_old and put this new vector as y
        //  but before the gaussian elimination we have u_old=u so we use u
        for (int i = 1; i < mpoints; i++)
        {
            u_old[i] = alpha * u[i - 1] + beta * u[i] + alpha * u[i + 1];
        }

        //  now we perform the gaussian elimination for (2I + alpha*B)*u = y
        tridia(mpoints, -alpha, gamma, -alpha, u, b, u_old);
        u_old = u;
    }

    // Printing output again
    output(folder, u, time_final);
}