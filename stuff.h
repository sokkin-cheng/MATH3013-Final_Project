#ifndef STUFF_h
#define STUFF_h
#pragma once

#include <vector>
#include <string>

void initial_conditions(std::vector<double>& u);
void initial_conditions(std::vector<double>& u, std::vector<double>& y);
void tridia(const int n, const double a, const double b_val,
    const double c, std::vector<double>& u, std::vector<double>& b, std::vector<double>& u_old);
void output(const std::string folder, const std::vector<double>& u, const double time_final);

void initial_conditions(std::vector<double>& u)
{
    //Initial conditions for the one-dimension system

    fill(u.begin(), u.end() - 1, 0);
    u.back() = 0.1;  //Suppose the u[t_end] is 0.1 which is the same for the following initial_conditions
}

void initial_conditions(std::vector<double>& u, std::vector<double>& u_old)
{
    //Initial conditions for the one-dimension system

    fill(u.begin(), u.end() - 1, 0);
    u.back() = 0.1;
    u_old = u;
}

void tridia(const int meshpoints, const double a,
    const double b_val, const double c, std::vector<double>& u,
    std::vector<double>& b, std::vector<double>& u_old)
{
    //  reset the matrix at each loop, for every scheme we need a power of the same matrix

    fill(b.begin(), b.end(), b_val);

    //  forward substitution for A*u = u_old
    //  remember that u(k) in the code is u(k-1) irl
    //  therefore we stop at u(n-1) which is u(n) irl
    //  we do not ever work on the last boundary point of the vector

    for (int i = 1; i < meshpoints; i++)
    {
        b[i] -= (a / b[i - 1]) * c;
        u_old[i] -= (a / b[i - 1]) * u_old[i - 1];
    }

    //  backward substitution

    for (int i = meshpoints - 1; i > 0; i--)
    {
        u[i] = (u_old[i] - c * u[i + 1]) / b[i];
    }
}
#endif