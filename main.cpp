#include "solvers.h"
#include <string>

using namespace std;

int main(int argc, const char* argv[])
{
    //Setup
    const string diff_exp = "Explicit";
    const string diff_imp = "Implicit";
    const string diff_cn = "CrankNicolson";
    const int mpoints = 100;
    const double time_final = 0.02;
    const int time_steps = 1e5;
    
    //Solving Diffusion Equaiton with 3 schemes
    DiffusionEquation de;
    de.expl(100, 0.02, 1e5, diff_exp);
    de.imp(100, 0.02, 1e5, diff_imp);
    de.cranknicolson(100, 0.02, 1e5, diff_cn);
    return 0;
}