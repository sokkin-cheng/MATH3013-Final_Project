#ifndef SOLVERS_H
#define SOLVERS_H
#pragma once
#include <string>

class DiffusionEquation {
public:
	void cranknicolson(const int meshpoints, const double time_final, const int time_steps, const std::string folder);
	void expl(const int meshpoints, const double time_final, const int time_steps, const std::string folder);
	void imp(const int meshpoints, const double time_final, const int time_steps, const std::string folder);
};
#endif 
