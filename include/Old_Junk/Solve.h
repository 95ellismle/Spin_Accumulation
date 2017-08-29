#ifndef SOLVE_H
#define SOLVE_H

#include <vector>
#include <iostream>

#include "math_funcs.h"

extern const unsigned int x_coords_len;
extern const double dt;
extern const double minf_magnitude;

//
// These Vectors are used in the calculations
extern std::vector<double> Uz;
extern std::vector<double> Uy;
extern std::vector<double> Ux;
extern std::vector<double> Ux_prev;
extern std::vector<double> Uy_prev;
extern std::vector<double> Uz_prev;

extern std::vector<double> vec_beta;
extern std::vector<double> vec_Mx;
extern std::vector<double> vec_My;
extern std::vector<double> vec_Mz;
extern std::vector<double> vec_D;
extern std::vector<double> vec_lambda_J;
extern std::vector<double> vec_lambda_phi;
extern std::vector<double> vec_lambda_sf;
extern std::vector<double> vec_beta_prime;

std::vector<double> SOLVEX(std::vector<double> Ux, std::vector<double> Jmx)
{
	std::vector<double> diff_Jmx (x_coords_len);
	diff_Jmx = Differentiator(Jmx);
	for (unsigned int i =0; i<x_coords_len; i++)
	{
		double fact1 = dt/(vec_lambda_J[i]*vec_lambda_J[i]);
		double fact2 = dt/(vec_lambda_phi[i]*vec_lambda_phi[i]);
		double first_term = diff_Jmx[i];
		double second_term = fact1*(-vec_My[i]*Uz_prev[i] + vec_Mz[i]*Uy_prev[i] );
		double third_term = fact2*(vec_My[i]*vec_My[i]*Ux_prev[i] + vec_Mz[i]*vec_Mz[i]*Ux_prev[i]);
		double fourth_term = fact3*(Ux_prev[i] - minf_magnitude*vec_Mx[i]);
		Ux[i] = Ux_prev[i] - (first_term + second_term + third_term + fourth_term);
	}
	return Ux;
}
std::vector<double> SOLVEY(std::vector<double> Uy, std::vector<double> Jmy)
{
	std::vector<double> diff_Jmy (x_coords_len);
	diff_Jmy = Differentiator(Jmy);
	for (unsigned int i =0; i<x_coords_len; i++)
	{
		double fact1 = dt/(vec_lambda_J[i]*vec_lambda_J[i]);
		double fact2 = dt/(vec_lambda_phi[i]*vec_lambda_phi[i]);
		double fact3 = dt/(vec_lambda_sf[i]*vec_lambda_sf[i]);
		double first_term = diff_Jmy[i]*dt;
		double second_term = fact1*(-vec_Mz[i]*Ux_prev[i]);
		double third_term = fact2*(-vec_My[i]*vec_Mz[i]*Uz_prev[i] + vec_Mz[i]*vec_Mz[i]*Uy_prev[i]);
		double fourth_term = fact3*(Uy_prev[i] - minf_magnitude*vec_My[i]);
		Uy[i] = Uy_prev[i] - (first_term +second_term + third_term + fourth_term);
	}
	return Uy;
}
std::vector<double> SOLVEZ(std::vector<double> Uz, std::vector<double> Jmz)
{
	std::vector<double> diff_Jmz (x_coords_len);
	diff_Jmz = Differentiator(Jmz);
	for (unsigned int i =0; i<x_coords_len; i++)
	{
		double fact1 = dt/(vec_lambda_J[i]*vec_lambda_J[i]);
		double fact2 = dt/(vec_lambda_phi[i]*vec_lambda_phi[i]);
		double fact3 = dt/(vec_lambda_sf[i]*vec_lambda_sf[i]);
		double first_term = diff_Jmz[i]*dt;
		double second_term = fact1*(vec_My[i]*Ux_prev[i]);
		double third_term = fact2*(vec_My[i]*vec_My[i]*Uz_prev[i] - vec_My[i]*vec_Mz[i]*Uy_prev[i]);
		double fourth_term = fact3*(Uz_prev[i] - minf_magnitude*vec_Mz[i]);
		Uz[i] = Uz_prev[i] - (first_term +second_term + third_term + fourth_term);
	}
	return Uz;
}


// void SOLVE(int which_code)
// {
// 	if (which_code == 0)
// 	{
// 		// First Step, Finding Uz_half
// 		A_Refresher(Cz, Cz, Cz, v);
// 		Uz_half = Solver(Uz,0);
// 		// Second Step, Finding Uy_half
// 		A_Refresher(Cy, Cy, Cy, v);
// 		Uy_half = Solver(Uy,2);
// 		// Third Step, Finding Ux_half
// 		A_Refresher(S2x, S2x, S2x, v);
// 		Ux = Solver(Ux,4);
// 		// Fourth Step, Finding Uz
// 		A_Refresher(Cz, Cz, Cz, v);
// 		Uz = Solver(Uz_half,1);
// 		// Fifth Step, Finding Uy
// 		A_Refresher(Cy, Cy, Cy, v);
// 		Uy = Solver(Uy_half,3);
// 		// Sixth Step, Finding Ux
// 		A_Refresher(S2x, S2x, S2x, v);
// 		Ux = Solver(Ux,5);
// 	}
// 	// else if (which_code == 1)
// 	// {
// 	// 	// First Step, Finding Uz_half
// 	// 	A_Refresher(Cz, Cz, Cz, Lz);
// 	// 	Uz_half = Solver(Uz,0);
// 	// 	// Second Step, Finding Uy_half
// 	// 	A_Refresher(Cy, Cy, Cy, Ly);
// 	// 	Uy_half = Solver(Uy,2);
// 	// 	// Third Step, Finding Ux_half
// 	// 	A_Refresher(S2x, S2x, S2x, v);
// 	// 	Ux = Solver(Ux,4);
// 	// 	// Fourth Step, Finding Uz
// 	// 	A_Refresher(Cz, Cz, Cz, Lz);
// 	// 	Uz = Solver(Uz_half,1);
// 	// 	// Fifth Step, Finding Uy
// 	// 	A_Refresher(Cy, Cy, Cy, Ly);
// 	// 	Uy = Solver(Uy_half,3);
// 	// 	// Sixth Step, Finding Ux
// 	// 	A_Refresher(S2x, S2x, S2x, v);
// 	// 	Ux = Solver(Ux,5);
// 	// }
// }

#endif
