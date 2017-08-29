#ifndef FILL_H
#define FILL_H

#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>

#include "math_funcs.h"


extern std::vector<double> vec_beta;
extern std::vector<double> vec_My;
extern std::vector<double> vec_Mz;
extern std::vector<double> vec_D;
extern std::vector<double> vec_beta_prime;

// A function to fill a vector with the beta values (to make the beta act as if it is changing as a function of x)
std::vector<double> tanh_fill(std::vector<double> vec, double max1, double max2=false, double zero=0)
{
	if (max2 == false) max2=max1;
	if (zero)
	{
		max1 = max1-zero;
		max2 = max1;
	}
	double mid_points[3] = {F1_Boundaries[1]-F1_Boundaries[0], F2_Boundaries[0]-F1_Boundaries[1], F2_Boundaries[1]-F2_Boundaries[0]};
	for( unsigned int i =0; i<x_coords_len; i++)
	{
		vec[i] = max1/2*( std::tanh( 2*( i*dx - F1_Boundaries[0]) / DL ) + 1 );
	}
	for( unsigned int i =0; i<x_coords_len; i++)
	{
		vec[i] = vec[i] - max1/2*( std::tanh( 2*( i*dx - F1_Boundaries[0] - mid_points[0] ) / DL ) + 1 ) + max1;
	}
	for( unsigned int i =0; i<x_coords_len; i++)
	{
		vec[i] = vec[i] + max2/2*( std::tanh( 2*( i*dx - F1_Boundaries[1] - mid_points[1] ) / DL ) + 1 );
	}
	for( unsigned int i =0; i<x_coords_len; i++)
	{
		vec[i] = (vec[i] - max2/2*( std::tanh( 2*( i*dx - F2_Boundaries[0] - mid_points[2] ) / DL ) + 1 ) + max2);
	}
	double min_element = vec[0];
	for( unsigned int i=0; i< x_coords_len; i++)
	{
		vec[i] = vec[i] - min_element + zero;
	}
	return vec;
}

std::vector<double> Jm (std::vector<double> vec, std::vector<double> M, std::vector<double> U)
{
	std::vector<double> diff_U (x_coords_len);
	diff_U = Differentiator(U);
	for(unsigned int i =0; i<x_coords_len; i++)
	{
		double first_term = vec_beta[i]*M[i]*j_e;
		double second_term_without_2D = diff_U[i] - (vec_beta[i]*vec_beta_prime[i]*M[i]*M[i]*diff_U[i]);
		vec[i] = first_term - (2*vec_D[i]*second_term_without_2D);
	}
	return vec;
}



#endif
