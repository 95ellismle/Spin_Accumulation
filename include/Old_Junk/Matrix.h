#ifndef MATRIX_H
#define  MATRIX_H


#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include "correctors.h"



//////////////////////////////////////////////////	SYSTEM SETUP PARAMETERS	/////////////////////////////////////////////////////////////////////////
// Parameters affecting accuracy and the dimensions of the system																										     //																																																																//
extern const double dx;													// The space step (m)																			        //
extern const double dt;													// The time step (s)																				        //
extern const double T;														// The total time (s)																			     //
extern const double L;														// The total length (m)																		        //
extern const double F1_Boundaries[2];		// Position of the first ferromagnetic layer															                 //
extern const double F2_Boundaries[2];		// Position of the first ferromagnetic layer														                    //
//  																																													        //
// Material Properties. First in the array is the Magnet, the Second is in a non-magnet.																	        //
//
extern const double DL;																																											        //
extern const double j_e;																																											        //
extern const double beta[2];	 									// The charge current polarisation															           //
extern const double beta_prime[2]; 						// The spin current polarisation																		           //
extern const double lambda_sf[2];					// The spin flip length																					              //
extern const double lambda_J[2];							// A damping coefficient (spin accumulation precession around Magnetisation)                   //
extern const double lambda_phi[2];						// Another damping coefficient (a dephasing lengthscale)										           //
extern const double D[2];								// Diffusivity 1st item is for the magnet and 2nd is for the non-magnet			                    //
extern const double minf_magnitude;						// The magnitude of the equilibrium spin accumlation vector 							              //
//  																																													        //
// The magnetisations of the layers																																			        //
//  																																													        //
extern const double M1[3];							// F1 Magnetisation, this code will only work with an x and y magnetisation 		     			        //
extern const double M2[3];										// F2 Magnetisation, this code will only work with an x and y magnetisation	 		        //
//  																																													 	     //
// The Number of Time Steps to Save																																			 	     //
//  																																													 	     //
extern unsigned int Saves;														// F1 Magnetisation, this code will only work with an x and y magnetisation     //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// The lengths of the data storage containers and the number of iterations.
extern const unsigned int x_coords_len;
extern const unsigned int Iterations;

// The Positions of the Magnetic Layers in the Storage Containers.
extern const unsigned int F1_Boundaries_Indices[2];
extern const unsigned int F2_Boundaries_Indices[2];

// Constants to use later
extern const double fact1;
extern const double fact2;
extern const double fact3[2];

// Vectors
extern double Byy;
extern double Bzz;
extern double B1yz;
extern double B2yz;
extern std::vector<double> vec_beta;
extern std::vector<double> diff_beta;
extern std::vector<double> Sx;
extern std::vector<double> v;
extern std::vector<double> S2x;
extern std::vector<double> Cy;
extern std::vector<double> Cz;
extern std::vector<unsigned int> steps_to_save;

// Data Storage Containers
//
// The Save Vectors Save the number of time steps requested evenly spaced.
extern std::vector < std::vector < double > > Save_Ux;
extern std::vector < std::vector < double > > Save_Uy;
extern std::vector < std::vector < double > > Save_Uz;
//
// These Vectors are used in the calculations
extern std::vector<double> Uz;
extern std::vector<double> Uz_half;
extern std::vector<double> Uy;
extern std::vector<double> Uy_half;
extern std::vector<double> Ux;

// The A Matrix
extern std::vector<double> A;






// A function to recreate the A Matrix
void A_Refresher(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> L)
{
	for (unsigned int i=3; i<(x_coords_len-1)*3;i+=3)
	{
		unsigned int x = (int) i/3;
		A[i] 		= -a[x];
		A[i+1] 	= 1+(2*b[x])+L[x];
		A[i+2]	= -c[x]+L[x];
	}
	A[0] = 0;
  A[1] = 1;
  A[2] = 0;
  A[3*x_coords_len-1] = 0;
  A[3*x_coords_len-2] = 1;
  A[3*x_coords_len-3] = 0;
}

// A function to solve the Matrix equation
std::vector<double> Solver(std::vector<double> & input_array, int xyz)
{
	std::vector <double> tmp_array (input_array);
	if (xyz == 0)
	{
		for(unsigned int i=1; i<F1_Boundaries_Indices[0]; i++)
		{
			tmp_array[i] = tmp_array[i] + Z_NM_corrector(i);
		}
		for (unsigned int i=F1_Boundaries_Indices[0];i<F1_Boundaries_Indices[1]; i++)
		{
			tmp_array[i] = tmp_array[i] + Z1_corrector(i);
		}
		for(unsigned int i=F1_Boundaries_Indices[1]; i<F2_Boundaries_Indices[0]; i++)
		{
			tmp_array[i] = tmp_array[i] + Z_NM_corrector(i);
		}
		for (unsigned int i=F2_Boundaries_Indices[0];i<F2_Boundaries_Indices[1]; i++)
		{
			tmp_array[i] = tmp_array[i] + Z2_corrector(i);
		}
		for(unsigned int i=F2_Boundaries_Indices[1]; i<x_coords_len-1; i++)
		{
			tmp_array[i] = tmp_array[i] + Z_NM_corrector(i);
		}
	}
	else if (xyz == 1)
	{
		for(unsigned int i=1; i<F1_Boundaries_Indices[0]; i++)
		{
			tmp_array[i] = tmp_array[i] + Z_NM_corrector(i);
		}
		for (unsigned int i=F1_Boundaries_Indices[0];i<F1_Boundaries_Indices[1]; i++)
		{
			tmp_array[i] = tmp_array[i] + Z1_half_corrector(i);
		}
		for(unsigned int i=F1_Boundaries_Indices[1]; i<F2_Boundaries_Indices[0]; i++)
		{
			tmp_array[i] = tmp_array[i] + Z_NM_corrector(i);
		}
		for (unsigned int i=F2_Boundaries_Indices[0];i<F2_Boundaries_Indices[1]; i++)
		{
			tmp_array[i] = tmp_array[i] + Z2_half_corrector(i);
		}
		for(unsigned int i=F2_Boundaries_Indices[1]; i<x_coords_len-1; i++)
		{
			tmp_array[i] = tmp_array[i] + Z_NM_corrector(i);
		}
	}
	else if (xyz == 2)
	{
		for(unsigned int i=1; i<F1_Boundaries_Indices[0]; i++)
		{
			tmp_array[i] = tmp_array[i] + Y_NM_corrector(i);
		}
		for (unsigned int i=F1_Boundaries_Indices[0];i<F1_Boundaries_Indices[1]; i++)
		{
			tmp_array[i] = tmp_array[i] + Y1_corrector(i);
		}
		for(unsigned int i=F1_Boundaries_Indices[1]; i<F2_Boundaries_Indices[0]; i++)
		{
			tmp_array[i] = tmp_array[i] + Y_NM_corrector(i);
		}
		for (unsigned int i=F2_Boundaries_Indices[0];i<F2_Boundaries_Indices[1]; i++)
		{
			tmp_array[i] = tmp_array[i] + Y2_corrector(i);
		}
		for(unsigned int i=F2_Boundaries_Indices[1]; i<x_coords_len-1; i++)
		{
			tmp_array[i] = tmp_array[i] + Y_NM_corrector(i);
		}
	}
	else if (xyz == 3)
	{
		for(unsigned int i=1; i<F1_Boundaries_Indices[0]; i++)
		{
			tmp_array[i] = tmp_array[i] + Y_NM_corrector(i);
		}
		for (unsigned int i=F1_Boundaries_Indices[0];i<F1_Boundaries_Indices[1]; i++)
		{
			tmp_array[i] = tmp_array[i] + Y1_half_corrector(i);
		}
		for(unsigned int i=F1_Boundaries_Indices[1]; i<F2_Boundaries_Indices[0]; i++)
		{
			tmp_array[i] = tmp_array[i] + Y_NM_corrector(i);
		}
		for (unsigned int i=F2_Boundaries_Indices[0];i<F2_Boundaries_Indices[1]; i++)
		{
			tmp_array[i] = tmp_array[i] + Y2_half_corrector(i);
		}
		for(unsigned int i=F2_Boundaries_Indices[1]; i<x_coords_len-1; i++)
		{
			tmp_array[i] = tmp_array[i] + Y_NM_corrector(i);
		}
	}
	else if (xyz == 4)
	{
		for(unsigned int i=1; i<F1_Boundaries_Indices[0]; i++)
		{
			tmp_array[i] = tmp_array[i] + X_NM_corrector(i);
		}
		for (unsigned int i=F1_Boundaries_Indices[0];i<F1_Boundaries_Indices[1]; i++)
		{
			tmp_array[i] = tmp_array[i] + X1_corrector(i);
		}
		for(unsigned int i=F1_Boundaries_Indices[1]; i<F2_Boundaries_Indices[0]; i++)
		{
			tmp_array[i] = tmp_array[i] + X_NM_corrector(i);
		}
		for (unsigned int i=F2_Boundaries_Indices[0];i<F2_Boundaries_Indices[1]; i++)
		{
			tmp_array[i] = tmp_array[i] + X2_corrector(i);
		}
		for(unsigned int i=F2_Boundaries_Indices[1]; i<x_coords_len-1; i++)
		{
			tmp_array[i] = tmp_array[i] + X_NM_corrector(i);
		}
	}
	else if (xyz == 5)
	{
		for(unsigned int i=1; i<F1_Boundaries_Indices[0]; i++)
		{
			tmp_array[i] = tmp_array[i] + X_NM_corrector(i);
		}
		for (unsigned int i=F1_Boundaries_Indices[0];i<F1_Boundaries_Indices[1]; i++)
		{
			tmp_array[i] = tmp_array[i] + X1_half_corrector(i);
		}
		for(unsigned int i=F1_Boundaries_Indices[1]; i<F2_Boundaries_Indices[0]; i++)
		{
			tmp_array[i] = tmp_array[i] + X_NM_corrector(i);
		}
		for (unsigned int i=F2_Boundaries_Indices[0];i<F2_Boundaries_Indices[1]; i++)
		{
			tmp_array[i] = tmp_array[i] + X2_half_corrector(i);
		}
		for(unsigned int i=F2_Boundaries_Indices[1]; i<x_coords_len-1; i++)
		{
			tmp_array[i] = tmp_array[i] + X_NM_corrector(i);
		}
	}

	int counter = 1;
	for(unsigned int i=3; i<3*(x_coords_len-1);i+=3)
	{
		double mul = A[i]/A[i-2];
		double new_value = (double) tmp_array[counter] - tmp_array[counter-1]*mul;
		tmp_array[counter] = new_value;
		A[i+1] = A[i+1] - A[i-1]*mul;
		A[i] = A[i] - A[i-2]*mul;
		counter++;
	}

	counter = x_coords_len-1;
	for(int i=3*(x_coords_len-1); i>3; i-=3)
	{
		double mul = A[i-1]/A[i+1];
		double new_value = (double) tmp_array[counter-1] - tmp_array[counter]*mul;
		tmp_array[counter-1] = new_value;
		A[i-1] = A[i-1] - A[i+1]*mul;
		counter--;
	}

	for(unsigned int i=1; i<x_coords_len-1; i++)
	{
		double new_value = (double) tmp_array[i]/A[3*i+1];
		tmp_array[i] = new_value;
	}
	//Neumann (floating) boundary conditions
	tmp_array[0] = tmp_array[1] - (tmp_array[2]-tmp_array[1])/2.718;
	tmp_array[x_coords_len-1] = tmp_array[x_coords_len-2]-(tmp_array[x_coords_len-3]-tmp_array[x_coords_len-2])/2.718;
	return tmp_array;
}



#endif
