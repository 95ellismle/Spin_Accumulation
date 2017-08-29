#ifndef DEBUG_H
#define  DEBUG_H

#include <iostream>
#include <vector>
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




// // The A Matrix
// extern std::vector<double> A;
//
//
// // A function to check the A Matrix
// // This can be useful for debugging.
// void A_Printer()
// {
// 	for (unsigned int i=0; i<(x_coords_len)*3;i+=3)
// 	{
// 		std::cout << A[i]<<","<<A[i+1]<<","<<A[i+2] << std::endl;
// 	}
// }

// A function to print out the parameters to the console
// This is very useful for debugging!
void Print_Params()
{
	std::cout << "\n############Parameters############\n\n"<< std::endl;

	std::cout << "Size Parameters:\n";
	std::cout << "L = " << L << std::endl;
	std::cout << "T = " << T << std::endl;
	std::cout << "dx = " << dx << std::endl;
	std::cout << "dt = " << dt << std::endl;
	std::cout << "F1_Boundary Indexes = " << F1_Boundaries_Indices[0] << "," << F1_Boundaries_Indices[1] << std::endl;
	std::cout << "F2_Boundary_Indices = " << F2_Boundaries_Indices[0] << "," << F2_Boundaries_Indices[1] << std::endl;

	std::cout << "\nMagnet Parameters:\n";
	std::cout << "Beta = " << beta[0] << std::endl;
	std::cout << "Beta' = " << beta_prime[0] << std::endl;
	std::cout << "D = " << D[0] << std::endl;
	std::cout << "lambda_sf = " << lambda_sf[0] << std::endl;
	std::cout << "lambda_phi = " << lambda_phi[0] << std::endl;
	std::cout << "lambda_J  = " << lambda_J[0] << std::endl;

	std::cout << "\nNon-Magnet Parameters:\n";
	std::cout << "Beta = " << beta[1] << std::endl;
	std::cout << "Beta' = " << beta_prime[1] << std::endl;
	std::cout << "D = " << D[1] << std::endl;
	std::cout << "lambda_sf = " << lambda_sf[1] << std::endl;
	std::cout << "lambda_phi = " << lambda_phi[1] << std::endl;
	std::cout << "lambda_J  = " << lambda_J[1] << std::endl;

	std::cout << "\nOptimising Parameters:\nParam\tNon-mag\tMag\n";
	std::cout << "\n#################################\n"<< std::endl;
}

#endif
