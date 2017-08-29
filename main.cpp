// Including standard libraries
#include <vector>
#include <iostream>
#include <string>
#include <string.h>

// Including my own libraries
#include "vec_fill.h"
#include "math_funcs.h"
#include "save.h"
#include "test.h"
#include "solve.h"
#include "parse_and_string_stuff.h"




//////////////////////////////////////////////////////////////
////////////////               PARAMETERS              ///////
double dt = 6e-21;														// Time step
const double dx = 2e-10;												// Space Step
const double L  = 150e-9;												// Simulation size
const double T = 600e-21;												// Simulation run time
																            //
const double boundaries[4] = {45e-9, 73e-9, 78e-9, 84e-9};  // Boundaries of magnetic layers (start of F1, end of F1, start of F2, end of F2)
unsigned int Saves = 1;                                     // How many steps to save (keep this as 1 -not yet implemented)
const double DL = 3e-9;                                     // Diffusiveness of interfaces
const double minf = 4e7;                                    // Equilibrium spin accumulation
const double je_peak = 1e11;                                // Peak of the ramp signal for charge Current
const double charge_current_teq = 1e-50;                    // equilibrium time until charge current reaches its peak
// The rest of the parameters are the peaks of each 			//
// The order is (non-magnet, F1, F2)								//
const double beta_peaks[3] = {0, 0.5, 0.5};                 // Peaks of beta
const double beta_prime_peaks[3] = {0, 0.7, 0.7};           // Peaks of beta prime
const double D_peaks[3] = {0.005, 0.003, 0.003};            // Peaks of D0
const double lambda_sf_peaks[3] = {20e-9, 5e-9, 5e-9};      // Peaks of lambda_sf
const double lambda_J_peaks[3] = {1e-20, 4e-9, 4e-9};       // Peaks of lambda_J
const double lambda_phi_peaks[3] = {1e-20, 4e-9, 4e-9};     // Peaks of lambda_phi
const double Mx_peaks[3] = {0, 0, 0};                  		// Peaks of X-component of magnetisation
const double My_peaks[3] = {0, 0, 1};                       // Peaks of Y-component of magnetisation
const double Mz_peaks[3] = {0, 1, 0};                       // Peaks of Z-component of magnetisation
//////////////////////////////////////////////////////////////


std::vector<double> x = range(0.0, L, dx);    // holds the x coordinates
const unsigned int xlen = x.size();           // total number of x coords (for use in loops)



/////// These vectors hold the diffuse parameters /////////////////////////////////////////////////////////////////////
/////// The function diffuse calc calculates the diffusive parameters /////////////////////////////////////////////////
/////// If the solution of Fick's law needs to be implemented this is where it should be done. 	///////////////////////
std::vector<double> beta = diffuse_calc(x, beta_peaks[0], beta_peaks[1], beta_peaks[2]);										//
std::vector<double> beta_prime = diffuse_calc(x, beta_prime_peaks[0], beta_prime_peaks[1], beta_prime_peaks[2]);		//
std::vector<double> D = diffuse_calc(x, D_peaks[0], D_peaks[1], D_peaks[2]);														//
std::vector<double> lambda_sf = diffuse_calc(x, lambda_sf_peaks[0], lambda_sf_peaks[1], lambda_sf_peaks[2]);			//
std::vector<double> lambda_J = diffuse_calc(x, lambda_J_peaks[0], lambda_J_peaks[1], lambda_J_peaks[2]);					//
std::vector<double> lambda_phi = diffuse_calc(x, lambda_phi_peaks[0], lambda_phi_peaks[1], lambda_phi_peaks[2]);     //
std::vector<double> Mx = diffuse_calc(x, Mx_peaks[0], Mx_peaks[1], Mx_peaks[2]);													//
std::vector<double> My = diffuse_calc(x, My_peaks[0], My_peaks[1], My_peaks[2]);													//
std::vector<double> Mz = diffuse_calc(x, Mz_peaks[0], Mz_peaks[1], Mz_peaks[2]);													//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// These are vectors to hold the dt/lambda factors.
std::vector<double> fact1 (xlen);
std::vector<double> fact2 (xlen);
std::vector<double> fact3 (xlen);

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// Data Storage Containers /////////////////////////////////////////////////////////////////
//																																	//
// The Save Vectors Save the number of time steps requested evenly spaced. (not yet implemented)   //
std::vector < std::vector < double > > Save_Ux (Saves, std::vector< double >(xlen,0));					//
std::vector < std::vector < double > > Save_Uy (Saves, std::vector< double >(xlen,0));					//
std::vector < std::vector < double > > Save_Uz (Saves, std::vector< double >(xlen,0));					//
std::vector < std::vector < double > > Save_Jmx (Saves, std::vector< double >(xlen,0));			   //
std::vector < std::vector < double > > Save_Jmy (Saves, std::vector< double >(xlen,0));				//
std::vector < std::vector < double > > Save_Jmz (Saves, std::vector< double >(xlen,0));				//
//																																	//
// These Vectors are used in the calculations																		//
std::vector<double> Ux (xlen, 0);																						//
std::vector<double> Uy (xlen, 0);																						//
std::vector<double> Uz (xlen, 0);																						//
/////// The Ui_new vectors are temporary storage containers to put the new data in after its     	//
/////// calculated. These are needed as the method requires the data from the previous time-step   //
/////// to calculate the next one and that would be overwritten without Ui_new...					   //
std::vector<double> Ux_new (xlen, 0);																					//
std::vector<double> Uy_new (xlen, 0);																					//
std::vector<double> Uz_new (xlen, 0);																					//
// Spin Current storage Containers																						//
std::vector<double> jmx (xlen, 0);																						//
std::vector<double> jmy (xlen, 0);																						//
std::vector<double> jmz (xlen, 0);																						//
// Stores the gradient of spin current 																				//
std::vector<double> grad_jmx (xlen, 0);																				//
std::vector<double> grad_jmy (xlen, 0);																				//
std::vector<double> grad_jmz (xlen, 0);																				//
/////////////////////////////////////////////////////////////////////////////////////////////////////

std::string filepath = "../../"; // Path to save the data to.


unsigned int iterations; // How many iterations will be carried out, this is normally determined by setting the T variable. However, a comand line argument can be used. E.g. -i=100.
int main(int argc, char* argv[])
{
	// Just prints a warning if the time step is too big.
	if ( (dt > 1.90980298e+02*(dx*dx) + (1.61831964e-09*dx) -4.90017167e-20) and (1.90980298e+02*(dx*dx) + (1.61831964e-09*dx) -4.90017167e-20) > 0) {
		std::cout << "---------------WARNING-----------------\n\nAre you sure you want to use a time step of "<< dt <<" with a space step of "<< dx <<
		" this will probably result in unstable code!"<< '\n';
		std::cout << "I recommend using a smaller time step. Something less than " << 1.90980298e+02*(dx*dx) + (1.61831964e-09*dx) -4.90017167e-20 << '\n';
		std::cout << "Please enter a more sensible time step: " << '\n';
		std::cin >> dt;
	}

	// Parsing the command line arguments, This doesn't affect how things run and is not important. Most of this can be removed without changing anything.
	std::string arg_buff = "i=";
	int temp = cmd_parse(arg_buff, argc, argv);
	if (temp>0) iterations = temp;
	else iterations = std::ceil(T/dt); // If there is no command line argument set iterations to stop at the desired time, T.


	// Filling the fact1, fact2 and fact3 vectors. (see the maths)
	for (unsigned int i=0; i<xlen; i++)
	{
		fact1[i] = dt/(lambda_J[i]*lambda_J[i]);
		fact2[i] = dt/(lambda_phi[i]*lambda_phi[i]);
		fact3[i] = dt/(lambda_sf[i]*lambda_sf[i]);
	}

	// Run the solver for the set number of iterations.
	for (unsigned int i=1; i<=iterations; i++)
	{
		std::cout << "iteration " << i << "/" << iterations << '\r';
		Solver(Ux, Uy, Uz, jmx, jmy, jmz);
	}

	// These fill the first index of the Save_Ux vector. This Vector holds all the steps that should be saved.
	// For example, if 100 steps have been completed and 2 steps should be saved the first and last steps are saved.
	// If it is set to 3 then 0, 50 and 100 will be saved. (Not yet implemented).
	Save_Ux[0] = Ux;
	Save_Uy[0] = Uy;
	Save_Uz[0] = Uz;
	Save_Jmx[0] = jmx;
	Save_Jmy[0] = jmy;
	Save_Jmz[0] = jmz;

	//std::cout << "\n" << '\n' << std::endl;
	//is_valid(Uz);


	Save_Params_V(); // Saves the vector parameters (beta, beta_prime etc...)
	Save_Params_S(); // Saves scalar parameters (beta_peaks, beta_prime_peaks etc...)
	//Save_Vector(jmz);// Will save vector to a file to be plotted.
	Save_Data(); // Will save the Ux, Uy, Uz, Jmx, Jmy and Jmz data to a csv file.
	return 0;
}
