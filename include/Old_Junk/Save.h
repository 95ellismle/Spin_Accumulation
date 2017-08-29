#ifndef SAVE_H
#define SAVE_H

#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>

#include "math_funcs.h"

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
extern const double fact3;

extern std::vector<unsigned int> steps_to_save;

// Data Storage Containers
//
// The Save Vectors Save the number of time steps requested evenly spaced.
extern std::vector < std::vector < double > > Save_Ux;
extern std::vector < std::vector < double > > Save_Uy;
extern std::vector < std::vector < double > > Save_Uz;

extern double Save_Time_Range[2];
extern std::string filepath;

std::vector<unsigned int> save_steps_calc(std::vector<unsigned int> steps_to_save, int Saves)
{
   if (Save_Time_Range[1] > T){ Save_Time_Range[1] = T;}
   if (Save_Time_Range[0] > Save_Time_Range[1]){ Save_Time_Range[0] = Save_Time_Range[1];}

   Save_Time_Range[0] = Save_Time_Range[0]/dt;
   Save_Time_Range[1] = Save_Time_Range[1]/dt;
   unsigned int counter = 0;
   if (Saves >  0)
   {
      if(almost_equal(Save_Time_Range[0],Save_Time_Range[1],1e-9))
      {
         steps_to_save.resize(2);
         steps_to_save[1] = Iterations;
         steps_to_save[0] = 0;
      }
      else
      {

         if (Saves < 2){ Saves = 2;}
         //Print_Params();
         // Deciding which steps are to be written to permanant storage
         for(unsigned int i=Save_Time_Range[0];i<Save_Time_Range[1];i++)
         {
            unsigned int x = std::ceil((Save_Time_Range[1]-Save_Time_Range[0])/(Saves));
            if (i%x == 0)
            {
               steps_to_save[counter] = i;
               counter++;
            }
         }
         steps_to_save[Saves-1] = Iterations;
      }
      for (size_t i = 0; i < steps_to_save.size(); i++)
      {
         std::cout << "steps_to_save[i] = " << steps_to_save[i] << '\n';
      }
   }
   else
   {
      steps_to_save.resize(1);
      steps_to_save[0] = 1e20;
   }
   return steps_to_save;
}


extern std::vector<double> Uz;
extern std::vector<double> Uy;
extern std::vector<double> Ux;

void Save_vec(std::vector<double> v)
{
   std::ofstream beta_file;
	std::string beta_filepath = filepath+"vec";
	beta_file.open(beta_filepath.c_str());
	beta_file << "vec\n";
	for (unsigned int i = 0; i < v.size(); i++ )
	{
		beta_file << v[i] << "\n";
	}
	beta_file.close();
}

void Save(std::string filepath)
{
	std::cout << "Just writing that data to a file\n\n";
	// Creating and Opening Files
	std::ofstream xfile;
	std::ofstream yfile;
	std::ofstream zfile;
	std::ofstream paramsfile;
	//Opening Files

	std::string x_filepath = filepath+"X_Output_Data";
	std::string y_filepath = filepath+"Y_Output_Data";
	std::string z_filepath = filepath+"Z_Output_Data";
	std::string p_filepath = filepath+"Parameters";

	xfile.open(x_filepath.c_str());
	yfile.open(y_filepath.c_str());
	zfile.open(z_filepath.c_str());
	paramsfile.open(p_filepath.c_str());
	for(unsigned int i=0;i<x_coords_len;i++)
	{
		if (i < x_coords_len - 1)
		{
			xfile <<  i << ",";
			yfile <<  i << ",";
			zfile <<  i << ",";
		}
		else
		{
			xfile <<  i <<  "\n";
			yfile <<  i <<  "\n";
			zfile <<  i <<  "\n";
		}
	}
	for (unsigned int n=0; n<1; n++)
	{
		for (unsigned int i=0;i<x_coords_len;i++)
		{
			if (i < x_coords_len - 1)
			{
				xfile <<  Ux[i] << ",";
				yfile <<  Uy[i] << ",";
				zfile <<  Uz[i] << ",";
			}
			else
			{
				xfile <<  Ux[i] <<  "\n";
				yfile <<  Uy[i] <<  "\n";
				zfile <<  Uz[i] <<  "\n";
			}
		}
	}
	//Writing the Parameter Data
	paramsfile << ":dx:" << dx;
	paramsfile << ":dt:" << dt;
	paramsfile << ":L:" << L;
	paramsfile << ":T:" << T;
	paramsfile << ":m_inf:" << minf_magnitude;
	paramsfile << ":Num_Time_Steps:" << Iterations;
	paramsfile << ":Num_Space_Steps:" << x_coords_len;
	paramsfile << ":M1:" << M1[0] << "," << M1[1] << "," << M1[2];
	paramsfile << ":M2:" << M2[0] << "," << M2[1] << "," << M2[2];
	paramsfile << ":F1_Boundaries:" << F1_Boundaries[0] << "," << F1_Boundaries[1];
	paramsfile << ":F2_Boundaries:" << F2_Boundaries[0] << "," << F2_Boundaries[1];
	paramsfile << ":F1_Boundaries_Indices:" << F1_Boundaries_Indices[0] << "," << F1_Boundaries_Indices[1];
	paramsfile << ":F2_Boundaries_Indices:" << F2_Boundaries_Indices[0] << "," << F2_Boundaries_Indices[1];
	paramsfile << ":lambda_sf:" << lambda_sf[0] << "," << lambda_sf[1];
	paramsfile << ":lambda_J:" << lambda_J[0] << "," << lambda_J[1];
	paramsfile << ":lambda_phi:" << lambda_phi[0] << "," << lambda_phi[1];
	paramsfile << ":Diffusivity:" << D[0] << "," << D[1];
	paramsfile << ":beta:" << beta[0] << "," << beta[1];
	paramsfile << ":beta_prime:" << beta_prime[0] << "," << beta_prime[1];

	//Closing Files
	xfile.close();
	yfile.close();
	zfile.close();
	paramsfile.close();
}
#endif
