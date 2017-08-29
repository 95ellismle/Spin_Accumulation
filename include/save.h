#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#ifndef SAVE_H
#define SAVE_H

// Importing external variables
extern std::string filepath;
extern const unsigned int xlen;
extern unsigned int Saves;

extern double dt;
extern const double dx;
extern const double L ;
extern const double T;
extern const double DL;
extern const double beta_peaks[3];
extern const double beta_prime_peaks[3];
extern const double D_peaks[3];
extern const double lambda_sf_peaks[3];
extern const double lambda_J_peaks[3];
extern const double lambda_phi_peaks[3];
extern const double Mx_peaks[3];
extern const double My_peaks[3];
extern const double Mz_peaks[3];
extern const double minf;
extern const double je_peak;
extern const double charge_current_teq;


extern std::vector<double> beta;
extern std::vector<double> beta_prime;
extern std::vector<double> D;
extern std::vector<double> lambda_J;
extern std::vector<double> lambda_sf;
extern std::vector<double> lambda_phi;
extern std::vector<double> Mx;
extern std::vector<double> My;
extern std::vector<double> Mz;
extern std::vector<double> x;

extern std::vector< std::vector<double> > Save_Ux;
extern std::vector< std::vector<double> > Save_Uy;
extern std::vector< std::vector<double> > Save_Uz;
extern std::vector< std::vector<double> > Save_Jmx;
extern std::vector< std::vector<double> > Save_Jmy;
extern std::vector< std::vector<double> > Save_Jmz;
////////////////////////////////////////////////////////////////



// A function to save a single vector in a specified location
void Save_Vector(std::vector<double> v, std::string path = filepath + "vec")
{
   std::ofstream vec_file;
	vec_file.open(path.c_str());
	vec_file << "vec\n";
	for (unsigned int i = 0; i < v.size(); i++ )
	{
		vec_file << v[i] << "\n";
	}
	vec_file.close();
}

// A function to save the parameters as a csv, scalars are saved in only the first row.
void Save_Params_V(std::string path = filepath + "params_vec")
{
   std::ofstream paramsfile;
   paramsfile.open(path.c_str());
   paramsfile << "Beta,Beta_Prime,D,Lambda_Sf,Lambda_J,Lambda_Phi,Mx,My,Mz,x\n";
   for (size_t i = 0; i < xlen; i++) {
      paramsfile << beta[i] << "," << beta_prime[i] << "," << D[i] << "," << lambda_sf[i] << "," << lambda_J[i] << "," << lambda_phi[i] << "," << Mx[i] << "," << My[i] << "," << Mz[i] << "," << x[i] << "\n";
   }
   paramsfile.close();
}

// Save the scalar parameters such as the peaks in the beta, beta_prime etc...
void Save_Params_S(std::string path = filepath + "params_scalar")
{
   std::ofstream paramsfile;
   paramsfile.open(path.c_str());
   paramsfile << "Beta,Beta_Prime,D,Lambda_Sf,Lambda_J,Lambda_Phi,Mx,My,Mz,dx,dt,je,Dl,T,L\n";
   paramsfile << beta_peaks[0] << ","
              << beta_prime_peaks[0] << ","
              << D_peaks[0] << ","
              << lambda_sf_peaks[0] << ","
              << lambda_J_peaks[0] << ","
              << lambda_phi_peaks[0] << ","
              << Mx_peaks[0] << ","
              << My_peaks[0] << ","
              << Mz_peaks[0] << ","
              << dx << ","
              << dt << ","
              << je_peak << ","
              << DL << ","
              << T << ","
              << L << "\n";
paramsfile << beta_peaks[1] << ","
          << beta_prime_peaks[1] << ","
          << D_peaks[1] << ","
          << lambda_sf_peaks[1] << ","
          << lambda_J_peaks[1] << ","
          << lambda_phi_peaks[1] << ","
          << Mx_peaks[1] << ","
          << My_peaks[1] << ","
          << Mz_peaks[1] << "\n";
paramsfile << beta_peaks[2] << ","
          << beta_prime_peaks[2] << ","
          << D_peaks[2] << ","
          << lambda_sf_peaks[2] << ","
          << lambda_J_peaks[2] << ","
          << lambda_phi_peaks[2] << ","
          << Mx_peaks[2] << ","
          << My_peaks[2] << ","
          << Mz_peaks[2] << "\n";
}

// A function to Save the Data in separate files for the separate components of spin current and spin accumulation.
void Save_Data(std::string path = filepath)
{
   // Creating and Opening Files
	std::ofstream Uxfile;
	std::ofstream Uyfile;
	std::ofstream Uzfile;
	std::ofstream Jmxfile;
	std::ofstream Jmyfile;
	std::ofstream Jmzfile;
	//Opening Files

	std::string Ux_filepath = path+"Ux";
	std::string Uy_filepath = path+"Uy";
	std::string Uz_filepath = path+"Uz";
	std::string Jmx_filepath = path+"Jmx";
	std::string Jmy_filepath = path+"Jmy";
	std::string Jmz_filepath = path+"Jmz";

	Uxfile.open(Ux_filepath.c_str());
	Uyfile.open(Uy_filepath.c_str());
	Uzfile.open(Uz_filepath.c_str());
	Jmxfile.open(Jmx_filepath.c_str());
	Jmyfile.open(Jmy_filepath.c_str());
	Jmzfile.open(Jmz_filepath.c_str());

   for (unsigned int i=0;i<xlen;i++)
   {
      if (i < xlen - 1)
      {
         Uxfile <<  i << ",";
         Uyfile <<  i << ",";
         Uzfile <<  i << ",";
         Jmxfile <<  i << ",";
         Jmyfile <<  i << ",";
         Jmzfile <<  i << ",";
      }
      else
      {
         Uxfile <<  i <<  "\n";
         Uyfile <<  i <<  "\n";
         Uzfile <<  i <<  "\n";
         Jmxfile <<  i <<  "\n";
         Jmyfile <<  i <<  "\n";
         Jmzfile <<  i <<  "\n";
      }
   }
   for (unsigned int n=0; n<Save_Ux.size(); n++)
   {
      for (unsigned int i=0;i<xlen;i++)
      {
         if (i < xlen - 1)
         {
            Uxfile <<  Save_Ux[n][i] << ",";
            Uyfile <<  Save_Uy[n][i] << ",";
            Uzfile <<  Save_Uz[n][i] << ",";
            Jmxfile <<  Save_Jmx[n][i] << ",";
            Jmyfile <<  Save_Jmy[n][i] << ",";
            Jmzfile <<  Save_Jmz[n][i] << ",";
         }
         else
         {
            Uxfile <<  Save_Ux[n][i] <<  "\n";
            Uyfile <<  Save_Uy[n][i] <<  "\n";
            Uzfile <<  Save_Uz[n][i] <<  "\n";
            Jmxfile <<  Save_Jmx[n][i] << "\n";
            Jmyfile <<  Save_Jmy[n][i] << "\n";
            Jmzfile <<  Save_Jmz[n][i] << "\n";
         }
      }
   }
   Uxfile.close();
   Uyfile.close();
   Uzfile.close();
   Jmxfile.close();
   Jmyfile.close();
   Jmzfile.close();
}


#endif
