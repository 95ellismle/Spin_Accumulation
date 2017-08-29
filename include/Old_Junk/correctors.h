#ifndef CORRECT_H
#define  CORRECT_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

// Constants to use later
extern const double fact1;
extern const double fact2;
extern const double fact3[2];
extern const double minf_magnitude;						// The magnitude of the equilibrium spin accumlation vector 							              //
extern const double M1[3];							// F1 Magnetisation, this code will only work with an x and y magnetisation 		     			        //
extern const double M2[3];										// F2 Magnetisation, this code will only work with an x and y magnetisation	 		        //

// Vectors
extern double Byy;
extern double Bzz;
extern double B1yz;
extern double B2yz;
extern std::vector<double> vec_beta;
extern std::vector<double> diff_beta;
extern std::vector<double> Sx;
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


// These are functions to adjust the Ux vector to sort out the cross products and things.
double Z_NM_corrector(unsigned int i){
	return - fact3[1]*(Uz[i]);
}
double Z1_corrector(unsigned int i){
	return -B1yz*(Uy[i+1]-2*Uy[i]+Uy[i-1]) - fact3[0]*(Uz[i]-minf_magnitude*M1[2]) - fact1*(M1[1]*Ux[i]) - fact2*(-M1[1]*M1[2]*Uy[i] + M1[1]*M1[1]*Uz[i]);
}
double Z2_corrector(unsigned int i){
	return -B2yz*(Uy[i+1]-2*Uy[i]+Uy[i-1]) - fact3[0]*(Uz[i]-minf_magnitude*M2[2]) - fact1*(M2[1]*Ux[i]) - fact2*(-M2[1]*M2[2]*Uy[i] + M2[1]*M2[1]*Uz[i]);
}
double Z1_half_corrector(unsigned int i){
	return -B1yz*(Uy_half[i+1]-2*Uy_half[i]+Uy_half[i-1])  - fact3[0]*(Uz_half[i]-minf_magnitude*M1[2]) - fact1*(M1[1]*Ux[i]) - fact2*(-M1[1]*M1[2]*Uy_half[i] + M1[1]*M1[1]*Uz_half[i]);
}
double Z2_half_corrector(unsigned int i){
	return -B2yz*(Uy_half[i+1]-2*Uy_half[i]+Uy_half[i-1]) - fact3[0]*(Uz_half[i]-minf_magnitude*M2[2]) - fact1*(M2[1]*Ux[i]) - fact2*(-M2[1]*M2[2]*Uy_half[i] + M2[1]*M2[1]*Uz_half[i]);
}
double Y_NM_corrector(unsigned int i){
	return  -fact3[1]*(Uy[i]);
}
double Y1_corrector(unsigned int i){
   return -B1yz*(Uz[i+1]-2*Uz[i]+Uz[i-1]) - fact3[0]*(Uy[i]-minf_magnitude*M1[1]) + fact1*(M1[2]*Ux[i]) - fact2*(M1[2]*M1[2]*Uy[i] - M1[1]*M1[2]*Uz[i]);
}
double Y2_corrector(unsigned int i){
	return -B2yz*(Uz[i+1]-2*Uz[i]+Uz[i-1]) - fact3[0]*(Uy[i]-minf_magnitude*M2[1]) + fact1*(M2[2]*Ux[i]) - fact2*(M2[2]*M2[2]*Uy[i] - M2[1]*M2[2]*Uz[i]);
}
double Y1_half_corrector(unsigned int i){
   return -B1yz*(Uz_half[i+1]-2*Uz_half[i]+Uz_half[i-1])  - fact3[0]*(Uy_half[i]-minf_magnitude*M1[1]) + fact1*(M1[2]*Ux[i]) - fact2*(M1[2]*M1[2]*Uy_half[i] - M1[1]*M1[2]*Uz_half[i]);
}
double Y2_half_corrector(unsigned int i){
   return -B2yz*(Uz_half[i+1]-2*Uz_half[i]+Uz_half[i-1])  - fact3[0]*(Uy_half[i]-minf_magnitude*M2[1]) + fact1*(M2[2]*Ux[i]) - fact2*(M2[2]*M2[2]*Uy_half[i] - M2[1]*M2[2]*Uz_half[i]);
}
double X_NM_corrector(unsigned int i){
	return - fact3[1]*(Ux[i]);
}
double X1_corrector(unsigned int i){
	return - fact1*(M1[2]*Uy[i] - M1[1]*Uz[i]) - fact2*(Ux[i]*(M1[2]*M1[2] +M1[1]*M1[1])) - fact3[0]*(Ux[i]);
}
double X2_corrector(unsigned int i){
	return - fact1*(M2[2]*Uy[i] - M2[1]*Uz[i]) - fact2*(Ux[i]*(M2[2]*M2[2] +M2[1]*M2[1])) - fact3[0]*(Ux[i]);
}
double X1_half_corrector(unsigned int i){
	return - fact1*(M1[2]*Uy_half[i] - M1[1]*Uz_half[i]) - fact2*(Ux[i]*(M1[2]*M1[2] +M1[1]*M1[1])) - fact3[0]*(Ux[i]);
}
double X2_half_corrector(unsigned int i){
	return - fact1*(M2[2]*Uy_half[i] - M2[1]*Uz_half[i]) - fact2*(Ux[i]*(M2[2]*M2[2] +M2[1]*M2[1])) - fact3[0]*(Ux[i]);
}

#endif
