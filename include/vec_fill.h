#ifndef VEC_FILL_H
#define VEC_FILL_H

#include <vector>
#include <cmath>

#include "math_funcs.h"

extern double dt;

extern const double je_peak;
extern const double charge_current_teq;
extern const double DL;
extern const double boundaries[4];
extern const unsigned int xlen;



std::vector<double> diffuse_calc(std::vector<double> x, double non_mag, double mag1, double mag2)
{
    mag1 = mag1 - non_mag;
    mag2 = mag2 - non_mag;

    double mid_points[3] = {boundaries[1] - boundaries[0], boundaries[2] - boundaries[1], boundaries[3] - boundaries[2]};
    double y10, y1;
    double y20, y2;
    double y30, y3;
    double y40, y4;
    y10 = mag1/2*( std::tanh( ( x[0] - boundaries[0] ) / (DL/2) ) + 1 );
    y20 = -mag1/2*( std::tanh( ( x[0] - boundaries[0] - mid_points[0] ) / (DL/2) ) + 1 ) + mag1;
    y30 = mag2/2*( std::tanh( ( x[0] - boundaries[1] - mid_points[1] ) / (DL/2) ) + 1 );
    y40 = -mag2/2*( std::tanh( ( x[0] - boundaries[2] - mid_points[2] ) / (DL/2) ) + 1 ) + mag2;
    std::vector<double> diff (xlen);
    for (unsigned int i=0; i<xlen; i++)
    {
       y1 = mag1/2*( std::tanh( ( x[i] - boundaries[0] ) / (DL/2) ) + 1 );
       y2 = -mag1/2*( std::tanh( ( x[i] - boundaries[0] - mid_points[0] ) / (DL/2) ) + 1 ) + mag1;
       y3 = mag2/2*( std::tanh( ( x[i] - boundaries[1] - mid_points[1] ) / (DL/2) ) + 1 );
       y4 = -mag2/2*( std::tanh( ( x[i] - boundaries[2] - mid_points[2] ) / (DL/2) ) + 1 ) + mag2;
       diff[i] = (y1 + y2 + y3 + y4) - (y10+y20+y30+y40)+non_mag;
   }
    return diff;
}

#endif
