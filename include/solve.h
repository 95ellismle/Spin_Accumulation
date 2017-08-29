#ifndef SOLVE_H
#define SOLVE_H

#include <vector>

#include "math_funcs.h"
#include "vec_fill.h"

extern double dt;
extern const double minf;

extern const unsigned int xlen;

extern std::vector<double> fact1;
extern std::vector<double> fact2;
extern std::vector<double> fact3;

extern std::vector<double> Mx;
extern std::vector<double> My;
extern std::vector<double> Mz;

extern std::vector<double> Ux;
extern std::vector<double> Ux_new;
extern std::vector<double> Uy;
extern std::vector<double> Uy_new;
extern std::vector<double> Uz;
extern std::vector<double> Uz_new;

extern std::vector<double> jmx;
extern std::vector<double> jmy;
extern std::vector<double> jmz;

// Fills a vector with the values of spin current
std::vector<double> spin_current(std::vector<double> m, std::vector<double> M)
{
   std::vector<double> jm (xlen);
   std::vector<double> diff_m = gradient(m);
   for (unsigned int i = 0; i < xlen; i++) {
      double first_term = beta[i] * M[i] * je_peak;//lin_increase_je(i*dt, charge_current_teq, je_peak);
      double second_term = diff_m[i];
      double third_term = beta[i] * beta_prime[i] * M[i] * M[i] * diff_m[i];
      //std::cout << diff_m[i] << '\n';
      jm[i] = first_term - 2 * D[i] * (second_term - third_term);
   }
   return jm;
}

// Carries out the forward Euler method on Ux
std::vector<double> Ux_solve(std::vector<double> Ux, std::vector<double> Ux_new)
{
   auto grad_jmx = gradient(jmx);
   for( unsigned int i=0; i<xlen; i++)
   {
      double first_term = Ux[i];
      double second_term = dt*grad_jmx[i];
      double third_term = fact1[i] * (-My[i]*Uz[i] + Mz[i]*Uy[i] );
      double fourth_term = fact2[i] * (-Mx[i]*My[i]*Uy[i] - Mx[i]*Mz[i]*Uz[i] + Ux[i]*(My[i]*My[i] + Mz[i]*Mz[i]) );
      double fifth_term = fact3[i] * (Ux[i] - minf*Mx[i]);
      Ux_new[i] = first_term - (second_term + third_term + fourth_term + fifth_term);
   }
   return Ux_new;
}

// Carries out the forward Euler method on Uy
std::vector<double> Uy_solve(std::vector<double> Uy, std::vector<double> Uy_new)
{
   auto grad_jmy = gradient(jmy);
   for( unsigned int i=0; i<xlen; i++)
   {
      double first_term = Uy[i];
      double second_term = dt*grad_jmy[i];
      double third_term = fact1[i] * ( Mx[i]*Uz[i] - Mz[i]*Ux[i] );
      double fourth_term = fact2[i] * (Mx[i]*Mx[i]*Uy[i] - Mx[i]*My[i]*Ux[i] - My[i]*Mz[i]*Uz[i] + Mz[i]*Mz[i]*Uy[i]);
      double fifth_term = fact3[i] * (Uy[i] - minf*My[i]);
      Uy_new[i] = first_term - (second_term + third_term + fourth_term + fifth_term);
   }
   return Uy_new;
}

// Carries out the forward Euler method on Uz
std::vector<double> Uz_solve(std::vector<double> Uz, std::vector<double> Uz_new)
{
   auto grad_jmz = gradient(jmz);
   for( unsigned int i=0; i<xlen; i++)
   {
      double first_term = Uz[i];
      double second_term = dt*grad_jmz[i];
      double third_term = fact1[i] * ( -Mx[i]*Uy[i] + My[i]*Ux[i] );
      double fourth_term = fact2[i] * (Mx[i]*Mx[i]*Uz[i] - Mx[i]*Mz[i]*Ux[i] + My[i]*My[i]*Uz[i] - My[i]*Mz[i]*Uy[i] );
      double fifth_term = fact3[i] * (Uz[i] - minf*Mz[i]);
      Uz_new[i] = first_term - (second_term + third_term + fourth_term + fifth_term);
   }
   return Uz_new;
}

//Wraps all the solving steps into 1 function, will hopefully make the main code more readable!
void Solver(std::vector<double> &Ux, std::vector<double> &Uy, std::vector<double> &Uz, std::vector<double> &jmx, std::vector<double> &jmy, std::vector<double> &jmz )
{
   // Solving for Ux
   Ux_new = Ux_solve(Ux, Ux_new);
   // Solving for Uy
   Uy_new = Uy_solve(Uy, Uy_new);
   // Solving for Uy
   Uz_new = Uz_solve(Uz, Uz_new);

   jmx = spin_current(Ux_new, Mx);
   jmy = spin_current(Uy_new, My);
   jmz = spin_current(Uz_new, Mz);
   Ux = Ux_new;
   Uy = Uy_new;
   Uz = Uz_new;
}

#endif
