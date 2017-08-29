#ifndef FUNCS_H
#define FUNCS_H

#include <vector>
#include <cmath>


extern const double dx;

// The A Matrix
extern std::vector<double> A;


// A function to check if 2 floats are equivalent
// a and b are the 2 floats to check for equivalency.
bool almost_equal(double a, double b, double tolerance)
{
   if ((a-b)*(a-b)<= tolerance*tolerance)
   {
      return true;
   }
   else
   {
      return false;
   }
}


// A function to differentiate a vector.
std::vector<double> Differentiator(std::vector<double> input_array)
{
	for(unsigned int i=0;i<input_array.size()-1;i++)
	{
		input_array[i] = (input_array[i+1]-input_array[i])/dx;
	}
	input_array[input_array.size()-1] = input_array[input_array.size()-2] - (input_array[input_array.size()-2]-input_array[input_array.size()-3])/dx;
	return input_array;
}





#endif
