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

// A function to linearly
double lin_increase_je(double t, double teq, double je)
{
   if (t < teq) return t * (je/teq);
   else return je;
}
//
// Equivalent to numpy.arange() in python
template<typename T>
std::vector<T> range(T start, T stop, T step = 1) {
    std::vector<T> values;
    for (T value = start; value < stop; value += step)
        values.push_back(value);
    return values;
}

// A function to differentiate a vector 1st order centered difference.
std::vector<double> gradient(std::vector<double> input_array, double h = dx)
{
   std::vector<double> grad_array (input_array.size());
	for(unsigned int i=1;i<input_array.size()-1;i++)
	{
		double m1 = (input_array[i+1]-input_array[i-1])/(2*h);
		//double m2 = (input_array[i+2]-input_array[i-2])/(4*h);
		//double m3 = (input_array[i+3]-input_array[i-3])/(6*h);
      grad_array[i] = m1;//1.5*m1 - 0.6*m2 + 0.1*m3;
	}
   //grad_array[2] = (1.5/(2*h))*(grad_array[1] - grad_array[3]) - (0.6/(3*h))*(grad_array[0] - grad_array[4]);
   //grad_array[1] = (grad_array[0] - grad_array[2])/(2*h);
   grad_array[0] = grad_array[1] - (grad_array[2]-grad_array[1]);
	return grad_array;
}





#endif
