#ifndef TEST_H
#define TEST_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>

extern const double minf;

template <typename T>
void print_vec(std::vector<T> v, char spacer='\n')
{
   for (unsigned int i=0; i<v.size(); i++)
   {
      std::cout << v[i] << spacer;
   }
}

// Checks if the data is within a valid range and (checks if there are no null values... -need to code up)
bool is_valid(std::vector<double> v)
{
   double max = *std::max_element(v.begin(), v.end());
   double min = *std::min_element(v.begin(), v.end());
   if (max < minf*2 and min > -minf*2)
   {
      for (unsigned int i=0; i<v.size(); i++)
      {
         if ( std::isnan(v[i]) ){
            return false;
         }
      }
      return true;
   }
   else
      return false;
}

// Tests whether a folder exists to save the data and if not will create one. (need to code up)
bool folder_exsists()
{
   return true;
}

#endif
