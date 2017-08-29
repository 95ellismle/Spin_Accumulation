#ifndef TEST_H
#define TEST_H

#include <iostream>
#include <vector>
#include <algorithm>

extern const double minf;

template <typename T>
void print_vec(std::vector<T> v, char spacer='\n')
{
   for (unsigned int i=0; i<v.size(); i++)
   {
      std::cout << v[i] << spacer;
   }
}

bool is_valid(std::vector<double> v)
{
   double max = *std::max_element(v.begin(), v.end());
   if (max < minf*2)
   {
      for (unsigned int i=1; i<v.size(); i++)
      {
         //if (v[i] or v[i] == 0) std::cout << "bob" << '\n';
         //else std::cout << " " << '\n';
      }
      return true;
   }
   else
      return false;
}

#endif
