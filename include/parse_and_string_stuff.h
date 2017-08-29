#ifndef TYPE_H
#define TYPE_H

#include <string>
#include <sstream>

// This function will convert an int variable to a string
// It takes an int variable and outputs a string
std::string Int2Str(int num)
{
   std::stringstream sstr_i;
   sstr_i << num;
   std::string str_i = sstr_i.str();
   return str_i;
}

// Converts a double to a string
// takes a double and outputs a string
std::string Doub2Str(double num)
{
   std::stringstream sstr_i;
   sstr_i << num;
   std::string str_i = sstr_i.str();
   return str_i;
}

// Will parse integer command line arguments
int cmd_parse(std::string buff, int argc, char* argv[])
{
   int var=0;
   for(int arg_num = 0; arg_num<argc; arg_num++)
   {
      std::string arg = argv[arg_num];
      int start_ind = arg.find(buff);
      if (start_ind != -1)
      {
         var = std::stoi(arg.substr(start_ind+2, arg.size()));
      }
   }
   return var;
}
#endif
