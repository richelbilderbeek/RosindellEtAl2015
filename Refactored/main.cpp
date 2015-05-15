// Code for the article
//
//   Rosindell, James, Luke J. Harmon, and Rampal S. Etienne.
//   "Unifying ecology and macroevolution with individual‐based theory."
//   Ecology letters 18.5 (2015): 472-482.
//
//
// Original version:
//   Version 1.2 Jan 2014
//   Author: James Rosindell
//
// Refactoring by Richel Bilderbeek


// NNIM = nearly neutral individual mutations

// positive or negative mutations at individual level,
// speciation is the result of a number of mutations together in build up

// parameters:

// JM (simulations size of metacommunity - non spatial)
// u probability of mutation ±
// s size of mutation
// n number of mutations necessary to count as a new species (accounted for here for fitness variations as a spectrum)

#include <cassert>
#include "ntsim.h"


int main(int argc, char* argv[])
{
  const bool debug{false};
  if (debug)
  {
    NTsim test;
    std::cout
      << "simulating M with fitness increase or decrease ... "
      << std::endl
    ;
    test.sim_all(2,10000,0.05,0.001,60);
    test.write_files();
  }

  // first get the command line arg data in
  const long task_iter = argc > 1 ? std::stoi(argv[1]) - 1 : 3;
  double task_max = argc > 2 ? std::stof(argv[2]) : 1; // the max time (minutes)


  std::cout << task_iter << " = task \n";
  std::cout << task_max << " = max time \n";

  std::vector<long> JM_vect;
  std::vector<double> nu_vect;
  std::vector<double> S_vect;

  JM_vect.clear();
  nu_vect.clear();
  S_vect.clear();

  JM_vect.push_back(100);
  JM_vect.push_back(1000);

  /*
  JM_vect.push_back(10000);
  JM_vect.push_back(100000);
  JM_vect.push_back(1000000);
  */
  nu_vect.push_back(0.0001);
  nu_vect.push_back(0.0003);
  nu_vect.push_back(0.001);
  nu_vect.push_back(0.003);
  nu_vect.push_back(0.01);
  nu_vect.push_back(0.03);
  nu_vect.push_back(0.1);
  nu_vect.push_back(0.3);

  S_vect.push_back(0.01);
  S_vect.push_back(0.001);
  S_vect.push_back(0.0001);
  S_vect.push_back(0.00001);

  std::vector<long> JM_list;
  std::vector<double> nu_list;
  std::vector<double> S_list;

  for (unsigned int i = 0 ; i < JM_vect.size() ; ++i)
  {
    for (unsigned int j = 0 ; j < nu_vect.size() ; j ++)
    {
      for (unsigned int k = 0 ; k < S_vect.size() ; k ++)
      {
        JM_list.push_back(JM_vect[i]);
        nu_list.push_back(nu_vect[j]);
        S_list.push_back(S_vect[k]);
      }
    }
  }

  NTsim sim;

  if (argc == 4)
  {
    // restore
    sim.sim_restore(argv[3],60*task_max);
  }
  else
  {
    std::cout << "simulating M with fitness increase or decrease ... ";
    std::cout << JM_list[task_iter] << " , " << nu_list[task_iter] << " , " <<  S_list[task_iter] << " \n";

    sim.sim_all(task_iter+10000 , JM_list[task_iter] , nu_list[task_iter] , S_list[task_iter] , 60*task_max);
  }
  sim.write_files();
}

// notes: copy files to $WORK where there is 150 gb of storage available
// restrict output to 1M lines (0.22 GB storage and 1.1 hrs per sorting) - done in main simulation routine loop
// 0.22 Mb per 1K lines and 3.95 seconds sorting per 1K lines

