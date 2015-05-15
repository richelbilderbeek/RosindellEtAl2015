#ifndef NTSIM_H
#define NTSIM_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include "nrrand.h"
#include "species.h"

/************************************************************
 NEUTRAL SIMULATION OBJECT
 ************************************************************/

struct NTsim {
  NTsim();

  void set_seed(long seedin);
  void setup(long seedin , long J_M_in , double mu_in , double s_in);

  // this function simulates a set of parameters for a set time and gets as many results as it can
  // all these reults are dropped in a common file named based on the seed and parameter compbination
  // the time is measured in minutes
  void sim_all(long seedin , long J_M_in , double mu_in , double s_in , double simtime);

  void sim_restore(char* fnamein,double simtime);
  void sim_step();

  void total_clear();

  // this function outputs the data including calculating the fitness variation data
  // it also sorts it so that no longer relevant parts are thrown out.
  void output();

  void write_files();

  #ifndef NDEBUG
  static void Test() noexcept;
  #endif

private:

  // *** simulation parameters ***

  long J_M;
  // metacommunity size
  double mu;
  // mutation rate (not the speciation rate any more) %^&
  double s;
  // size of selective advantage
  // note: that the parameter n that decides how many mutations apart a species is does not have to be defined here
  // it can be accounted for later from some outputs.  Also here we can calculate a spectrum of different values of n up to the value where everything becomes one species.

  // *** suspend / resume variables ***

  long num_suspends;
  // the number of times this simulation has been suspended before
  long num_file_outputs;
  // the number of times this simulation has yielded file outputs

  // *** general simulation data ***

  double generation;
  // current generation relative to others
  double true_generation;
  // the true number of generations since the simulation started
  double gen_burned_in;
  // the generation at which the system was approximately burned in - all species have a common ancestor in the simulation.
  long simnum;
  // the number of readings taken from this simulation so far - bearing in mind they are not quite independent the have gen_burned_in generations between them
  long max_selection_level;
  // what is the highest selection level in the system (larger means fitter - these are integers)
  long richness;
  // the current (sub) species richness in the system

  // *** recycling utility for (sub) species numbers ***

  long upto;
  // what number are we up to with species labeling
  std::vector<long> recycler;
  // a list of indices in the per-species data that can now be used because the species are extinct
  long end_recycler;
  // the end of the recycler std::vector (in case it is required to contract) so it's equal to 0 when recyler is empty

  // *** burn in checker and per individual data - becomes irrelevant after first output ***

  std::vector<long> init_ancestor;
  // the ancestor of this individual in the initial state - index referres to currently living individuals in the metacommunity
  std::vector<long> num_descend;
  //the number of current descendents of this individual - index referrs to the index of individuals initialised in the metacommunity
  long num_remain;
  // the number of values in num_descend that are non zero (when this decreases to one we're burned in)

  // *** per individual data that remains relevant after burn in

  std::vector<long> metacommunity;
  // metacommunity - contains the index in per-species data of the species of each individual (for ease of simulation) - the index of this std::vector is the position in the metacommunity

  // *** per-species simulation data ***

  std::vector<long> abundances;
  // the abundnance of each species
  std::vector<long> parents;
  // the parent species for each species
  std::vector<double> born;
  // the date of species birth for each species
  std::vector<double> died;
  // the date of species death for each species (-1 means not dead)
  std::vector<long> selection_level;
  // the selection level of each species (larger numbers mean greater advantage)
  std::vector<long> num_mut;
  // number of mutations along the branch that connects this to the parent node, these could be positive or negative mutations

  std::vector<long> max_abundance;
  // the maximum abundance ever experienced by this species
  std::vector<double> maxabtime1;
  // the earliest date when that maximum abundance was hit
  std::vector<double> maxabtime2;
  // the latest date when that maximum abundance was hit

  std::vector<long> num_descend_spec;
  // the number of species (lineages) that descend from this lineage and are still alive (including the species itself if it's still alive, not if it's dead)

  long iter; // keeps count of the output lines
  long iter2; // keeps count of number output lines from this current simulation only

  // *** random number generator ***

  NRrand NR;
  // generator object
  bool seeded;
  // is it seeded
  long the_seed;
  // the seed

  // *** data that will be outputted to files *** //

  // this data is stored in RAM and then outputted to files later.
  // there will be four files: newick, out, res, vars and fitvar

  // newick
  std::vector<std::string> FILE_final_newick;

  // out (per species)

  // header of the filal file will read
  // generation,iname,iname2,index,abundances,parents,born,died,selection_level,num_mut,max_abundance,maxabtime1,maxabtime2,num_descend_spec,divergence_date,relgen,tage,page

  std::vector<double> FILE_generation;
  std::vector<std::string> FILE_iname;
  std::vector<long> FILE_iname2;
  std::vector<long> FILE_index;
  std::vector<long> FILE_abundances;
  std::vector<long> FILE_parents;
  std::vector<double> FILE_born;
  std::vector<double> FILE_died;
  std::vector<long> FILE_selection_level;
  std::vector<long> FILE_num_mut;
  std::vector<long> FILE_max_abundance;
  std::vector<double> FILE_maxabtime1;
  std::vector<double> FILE_maxabtime2;
  std::vector<long> FILE_num_descend_spec;
  std::vector<double> FILE_divergence_date;
  std::vector<double> FILE_relgen;
  std::vector<double> FILE_tage;
  std::vector<double> FILE_page;

  // res and vars are based on overall sim parameters at end of sim and don't need to be saved.

  // fitvar shows the fitness variation overall and within each species for a variety of differnet values of n
  // header of the final file will read: generation_var, meanc_all , varc_all , n , num_subspec , meanc , varc , maxc , minc

  std::vector<double> FILE_generation_var; // will line up with the variance output rows for frequency of output here
  std::vector<long> FILE_n; // value of n that applies for the rest of the readings (the first three above will probably repeat a lot)
  std::vector<long> FILE_n_richness; // the species richness overall at this value of n
  std::vector<double> FILE_num_subspec; // number of sub species - average over all species
  std::vector<double> FILE_tot_abund; // total abundance of all sub species - average over all species
  std::vector<double> FILE_meanc; // mean value of c - (within species) average over all species
  std::vector<double> FILE_varc; // variange in c - (within species) average over all species
  std::vector<double> FILE_rangec; // range in value of c - (within species) average over all species
  std::vector<double> FILE_meanc_w; // mean value of c - (within species) weighted average over all species
  std::vector<double> FILE_varc_w; // variange in c - (within species) weighted average over all species
  std::vector<double> FILE_rangec_w; // range in value of c - (within species) weighted average over all species

};

#endif // NTSIM_H
