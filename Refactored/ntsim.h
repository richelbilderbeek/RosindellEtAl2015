// Code for the article
//
//   Rosindell, James, Luke J. Harmon, and Rampal S. Etienne.
//   "Unifying ecology and macroevolution with individual‐based theory."
//   Ecology letters 18.5 (2015): 472-482.
//
// Original version:
//   Version 1.2 Jan 2014
//   Author: James Rosindell
//
// Refactoring by Richel Bilderbeek

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

  const auto& GetFitnessCategories() const noexcept { return m_fitness_categories; }


  void set_seed(long seedin);
  void setup(
    const long seedin,
    const long J_M_in,
    const double mu_in,
    const double s_in
  );

  // this function simulates a set of parameters for a set time and gets as many results as it can
  // all these reults are dropped in a common file named based on the seed and parameter compbination
  // the time is measured in minutes
  void sim_all(
    const long seedin,
    const long J_M_in,
    const double mu_in,
    const double s_in,
    const double simtime
  );

  void sim_restore(const char* const fnamein, const double simtime);
  void sim_step();

  void total_clear();

  // this function outputs the data including calculating the fitness variation data
  // it also sorts it so that no longer relevant parts are thrown out.
  void output() const;

  void write_files() const;

  #ifndef NDEBUG
  static void Test() noexcept;
  #endif

private:

  ///The size of the metacommunity size
  ///Symbol used in article: J_M
  long m_metacommunity_size;

  ///Mutation rate (not the speciation rate any more)
  ///Symbol used in article: mu
  double m_mutation_rate;

  ///Selection strength. Also: size of selective advantage
  ///Symbol used in article: s
  // note: that the parameter n that decides how many mutations apart a species is does not have to be defined here
  // it can be accounted for later from some outputs.
  // Also here we can calculate a spectrum of different values of n up to the value where everything becomes one species.
  double m_selection_strength;


  // *** general simulation data ***

  mutable double m_current_generation;
  // current generation relative to others
  double m_true_generation;
  // the true number of generations since the simulation started
  double m_gen_burned_in;
  // the generation at which the system was approximately burned in - all species have a common ancestor in the simulation.
  long m_simnum;
  // the number of readings taken from this simulation so far - bearing in mind they are not quite independent the have gen_burned_in generations between them
  long m_max_selection_level;
  // what is the highest selection level in the system (larger means fitter - these are integers)
  long m_current_richness;
  // the current (sub) species richness in the system

  // *** recycling utility for (sub) species numbers ***

  long m_upto_which_label;
  // what number are we up to with species labeling
  mutable std::vector<long> m_recycler;
  // a list of indices in the per-species data that can now be used because the species are extinct
  mutable long m_end_recycler;
  // the end of the recycler std::vector (in case it is required to contract) so it's equal to 0 when recyler is empty

  // *** burn in checker and per individual data - becomes irrelevant after first output ***

  std::vector<long> m_init_ancestor;
  // the ancestor of this individual in the initial state - index referres to currently living individuals in the metacommunity
  std::vector<long> m_n_descend;
  //the number of current descendents of this individual - index referrs to the index of individuals initialised in the metacommunity
  long m_n_remain;
  // the number of values in num_descend that are non zero (when this decreases to one we're burned in)

  // *** per individual data that remains relevant after burn in

  std::vector<long> m_metacommunity;
  // metacommunity - contains the index in per-species data of the species of each individual (for ease of simulation) - the index of this std::vector is the position in the metacommunity

  // *** per-species simulation data ***

  mutable std::vector<long> m_abundances;
  // the abundnance of each species
  mutable std::vector<long> m_parents;
  // the parent species for each species
  mutable std::vector<double> m_born;
  // the date of species birth for each species
  mutable std::vector<double> m_died;
  // the date of species death for each species (-1 means not dead)

  ///The fitness category of each species (also: the selection level of each species)
  ///A higher fitness category means being a better competitor
  mutable std::vector<long> m_fitness_categories;

  mutable std::vector<long> m_n_mut;
  // number of mutations along the branch that connects this to the parent node, these could be positive or negative mutations

  mutable std::vector<long> m_max_abundance;
  // the maximum abundance ever experienced by this species
  mutable std::vector<double> m_maxabtime1;
  // the earliest date when that maximum abundance was hit
  mutable std::vector<double> m_maxabtime2;
  // the latest date when that maximum abundance was hit

  mutable std::vector<long> m_n_descend_spec;
  // the number of species (lineages) that descend from this lineage and are still alive (including the species itself if it's still alive, not if it's dead)

  mutable long m_iter; // keeps count of the output lines
  mutable long m_iter2; // keeps count of number output lines from this current simulation only

  // *** random number generator ***

  mutable NRrand m_rng;
  // generator object
  bool m_is_seeded;
  // is it seeded
  long m_seed_used;
  // the seed

  // *** suspend / resume variables ***
  /// the number of times this simulation has been suspended before
  mutable long m_n_suspends;
  /// the number of times this simulation has yielded file outputs
  mutable long m_n_file_outputs;

  // *** data that will be outputted to files *** //

  // this data is stored in RAM and then outputted to files later.
  // there will be four files: newick, out, res, vars and fitvar

  // newick
  mutable std::vector<std::string> m_FILE_final_newick;

  // out (per species)

  // header of the filal file will read
  // generation,iname,iname2,index,abundances,parents,born,died,selection_level,num_mut,max_abundance,maxabtime1,maxabtime2,num_descend_spec,divergence_date,relgen,tage,page

  mutable std::vector<double> m_FILE_generation;
  mutable std::vector<std::string> m_FILE_iname;
  mutable std::vector<long> m_FILE_iname2;
  mutable std::vector<long> m_FILE_index;
  mutable std::vector<long> m_FILE_abundances;
  mutable std::vector<long> m_FILE_parents;
  mutable std::vector<double> m_FILE_born;
  mutable std::vector<double> m_FILE_died;
  mutable std::vector<long> m_FILE_selection_level;
  mutable std::vector<long> m_FILE_num_mut;
  mutable std::vector<long> m_FILE_max_abundance;
  mutable std::vector<double> m_FILE_maxabtime1;
  mutable std::vector<double> m_FILE_maxabtime2;
  mutable std::vector<long> m_FILE_num_descend_spec;
  mutable std::vector<double> m_FILE_divergence_date;
  mutable std::vector<double> m_FILE_relgen;
  mutable std::vector<double> m_FILE_tage;
  mutable std::vector<double> m_FILE_page;

  // res and vars are based on overall sim parameters at end of sim and don't need to be saved.

  // fitvar shows the fitness variation overall and within each species for a variety of differnet values of n
  // header of the final file will read: generation_var, meanc_all , varc_all , n , num_subspec , meanc , varc , maxc , minc

  mutable std::vector<double> m_FILE_generation_var; // will line up with the variance output rows for frequency of output here
  mutable std::vector<long> m_FILE_n; // value of n that applies for the rest of the readings (the first three above will probably repeat a lot)
  mutable std::vector<long> m_FILE_n_richness; // the species richness overall at this value of n
  mutable std::vector<double> m_FILE_num_subspec; // number of sub species - average over all species
  mutable std::vector<double> m_FILE_tot_abund; // total abundance of all sub species - average over all species
  mutable std::vector<double> m_FILE_meanc; // mean value of c - (within species) average over all species
  mutable std::vector<double> m_FILE_varc; // variange in c - (within species) average over all species
  mutable std::vector<double> m_FILE_rangec; // range in value of c - (within species) average over all species
  mutable std::vector<double> m_FILE_meanc_w; // mean value of c - (within species) weighted average over all species
  mutable std::vector<double> m_FILE_varc_w; // variange in c - (within species) weighted average over all species
  mutable std::vector<double> m_FILE_rangec_w; // range in value of c - (within species) weighted average over all species

};

#endif // NTSIM_H
