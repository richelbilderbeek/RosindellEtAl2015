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

#include "ntsim.h"

#include <cassert>

#include "helper.h"

NTsim::NTsim()
  : m_metacommunity_size{0},
    m_mutation_rate{0.0},
    m_selection_strength{0.0},
    m_current_generation{0.0},
    m_true_generation{0.0},
    m_gen_burned_in{0.0},
    m_simnum{0},
    m_max_selection_level{0},
    m_current_richness{0},
    m_upto_which_label{0},
    m_recycler{},
    m_end_recycler{0},
    m_init_ancestor{},
    m_n_descend{},
    m_n_remain{0},
    m_metacommunity{},
    m_abundances{},
    m_parents{},
    m_born{},
    m_died{},
    m_fitness_categories{},
    m_n_mut{},
    m_max_abundance{},
    m_maxabtime1{},
    m_maxabtime2{},
    m_n_descend_spec{},
    m_iter{0},
    m_iter2{0},
    m_rng{},
    m_is_seeded{false},
    m_seed_used{0},
    m_n_suspends{0},
    m_n_file_outputs{0},
    m_FILE_final_newick{},
    m_FILE_generation{},
    m_FILE_iname{},
    m_FILE_iname2{},
    m_FILE_index{},
    m_FILE_abundances{},
    m_FILE_parents{},
    m_FILE_born{},
    m_FILE_died{},
    m_FILE_selection_level{},
    m_FILE_num_mut{},
    m_FILE_max_abundance{},
    m_FILE_maxabtime1{},
    m_FILE_maxabtime2{},
    m_FILE_num_descend_spec{},
    m_FILE_divergence_date{},
    m_FILE_relgen{},
    m_FILE_tage{},
    m_FILE_page{},
    m_FILE_generation_var{},
    m_FILE_n{},
    m_FILE_n_richness{},
    m_FILE_num_subspec{},
    m_FILE_tot_abund{},
    m_FILE_meanc{},
    m_FILE_varc{},
    m_FILE_rangec{},
    m_FILE_meanc_w{},
    m_FILE_varc_w{},
    m_FILE_rangec_w{}
{
  #ifndef NDEBUG
  Test();
  #endif

  total_clear();
}



void NTsim::output() const
{
  // std::cout << "outputting \n";

  const bool output_phy = true; // output the phylogeny?
  const bool output_var = true; // output the variance?

  // first of all work out if each line should appear in the final output
  // extinct species that were not a most recent common ancestor between two species are not outputted
  // the date of divergence from the (grand)parent species is noted in another column of data now
  std::vector<bool> keep_data;
  std::vector<double> divergence_date;
  std::vector<long> intermediate_parents;
  std::vector<long> num_mut2 = m_n_mut;
  keep_data.clear();
  divergence_date.clear();
  intermediate_parents.clear();
  for (unsigned int i = 0 ; i < m_parents.size() ; ++i)
  {
    keep_data.push_back(true);
    divergence_date.push_back(-1.0);
    intermediate_parents.push_back(m_parents[i]);
  }
  for (unsigned int i = 0 ; i < m_parents.size() ; ++i)
  {
    if (m_abundances[i] < 0)
    {
      keep_data[i] = false;
    }
    if (keep_data[i])
    {
      long current_index = i;
      while(((current_index >= 0)&&(m_parents[current_index] >= 0))&&(m_n_descend_spec[m_parents[current_index]]==m_n_descend_spec[i]))
      {
        current_index = m_parents[current_index];
        if (current_index >= 0)
        {
          keep_data[current_index] = false;
          m_n_mut[i] += num_mut2[current_index]; // accumulate the mutations of all the intermediate link species that only have one descendent species %^&
        }
      }
      if (m_parents[current_index] >= 0)
      {
        divergence_date[i] = m_born[current_index];
        intermediate_parents[i] = m_parents[current_index];
      }
      else
      {
        divergence_date[i] = -1;
        intermediate_parents[i] = -1;
      }
    }
  }

  // now relabel all the indices

  std::vector<long> newi;
  newi.clear();

  for (unsigned int i = 0 ; i < m_parents.size() ; ++i)
  {
    if (keep_data[i])
    {
      newi.push_back(m_iter);
      m_iter ++;
    }
    else
    {
      newi.push_back(-1);
    }
  }

  m_iter2 = 1; // note - I added iter 2 above to count the number of lines in this simulation on its own
  std::vector<long> newi2; // counts from 1 up in each simulation on its own
  newi2.clear();

  for (unsigned int i = 0 ; i < m_parents.size() ; ++i)
  {
    if (keep_data[i])
    {
      newi2.push_back(m_iter2);
      m_iter2 ++;
    }
    else
    {
      newi2.push_back(-1);
    }
  }

  // first sort out the data to remove lineages that are extinct and therfore free up memory.
  // NOTE - THE OUTPUT CODE COULD BE MOVED HERE IF WE WANT TO OUTPUT A FULL TREE INCLUDING ALL EXTINCT LINEAGES

  // now sort out the newick output
  // each node stores a num descendents that are accounted for so far and two std::vectors
  // the std::vectors are of child nodes and store their newick std::strings and dates

  std::string final_newick; // the final answer
  std::vector<long> num_descend_done;
  std::vector< std::vector<std::string> > newick_parts;
  std::vector< std::vector<double> > join_dates;
  std::vector< std::vector<double> > last_split_dates; // the last spit of a lineage - needed so that joindates-lastsplit generally gives us a branch length

  std::vector<double> page; // the phylogenetic age of the species for outputting.

  num_descend_done.clear();
  newick_parts.clear();
  join_dates.clear();
  last_split_dates.clear();

  page.clear();


  if (output_phy)
  {

    for (unsigned int i = 0 ; i < m_parents.size() ; ++i)
    {
      std::string temps = "";
      std::vector<std::string> tempvs;
      tempvs.clear();
      std::vector<double> tempvd;
      tempvd.clear();
      if (keep_data[i])
      {
        num_descend_done.push_back(0); // negative number indicates completed, positive number indicates number of descendents added to the newick_parts and join_dates std::vectors
      }
      else
      {
        num_descend_done.push_back(-2);
      }
      newick_parts.push_back(tempvs);
      join_dates.push_back(tempvd);
      last_split_dates.push_back(tempvd);

      page.push_back(-2);
    }

    for (unsigned int i = 0 ; i < m_parents.size() ; ++i)
    {
      if (m_abundances[i] > 0) // automatically excludes any lineages that might not be kept
      {
        // we're at the end
        std::ostringstream ss;
        ss.precision(10);
        ss << alpha_format_auto(newi[i]);
        ss << "["; // %^&
        ss << newi2[i]; // %^&
        ss << "]"; // %^&
        newick_parts[i].push_back(ss.str());
        last_split_dates[i].push_back(m_current_generation);
        join_dates[i].push_back(m_current_generation);
        num_descend_done[i] += 1;
      }
    }

    bool keep_looping = true;
    while(keep_looping)
    {

      keep_looping = false;
      for (unsigned long i = 0 ; i < m_parents.size() ; ++i)
      {
        // check if we're done
        if ((num_descend_done[i] >= 0)&&(keep_data[i])) // num_descend done goes negative when we're done
        {
          keep_looping = true;

          if (num_descend_done[i] == m_n_descend_spec[i])
          {

            // can sort out, otherwise can do nothing
            // first order the two std::vectors
            if ((newick_parts[i]).size() > 1) // only if ordering needed
            {
              bool unsorted = true;
              while(unsorted)
              {
                unsorted = false;
                for (unsigned int j = 1 ; j < (newick_parts[i]).size() ; j ++)
                {
                  if (join_dates[i][j-1] > join_dates[i][j])
                  {
                    unsorted = true;

                    double tempd = join_dates[i][j-1];
                    join_dates[i][j-1] = join_dates[i][j];
                    join_dates[i][j] = tempd;

                    tempd = last_split_dates[i][j-1];
                    last_split_dates[i][j-1] = last_split_dates[i][j];
                    last_split_dates[i][j] = tempd;

                    std::string temps = newick_parts[i][j-1];
                    newick_parts[i][j-1] = newick_parts[i][j];
                    newick_parts[i][j] = temps;
                  }
                }
              }
            }

            // this is the full std::string
            std::ostringstream ss2;
            ss2.precision(10);
            double thislastsplit = -1.0;

            if ((newick_parts[i]).size() > 1) // only if ordering needed
            {
              // this is the interior node name
              std::string int_node;
              std::ostringstream ss;
              ss.precision(10);
              ss << alpha_format_auto(newi[i]);
              ss << "["; // %^&
              ss << newi2[i]; // %^&
              ss << "]"; // %^&
              int_node = ss.str();

              ss2 << "(";

              if((newick_parts[i]).size() >2)
              {
                for (unsigned int j = 0 ; j < ((newick_parts[i]).size()-2) ; j ++)
                {
                  ss2 << newick_parts[i][j] << ":" << last_split_dates[i][j] - join_dates[i][j] << ",(";
                }
              }

              thislastsplit = join_dates[i][0];

              ss2 << newick_parts[i][((newick_parts[i]).size()-2)] << ":" << (last_split_dates[i][((newick_parts[i]).size()-2)] - join_dates[i][((join_dates[i]).size()-2)]);
              ss2 << ",";
              ss2 << newick_parts[i][((newick_parts[i]).size()-1)] << ":" << last_split_dates[i][((newick_parts[i]).size()-1)] - join_dates[i][((join_dates[i]).size()-2)];

              for (int j = ((newick_parts[i]).size()-3) ; j >= 0  ; j --)
              {
                ss2 << ")";
                ss2 << int_node;
                ss2 << ":";
                ss2 << ((join_dates[i][j+1])-(join_dates[i][j]));
              }
              ss2 << ")" << int_node ;
            }
            else
            {
              ss2 << newick_parts[i][0];
              thislastsplit = m_current_generation; // this can only really happen when were're at the tip of a tree
            }

            num_descend_done[i] = -1;
            // next sort out and move up - hand to parents
            // if it's not possible to move up return
            if (intermediate_parents[i] == -1)
            {
              final_newick = ss2.str();
            }
            else
            {
              newick_parts[intermediate_parents[i]].push_back(ss2.str());
              join_dates[intermediate_parents[i]].push_back(divergence_date[i]);
              num_descend_done[intermediate_parents[i]] += m_n_descend_spec[i];
              last_split_dates[intermediate_parents[i]].push_back(thislastsplit);
            }

          }
        }
      }
    }


    for (unsigned long i = 0 ; i < m_parents.size() ; ++i)
    {
      if (m_abundances[i] > 0)
      {
        // page needs updating
        if (m_n_descend_spec[i] > 1)
        {
          // it's the minimim of all the other species that join it
          page[i] = (m_current_generation - join_dates[i][((join_dates[i]).size()-2)]);
        }
        else
        {
          page[i] = m_current_generation - divergence_date[i];
          if (m_abundances[intermediate_parents[i]]<=0)
          {
            // it joins an exitinct species
            if (last_split_dates[i][((newick_parts[i]).size()-1)] == divergence_date[i])
            {
              // it was the last to join need to make a change
              page[i] = m_current_generation - (last_split_dates[i][((newick_parts[i]).size()-2)]);
            }
          }
        }
      }
    }

    m_FILE_final_newick.push_back(final_newick);

  }

  long topindex = -1;

  if (output_var)
  {

    // this code will produce the species variance data and store it

    // DATA TO BE PRODUCED

    // std::vector<double> m_FILE_generation_var; // will line up with the variance output rows for frequency of output here
    // std::vector<long> m_FILE_n; // value of n that applies for the rest of the readings (the first three above will probably repeat a lot)
    // std::vector<long> m_FILE_n_richness; // the species richness overall at this value of n
    // std::vector<long> m_FILE_num_subspec; // number of sub species - for this good species only
    //vector<long> m_FILE_tot_abund; // total abundance of all sub species - for this good species only
    // std::vector<double> m_FILE_meanc; // mean value of c - for this good species only
    // std::vector<double> m_FILE_varc; // variange in c - for this good species only
    // std::vector<long> m_FILE_maxc; // max value in c - for this good species only
    // std::vector<long> m_FILE_minc; // min value in c - for this good species only


    /*
           new data:

           std::vector<double> m_FILE_generation_var; // will line up with the variance output rows for frequency of output here
           std::vector<long> m_FILE_n; // value of n that applies for the rest of the readings (the first three above will probably repeat a lot)
           std::vector<long> m_FILE_n_richness; // the species richness overall at this value of n
           std::vector<long> m_FILE_num_subspec; // number of sub species - average over all species
           std::vector<long> m_FILE_tot_abund; // total abundance of all sub species - average over all species
           std::vector<double> m_FILE_meanc; // mean value of c - (within species) average over all species
           std::vector<double> m_FILE_varc; // variange in c - (within species) average over all species
           std::vector<long> m_FILE_rangec; // range in value of c - (within species) average over all species
           std::vector<double> m_FILE_meanc_w; // mean value of c - (within species) weighted average over all species
           std::vector<double> m_FILE_varc_w; // variange in c - (within species) weighted average over all species
           std::vector<long> m_FILE_rangec_w; // range in value of c - (within species) weighted average over all species


           */


    // DATA WE HAVE AVAILABLE TO US

    // std::vector<long> abundances;
    // the abundnance of each species
    // std::vector<long> parents;
    // the parent species for each species
    // std::vector<long> selection_level;
    // the selection level of each species (larger numbers mean greater advantage)
    // std::vector<long> num_mut;
    // number of mutations along the branch that connects this to the parent node, these could be positive or negative mutations

    // THE ALGORITHM

    // to do these calculations I'm going to use the species_obj to record details of good species
    // I will prune the tree and maintain a std::vector of good species to one side
    // I will repeat the process for different values of n including infinity

    // to perform this I'll store a number of useful std::vectors that will line up with the std::vectors of all sub species

    // choose a value of n to do this for
    //long this_n = -1; // loop over more values of n



    for (int this_n = -1 ; this_n < 70 ; )
    {

      // initialise a min_num_mutations (not including mutations up to parent) to living species std::vector for each position so that when working up the tree lumping can be done easily
      std::vector<long> min_num_mutations;
      min_num_mutations.clear();

      // initialise a min_num_mutations to living species std::vector for each position including the possibility that the fewest mutations is through the parent
      std::vector<long> min_num_dir;
      min_num_dir.clear();

      // initialise min_num_which which shows us which path gives the route to a living species with the minimum - this takes on the value of the node with the minimum number of mutations path to a living species.
      // Note that a draw doesn't matter either can be chosen as they'll be treated just the same later anyway - either both lumped or both split.
      // Note that the route could be through the parent
      //vector<long> min_num_which; // I don't think we actually need this after all commented out everywhere now
      //min_num_which.clear();

      // initialise a num_child std::vector (down to all direct chidren nodes) for all nodes as a useful tool
      std::vector<long> num_child;
      num_child.clear();

      // initilise a num_done std::vector (tells us how many of the children we have dealt with so far with pruning - a completely pruned branch should have this equal to num_child)
      // after a completely pruned branch has had it's details dealt with and pushed up to its parent the num_done tag changes to -1
      std::vector<long> num_done;
      num_done.clear();

      // initialise a std::vector pruned_spec of  species objects that gets passed up the tree and amalgamated as appropriate - this is stored on the parent of the species
      std::vector<Species> pruned_spec;
      pruned_spec.clear();

      // here is a std::vector that will store the final result
      std::vector<Species> final_spec;
      final_spec.clear();
      // species object std::vector contains the abundance and c values of all sub species and can produce the necessary outputs.



      // now we fill all these std::vectors with their initial conditions.
      for (unsigned long i = 0 ; i < m_abundances.size() ; ++i)
      {
        Species newspec;
        newspec.Setup(m_abundances[i],m_fitness_categories[i]);
        pruned_spec.push_back(newspec);
        num_done.push_back(0);
        num_child.push_back(0);
        if (m_abundances[i] >0)
        {
          min_num_mutations.push_back(0);
          min_num_dir.push_back(0);
          //min_num_which.push_back(-1000); // used to indicate direction is down
        }
        else
        {
          min_num_mutations.push_back(-1);
          min_num_dir.push_back(-1);
          //min_num_which.push_back(-1); // used to indicate unknown so far
        }
      }







      // we need to iterate through the data to set num_child
      for (unsigned long i = 0 ; i < m_abundances.size() ; ++i)
      {
        if (m_parents[i] >= 0)
        {
          num_child[m_parents[i]] += 1;
        }
      }



      // to complete the initialisation we need to set min_num_mutations
      // an easy way to do this is to keep going until there's no change from a whole sweep of the tree
      bool looping = true;
      while(looping)
      {
        looping = false;
        for (unsigned long i = 0 ; i < m_abundances.size() ; ++i)
        {
          if (m_parents[i] >= 0)
          {
            // loop over whole tree and see if parent's min is set correctly
            long parent_min = min_num_mutations[m_parents[i]];
            long this_min = min_num_mutations[i];
            if (parent_min < 0)
            {
              looping = true;
              if (this_min >= 0)
              {
                min_num_mutations[m_parents[i]] = (this_min+m_n_mut[i]);
                min_num_dir[m_parents[i]] = (this_min+m_n_mut[i]);
                //min_num_which[parents[i]] = i;
              }
            }
            else
            {
              if (this_min < 0)
              {
                looping = true;
              }
              else
              {
                if ((this_min+m_n_mut[i]) < parent_min)
                {
                  min_num_mutations[m_parents[i]] = (this_min+m_n_mut[i]);
                  min_num_dir[m_parents[i]] = (this_min+m_n_mut[i]);
                  //min_num_which[parents[i]] = i;
                  looping = true;
                }
              }
            }
          }
        }
      }



      // now we have to sort out min_num_dir and with the possibility that min_num_which will change too
      // min_num_dir could only decrease really.
      looping = true;
      while(looping)
      {
        looping = false;
        for (unsigned long i = 0 ; i < m_abundances.size() ; ++i)
        {
          if (m_parents[i] >= 0)
          {
            if ((min_num_dir[m_parents[i]]+m_n_mut[i]) < min_num_dir[i])
            {
              min_num_dir[i] = (min_num_dir[m_parents[i]]+m_n_mut[i]);
              //min_num_which[i] = parents[i];
              looping = true;
            }
          }
        }
      }


      // now do the actual pruning step after the std::vectors have been set up
      // we're going to keep looping until every node is dealt with
      looping = true;
      while(looping)
      {
        looping = false;
        for (unsigned long i = 0 ; i < m_abundances.size() ; ++i)
        {
          if (m_parents[i] >= 0)
          {
            if (num_done[i] != -1)
            {
              looping = true;
              // we're not done yet keep going
              if (num_done[i] == num_child[i])
              {
                // ready to be pushed up
                // here we deal with a branch and prune it pushing up the summary information
                // we need to know will it be a good species on its own or will it link up to its parent
                if (((min_num_dir[m_parents[i]]+min_num_mutations[i]+m_n_mut[i])  <  this_n)||(0 > this_n )) // negative values mean infinity
                {
                  // lumping

                  pruned_spec[m_parents[i]].Add(pruned_spec[i]);
                  num_done[i] = -1;
                  num_done[m_parents[i]] ++;
                  pruned_spec[i].Clear();
                }
                else
                {
                  // splitting
                  if ((pruned_spec[i]).GetSumAbundances() > 0)
                  {
                    final_spec.push_back(pruned_spec[i]);
                  }
                  num_done[i] = -1;
                  num_done[m_parents[i]] ++;
                  pruned_spec[i].Clear();
                }
              }
            }
          }
          // else nothing to do for this loop around
        }
      }



      // could be that the final parent needs to be added still
      for (unsigned long i = 0 ; i < m_abundances.size() ; ++i)
      {
        if ((m_parents[i] == -1)&&((pruned_spec[i]).GetNumberOfSpecies()>0))
        {
          final_spec.push_back(pruned_spec[i]);
        }
      }


      long num_spec_total = 0;
      long num_ind_total = 0;
      for (unsigned long i = 0 ; i < final_spec.size() ; ++i)
      {
        if ((final_spec[i].GetSumAbundances())>0)
        {
          num_spec_total ++;
          num_ind_total += (final_spec[i].GetSumAbundances());
        }
      }

      double temp_num_subspec = 0;
      double temp_tot_abund = 0;

      double temp_meanc = 0;
      double temp_varc = 0;
      double temp_rangec = 0;

      double temp_meanc_w = 0;
      double temp_varc_w = 0;
      double temp_rangec_w = 0;

      for (unsigned long i = 0 ; i < final_spec.size() ; ++i)
      {
        if ((final_spec[i].GetSumAbundances())>0)
        {

          temp_num_subspec += (double(final_spec[i].GetNumberOfSpecies()));
          temp_tot_abund += (double(final_spec[i].GetSumAbundances()));

          temp_meanc += (double(final_spec[i].GetFitnessMean()));
          temp_varc += (double(final_spec[i].GetFitnessVariance()));
          temp_rangec += (double(final_spec[i].GetFitnessMax()-final_spec[i].GetFitnessMin()));

          temp_meanc_w += (double(final_spec[i].GetSumAbundances())*double(final_spec[i].GetFitnessMean()));
          temp_varc_w += (double(final_spec[i].GetSumAbundances())*double(final_spec[i].GetFitnessVariance()));
          temp_rangec_w += (double(final_spec[i].GetSumAbundances())*double(final_spec[i].GetFitnessMax()-final_spec[i].GetFitnessMin()));

        }

      }

      m_FILE_generation_var.push_back(m_true_generation); // will line up with the variance output rows for frequency of output here
      m_FILE_n.push_back(this_n); // value of n that applies for the rest of the readings (the first three above will probably repeat a lot)
      m_FILE_n_richness.push_back(num_spec_total); // the species richness overall at this value of n

      m_FILE_num_subspec.push_back(temp_num_subspec/double(num_spec_total));
      m_FILE_tot_abund.push_back(temp_tot_abund/double(num_spec_total));
      m_FILE_meanc.push_back(temp_meanc/double(num_spec_total));
      m_FILE_varc.push_back(temp_varc/double(num_spec_total));
      m_FILE_rangec.push_back(temp_rangec/double(num_spec_total));

      m_FILE_meanc_w.push_back(temp_meanc_w/double(num_ind_total));
      m_FILE_varc_w.push_back(temp_varc_w/double(num_ind_total));
      m_FILE_rangec_w.push_back(temp_rangec_w/double(num_ind_total));

      if (this_n == -1)
      {
        this_n = 1;
      }
      else
      {
        if (this_n >= 30)
        {
          this_n += 5;
        }
        else
        {
          if (this_n >= 10)
          {
            this_n += 2;
          }
          else
          {
            this_n ++;
          }

        }
      }
    }
  }

  for (unsigned int i = 0 ; i < m_parents.size() ; ++i)
  {
    if (keep_data[i])
    {
      m_FILE_generation.push_back(m_true_generation); // when the data was taken - so replaces simnum
      m_FILE_iname.push_back(alpha_format_auto(newi[i])); // can remove this line later
      m_FILE_iname2.push_back(newi2[i]); // can remove this line later
      m_FILE_index.push_back(newi[i]);
      m_FILE_abundances.push_back(m_abundances[i]);

      if (intermediate_parents[i] >= 0)
      {
        m_FILE_parents.push_back(newi[intermediate_parents[i]]);
      }
      else
      {
        m_FILE_parents.push_back(-1);
        topindex = i;
      }
      m_FILE_born.push_back(m_born[i]);
      m_FILE_died.push_back(m_died[i]);
      m_FILE_selection_level.push_back(m_fitness_categories[i]);
      m_FILE_num_mut.push_back(m_n_mut[i]);
      m_FILE_max_abundance.push_back(m_max_abundance[i]);
      m_FILE_maxabtime1.push_back(m_maxabtime1[i]);
      m_FILE_maxabtime2.push_back(m_maxabtime2[i]);
      m_FILE_num_descend_spec.push_back(m_n_descend_spec[i]);
      m_FILE_divergence_date.push_back(divergence_date[i]);
      m_FILE_relgen.push_back(m_current_generation);
      if (m_abundances[i] > 0)
      {
        m_FILE_tage.push_back(m_current_generation-m_born[i]);
      }
      else
      {
        m_FILE_tage.push_back(-1);
      }
      m_FILE_page.push_back(page[i]) ;
    }
  }

  long tempindex = m_parents[topindex];
  // topindex is the index of the root node %^&

  double diffgen = m_born[topindex];
  m_parents[topindex] = -1;
  m_current_generation -= diffgen;
  for (unsigned int i = 0 ; i < m_parents.size() ; ++i)
  {
    if (m_abundances[i] >= 0)
    {
      m_born[i] -= diffgen;
      if (m_died[i] >= 0)
      {
        m_died[i] -= diffgen;
      }
      m_maxabtime1[i] -= diffgen;
      m_maxabtime2[i] -= diffgen;
    }
  }

  m_born[topindex] = 0.0;

  while (tempindex >= 0)
  {
    long   newtempindex = m_parents[tempindex];

    m_abundances[tempindex] = -2;
    m_max_abundance[tempindex] = -2;
    m_maxabtime1[tempindex] = -2.0;
    m_maxabtime2[tempindex] = -2.0;
    m_parents[tempindex] = -3;
    m_born[tempindex] = -2.0;
    m_died[tempindex] = -2.0;
    m_fitness_categories[tempindex] = -2000000;
    m_n_mut[tempindex]  = -2.0; // %^&
    m_n_descend_spec[tempindex] = -2;

    if (m_end_recycler < static_cast<int>(m_recycler.size()))
    {
      m_recycler[m_end_recycler] = tempindex;
    }
    else
    {
      m_recycler.push_back(tempindex);
    }
    m_end_recycler ++;

    tempindex = newtempindex;

  }

}

void NTsim::set_seed(long seedin)
{
  if (!m_is_seeded)
  {
    m_seed_used = seedin;
    m_is_seeded = true;
    m_rng.set_seed(seedin);
  }
}

void NTsim::setup(
  const long seed,
  const long metacommunity_size,
  const double mutation_rate,
  const double selection_strength
)
{
  // Prepare for a new simulation with the given parameters
  // note there is also a system restore funciton lower down in this file
  total_clear();
  set_seed(seed);

  // Set simulation parmaeters correcty first
  m_metacommunity_size = metacommunity_size;
  m_mutation_rate = mutation_rate;
  m_selection_strength = selection_strength;

  // setup metacommunity
  for (long i = 0 ; i < m_metacommunity_size ; ++i)
  {
    m_metacommunity.push_back(m_upto_which_label);
  }
  m_upto_which_label ++;
  m_current_richness = 1;

  // setup burnin counting variables
  for (long i = 0 ; i < m_metacommunity_size ; ++i)
  {
    m_init_ancestor.push_back(i);
    m_n_descend.push_back(1);
  }
  m_n_remain = m_metacommunity_size;

  // species level data std::vectors to be set with initial conditions
  m_abundances.push_back(m_metacommunity_size);
  m_parents.push_back(-1);
  m_born.push_back(0.0);
  m_died.push_back(-1.0);
  m_fitness_categories.push_back(0);
  m_n_mut.push_back(0); // %^&

  m_max_abundance.push_back(m_metacommunity_size);
  m_maxabtime1.push_back(0.0);
  m_maxabtime2.push_back(0.0);

  m_n_descend_spec.push_back(1);

  // now setup the output file vars

  m_FILE_final_newick.clear();

  m_FILE_generation.clear();
  m_FILE_iname.clear();
  m_FILE_iname2.clear();
  m_FILE_index.clear();
  m_FILE_abundances.clear();
  m_FILE_parents.clear();
  m_FILE_born.clear();
  m_FILE_died.clear();
  m_FILE_selection_level.clear();
  m_FILE_num_mut.clear();
  m_FILE_max_abundance.clear();
  m_FILE_maxabtime1.clear();
  m_FILE_maxabtime2.clear();
  m_FILE_num_descend_spec.clear();
  m_FILE_divergence_date.clear();
  m_FILE_relgen.clear();
  m_FILE_tage.clear();
  m_FILE_page.clear();

  m_FILE_generation_var.clear();
  m_FILE_n.clear();
  m_FILE_n_richness.clear();
  m_FILE_num_subspec.clear();
  m_FILE_tot_abund.clear();
  m_FILE_meanc.clear();
  m_FILE_varc.clear();
  m_FILE_rangec.clear();
  m_FILE_meanc_w.clear();
  m_FILE_varc_w.clear();
  m_FILE_rangec_w.clear();

}

void NTsim::sim_all(
  const long seed,
  const long metacommunity_size ,
  const double mutation_rate ,
  const double selection_strength ,
  const double sim_time_secs)
{
  // initialise timing points
  time_t start , end;
  time(&start);

  bool timeout = false;
  long steps = 0;

  setup(seed,metacommunity_size,mutation_rate,selection_strength);

  while ((m_n_remain > 1 || m_abundances[0] > 0) && !timeout)
  {
    sim_step();
    steps ++;

    if (steps%10000 == 0) //£$% - removed 00
    {

      time(&end);
      if (difftime(end,start) >= (sim_time_secs))
      {
        timeout = true;
      }
    }
  }

  if (!timeout)
  {
    m_gen_burned_in = m_true_generation; // inside the if statement beacuse if we have run out of time then this is not the correct gen_burned_in

    while (!timeout && m_iter <= 1000000) // checks in time and also that not too much HDD space has been used
    {
      sim_step();
      steps ++;

      if (steps%1000000 == 0)
      {
        time(&end);
        if (difftime(end,start) >= (sim_time_secs))
        {
          timeout = true;
        }
      }

      if (m_true_generation >= (m_gen_burned_in*double(m_simnum)/10.0 +4.0)) // 4 times the expected minimum burnin (hence the +4 here - up from 3) the /10 increases the frequency of outputting
      {

        if (m_simnum == 1)
        {
          std::cout << "simulation fully burned in\n";
        }
        output();
        m_simnum ++;
      }
    }
  }
  // do a final output just in case the system needs to be suspended here
  output();
}

void NTsim::sim_restore(const char * const  fnamein, const double simtime)
{
  std::cout << "restoring NNIM from " << fnamein << "\n";
  std::ifstream in;
  in.open(fnamein);

  // scalars

  in >> m_metacommunity_size ;
  in >> m_mutation_rate ;
  in >> m_selection_strength ;
  in >> m_n_suspends ;
  in >> m_n_file_outputs ;
  in >> m_current_generation ;
  in >> m_true_generation ;
  in >> m_gen_burned_in ;
  in >> m_simnum ;
  in >> m_max_selection_level ;
  in >> m_current_richness ;
  in >> m_upto_which_label ;
  in >> m_end_recycler ;
  in >> m_n_remain ;
  in >> m_iter ;
  in >> m_iter2 ;
  in >> m_is_seeded ;
  in >> m_seed_used ;

  // std::vectors

  long temp_size;
  long temp_long;
  double temp_double;

  in >> temp_size;
  m_recycler.clear();
  for (long i = 0 ; i < temp_size ; ++i)
  {
    in >> temp_long;
    m_recycler.push_back(temp_long);
  }

  in >> temp_size;
  m_init_ancestor.clear();
  for (long i = 0 ; i < temp_size ; ++i)
  {
    in >> temp_long;
    m_init_ancestor.push_back(temp_long);
  }

  in >> temp_size;
  m_n_descend.clear();
  for (long i = 0 ; i < temp_size ; ++i)
  {
    in >> temp_long;
    m_n_descend.push_back(temp_long);
  }

  in >> temp_size;
  m_metacommunity.clear();
  for (long i = 0 ; i < temp_size ; ++i)
  {
    in >> temp_long;
    m_metacommunity.push_back(temp_long);
  }

  in >> temp_size;
  m_abundances.clear();
  for (long i = 0 ; i < temp_size ; ++i)
  {
    in >> temp_long;
    m_abundances.push_back(temp_long);
  }

  in >> temp_size;
  m_parents.clear();
  for (long i = 0 ; i < temp_size ; ++i)
  {
    in >> temp_long;
    m_parents.push_back(temp_long);
  }

  in >> temp_size;
  m_born.clear();
  for (long i = 0 ; i < temp_size ; ++i)
  {
    in >> temp_double;
    m_born.push_back(temp_double);
  }

  in >> temp_size;
  m_died.clear();
  for (long i = 0 ; i < temp_size ; ++i)
  {
    in >> temp_double;
    m_died.push_back(temp_double);
  }

  in >> temp_size;
  m_fitness_categories.clear();
  for (long i = 0 ; i < temp_size ; ++i)
  {
    in >> temp_long;
    m_fitness_categories.push_back(temp_long);
  }

  in >> temp_size;
  m_n_mut.clear();
  for (long i = 0 ; i < temp_size ; ++i)
  {
    in >> temp_long;
    m_n_mut.push_back(temp_long);
  }

  in >> temp_size;
  m_max_abundance.clear();
  for (long i = 0 ; i < temp_size ; ++i)
  {
    in >> temp_long;
    m_max_abundance.push_back(temp_long);
  }

  in >> temp_size;
  m_maxabtime1.clear();
  for (long i = 0 ; i < temp_size ; ++i)
  {
    in >> temp_double;
    m_maxabtime1.push_back(temp_double);
  }

  in >> temp_size;
  m_maxabtime2.clear();
  for (long i = 0 ; i < temp_size ; ++i)
  {
    in >> temp_double;
    m_maxabtime2.push_back(temp_double);
  }

  in >> temp_size;
  m_n_descend_spec.clear();
  for (long i = 0 ; i < temp_size ; ++i)
  {
    in >> temp_long;
    m_n_descend_spec.push_back(temp_long);
  }

  std::vector<long> NR_suspend;

  in >> temp_size;
  NR_suspend.clear();
  for (long i = 0 ; i < temp_size ; ++i)
  {
    in >> temp_long;
    NR_suspend.push_back(temp_long);
  }

  m_rng.resume(NR_suspend);

  in.close();

  std::cout << "system restored \n";

  // initialise timing points
  time_t start , end;
  time(&start);

  bool timeout = false;
  long steps = 0;

  if (m_gen_burned_in < 0)
  {
    while (((m_n_remain > 1)||(m_abundances[0] > 0))&&(!timeout))
    {
      sim_step();
      steps ++;

      if (steps%1000000 == 0)
      {
        time(&end);
        if (difftime(end,start) >= (simtime))
        {
          timeout = true;
        }
      }

    }
    if (!timeout)
    {
      m_gen_burned_in = m_true_generation;
    }
  }

  std::cout << "simulation part burned in\n";

  if (!timeout)
  {

    while ((!timeout)&&(m_iter <= 1000000)) // checks in time and also that not too much HDD space has been used
    {
      sim_step();
      steps ++;

      if (steps%1000000 == 0)
      {
        time(&end);
        if (difftime(end,start) >= (simtime))
        {
          timeout = true;
        }
      }

      if (m_true_generation >= (m_gen_burned_in*double(m_simnum)/10.0 +4.0)) // 4 times the expected minimum burnin (hence the +3 here) the /10 increase the frequency of outputting
      {

        if (m_simnum == 1)
        {
          std::cout << "simulation fully burned in\n";
        }

        output();
        m_simnum ++;

      }
    }
  }

}

void NTsim::sim_step()
{

  // increment generation counter
  m_current_generation += 2.0/m_metacommunity_size;
  m_true_generation += 2.0/m_metacommunity_size;

  // choose individual to die (and be replaced)
  long chosen2;
  chosen2 = m_rng.i0(m_metacommunity_size-1);

  // decide if speciation happened
  bool speciation = false;
  if (m_rng.d01() <= m_mutation_rate)
  {
    speciation = true;
  }

  // chooose replacement individual
  long chosen;
  bool repeat_selection = true;
  while(repeat_selection)
  {
    // keep going until selection criteria are satisfied
    chosen = m_rng.i0(m_metacommunity_size-1);
    if (chosen != chosen2)
    {
      if (m_selection_strength > 0)
      {
        // deals with selection
        long temp_selection = m_fitness_categories[m_metacommunity[chosen]];
        if (speciation)
        {
          // this section is for additive selective advantage
          /*
                       if (NR.d01()<=(1.0/(1.0+double(max_selection_level - temp_selection -1)*s)))
                       {
                       repeat_selection = false;
                       }
                       // note: if max selection is about to be increased because the focal species is the fittest
                       // the if statement will be true anyway and the loop will quit - the RHS of the inequality will be greater than 1
                       //*/

          // this section is for multiplicative selective advantage
          //*
          if (m_rng.d01()<=(pow((1.0+m_selection_strength),-1.0*double(m_max_selection_level - temp_selection -1))))
          {
            repeat_selection = false;
          }
          // note: if max selection is about to be increased because the focal species is the fittest
          // the if statement will be true anyway and the loop will quit - the RHS of the inequality will be greater than 1
          //*/
        }
        else
        {
          // this section is for additive selective advantage
          /*
                       if (NR.d01()<=(1.0/(1.0+double(max_selection_level - temp_selection)*s)))
                       {
                       repeat_selection = false;
                       }
                       //*/

          // this section is for multiplicative selective advantage
          //*
          if (m_rng.d01()<=(pow((1.0+m_selection_strength),-1.0*double(m_max_selection_level - temp_selection))))
          {
            repeat_selection = false;
          }
          //*/
        }
      }
      else
      {
        repeat_selection = false;
      }
    }
  }

  // update burn in checker only if we're not burned in yet
  if (m_gen_burned_in<0)
  {
    m_n_descend[m_init_ancestor[chosen]]++;
    m_n_descend[m_init_ancestor[chosen2]]--;
    if (m_n_descend[m_init_ancestor[chosen2]] == 0)
    {
      m_n_remain --;
    }
    m_init_ancestor[chosen2] = m_init_ancestor[chosen];
  }

  // now update the reaminder of the variables based on
  // speciation, chosen and chosen2

  // update main data

  m_abundances[m_metacommunity[chosen2]]--;
  if (m_abundances[m_metacommunity[chosen2]] == 0)
  {

    m_died[m_metacommunity[chosen2]] = m_current_generation;
    m_current_richness --;

    long tempindex = m_metacommunity[chosen2]; // the species index of the now dead species
    bool keeplooping = true;
    while(keeplooping)
    {
      m_n_descend_spec[tempindex] --;

      long newtempindex = m_parents[tempindex]; // move up to see if more can be recycled now

      if (m_n_descend_spec[tempindex] == 0)
      {
        // recycle

        // set flags just to be sure it's recycled (good for redundancy in case of bugs so data that slipps through is obviously wrong)
        m_abundances[tempindex] = -2;
        m_max_abundance[tempindex] = -2;
        m_maxabtime1[tempindex] = -2.0;
        m_maxabtime2[tempindex] = -2.0;
        m_parents[tempindex] = -3;
        m_born[tempindex] = -2.0;
        m_died[tempindex] = -2.0;
        m_fitness_categories[tempindex] = -2000000;
        m_n_mut[tempindex] = -2;
        m_n_descend_spec[tempindex] = -2;

        if (m_end_recycler < static_cast<int>(m_recycler.size()))
        {
          m_recycler[m_end_recycler] = tempindex;
        }
        else
        {
          m_recycler.push_back(tempindex);
        }
        m_end_recycler ++;

      }

      tempindex = newtempindex;
      if (tempindex < 0)
      {
        keeplooping = false;
      }
    }
  }

  if (speciation)
  {
    if (m_end_recycler == 0)
    {
      // no species numbers available to recycle
      m_metacommunity[chosen2] = m_upto_which_label;
      m_abundances.push_back(1);
      m_max_abundance.push_back(1);
      m_maxabtime1.push_back(m_current_generation);
      m_maxabtime2.push_back(m_current_generation);

      m_parents.push_back(m_metacommunity[chosen]);
      m_born.push_back(m_current_generation);
      m_died.push_back(-1.0);

      m_n_mut.push_back(1);

      if (m_rng.d01() < 0.5)
      {
        m_fitness_categories.push_back(m_fitness_categories[m_metacommunity[chosen]]-1);
      }
      else
      {
        m_fitness_categories.push_back(m_fitness_categories[m_metacommunity[chosen]]+1);
        if ((m_fitness_categories[m_metacommunity[chosen]]+1)>m_max_selection_level)
        {
          m_max_selection_level = (m_fitness_categories[m_metacommunity[chosen]]+1);
        }
      }

      m_n_descend_spec.push_back(1);

      m_upto_which_label ++;
    }
    else
    {
      // species are available to recycle
      long tempindex = m_recycler[m_end_recycler-1];

      m_metacommunity[chosen2] = tempindex;
      m_abundances[tempindex] = 1;
      m_max_abundance[tempindex] = 1;
      m_maxabtime1[tempindex] = m_current_generation;
      m_maxabtime2[tempindex] = m_current_generation;

      m_parents[tempindex] = m_metacommunity[chosen];
      m_born[tempindex] = m_current_generation;
      m_died[tempindex] = -1.0;

      m_n_mut[tempindex] = 1;

      if (m_rng.d01() < 0.5)
      {
        m_fitness_categories[tempindex] = (m_fitness_categories[m_metacommunity[chosen]]-1);
      }
      else
      {
        m_fitness_categories[tempindex] = (m_fitness_categories[m_metacommunity[chosen]]+1);

        if ((m_fitness_categories[tempindex])>m_max_selection_level)
        {
          m_max_selection_level = m_fitness_categories[tempindex];
        }

      }

      m_n_descend_spec[tempindex] = 1;

      m_end_recycler --;
    }

    long tempindex = m_metacommunity[chosen];
    while(tempindex >= 0)
    {
      m_n_descend_spec[tempindex] ++;
      tempindex = m_parents[tempindex]; // move up to see if more can be recycled now
    }

    m_current_richness ++;

  }
  else
  {
    // no speciation
    m_abundances[m_metacommunity[chosen]]++;
    if ((m_abundances[m_metacommunity[chosen]]) == m_max_abundance[m_metacommunity[chosen]])
    {
      m_maxabtime2[m_metacommunity[chosen]] = m_current_generation;
    }
    if ((m_abundances[m_metacommunity[chosen]]) > m_max_abundance[m_metacommunity[chosen]])
    {
      m_maxabtime2[m_metacommunity[chosen]] = m_current_generation;
      m_maxabtime1[m_metacommunity[chosen]] = m_current_generation;
      m_max_abundance[m_metacommunity[chosen]] = m_abundances[m_metacommunity[chosen]];
    }

    m_metacommunity[chosen2] = m_metacommunity[chosen];
  }
}

void NTsim::total_clear()
{
  m_current_generation = 0.0;
  m_true_generation = 0.0;
  m_gen_burned_in = -1.0;
  m_simnum = 1;
  m_max_selection_level = 0;
  m_current_richness = 0;

  m_upto_which_label = 0;
  m_recycler.clear();
  m_end_recycler = 0;

  m_init_ancestor.clear();
  m_n_descend.clear();
  m_n_remain = -1;

  m_metacommunity.clear();

  m_abundances.clear();
  m_parents.clear();
  m_born.clear();
  m_died.clear();
  m_fitness_categories.clear();
  m_n_mut.clear();

  m_max_abundance.clear();
  m_maxabtime1.clear();
  m_maxabtime2.clear();

  m_n_descend_spec.clear();
  // random number generator, simulation parameters are all set elsewhere
}



void NTsim::write_files() const
{
  std::cout << "writing files \n";

  std::ofstream out;
  out.precision(10);

  if (m_FILE_generation.size() > 0) // check that there is anything to output
  {
    m_n_file_outputs ++;

    // now setup the output file vars
    // species level file
    char filename_ab[100];
    sprintf (filename_ab, "Data_%i_FV_Out_%i.csv", int(m_seed_used),int(m_n_file_outputs)); // edited here to make csv &&&&&&&
    out.open(filename_ab);
    // output just the header now
    out <<  "generation,iname,iname2,index,abundances,parents,born,died,selection_level,num_mut,max_abundance,maxabtime1,maxabtime2,num_descend_spec,divergence_date,relgen,tage,page\n";

    for (unsigned long i = 0 ; i < m_FILE_generation.size() ; ++i)
    {
      out <<  m_FILE_generation[i] <<",";
      out <<  m_FILE_iname[i] <<",";
      out <<  m_FILE_iname2[i] <<",";
      out <<  m_FILE_index[i] <<",";
      out <<  m_FILE_abundances[i] <<",";
      out <<  m_FILE_parents[i] <<",";
      out <<  m_FILE_born[i] <<",";
      out <<  m_FILE_died[i] <<",";
      out <<  m_FILE_selection_level[i] <<",";
      out <<  m_FILE_num_mut[i] <<",";
      out <<  m_FILE_max_abundance[i] <<",";
      out <<  m_FILE_maxabtime1[i] <<",";
      out <<  m_FILE_maxabtime2[i] <<",";
      out <<  m_FILE_num_descend_spec[i] <<",";
      out <<  m_FILE_divergence_date[i] <<",";
      out <<  m_FILE_relgen[i] <<",";
      out <<  m_FILE_tage[i] <<",";
      out <<  m_FILE_page[i] <<",\n";
    }
    out.close();

    // tree file
    char filename_newick[100];
    sprintf (filename_newick, "Data_%i_FV_Newick_%i.txt", int(m_seed_used),int(m_n_file_outputs));
    out.open(filename_newick);
    for (unsigned long i = 0 ; i < m_FILE_final_newick.size() ; ++i)
    {
      out << m_FILE_final_newick[i] << ";\n";
    }
    out.close();

    // variance file
    sprintf (filename_ab, "Data_%i_FV_FitVar_%i.csv", int(m_seed_used),int(m_n_file_outputs)); // edited here to make csv &&&&&&&
    out.open(filename_ab);
    // output just the header now
    out <<  "generation_var,n,n_richness,num_subspec,tot_abund,meanc,varc,rangec,meanc_w,varc_w,rangec_w\n";

    for (unsigned long i = 0 ; i < m_FILE_generation_var.size() ; ++i)
    {
      out <<  m_FILE_generation_var[i] <<",";
      out <<  m_FILE_n[i] <<",";
      out <<  m_FILE_n_richness[i] <<",";
      out <<  m_FILE_num_subspec[i] <<",";
      out <<  m_FILE_tot_abund[i] <<",";
      out <<  m_FILE_meanc[i] <<",";
      out <<  m_FILE_varc[i] <<",";
      out <<  m_FILE_rangec[i] <<",";
      out <<  m_FILE_meanc_w[i] <<",";
      out <<  m_FILE_varc_w[i] <<",";
      out <<  m_FILE_rangec_w[i] <<",\n";

    }
    out.close();
  }

  if (m_n_suspends == 0)
  {
    char filename_vars[100];
    sprintf (filename_vars, "Data_%i_FV_Vars.txt", int(m_seed_used));
    out.open(filename_vars);
    out << m_seed_used << " \n";
    out << m_metacommunity_size << " \n";
    out << m_mutation_rate << " \n";
    out << m_selection_strength << " \n";
    out.close();
  }

  m_n_suspends ++;

  char filename_vars[100];
  sprintf (filename_vars, "Data_%i_FV_Res_%i.txt", int(m_seed_used) ,int(m_n_suspends));
  out.open(filename_vars);
  // strore the simulation state so we can continue another time if necessary

  out.precision(15);

  // scalars

  out << m_metacommunity_size << " ";
  out << m_mutation_rate << " ";
  out << m_selection_strength << " ";
  out << m_n_suspends << " ";
  out << m_n_file_outputs << " ";
  out << m_current_generation << " ";
  out << m_true_generation << " ";
  out << m_gen_burned_in << " ";
  out << m_simnum << " ";
  out << m_max_selection_level << " ";
  out << m_current_richness << " ";
  out << m_upto_which_label << " ";
  out << m_end_recycler << " ";
  out << m_n_remain << " ";
  out << m_iter << " ";
  out << m_iter2 << " ";
  out << m_is_seeded << " ";
  out << m_seed_used << " ";

  // std::vectors

  out << m_recycler.size() << " ";
  for (unsigned long i = 0 ; i < m_recycler.size() ; ++i)
  {
    out << m_recycler[i] << " ";
  }

  out << m_init_ancestor.size() << " ";
  for (unsigned long i = 0 ; i < m_init_ancestor.size() ; ++i)
  {
    out << m_init_ancestor[i] << " ";
  }

  out << m_n_descend.size() << " ";
  for (unsigned long i = 0 ; i < m_n_descend.size() ; ++i)
  {
    out << m_n_descend[i] << " ";
  }

  out << m_metacommunity.size() << " ";
  for (unsigned long i = 0 ; i < m_metacommunity.size() ; ++i)
  {
    out << m_metacommunity[i] << " ";
  }

  out << m_abundances.size() << " ";
  for (unsigned long i = 0 ; i < m_abundances.size() ; ++i)
  {
    out << m_abundances[i] << " ";
  }

  out << m_parents.size() << " ";
  for (unsigned long i = 0 ; i < m_parents.size() ; ++i)
  {
    out << m_parents[i] << " ";
  }

  out << m_born.size() << " ";
  for (unsigned long i = 0 ; i < m_born.size() ; ++i)
  {
    out << m_born[i] << " ";
  }

  out << m_died.size() << " ";
  for (unsigned long i = 0 ; i < m_died.size() ; ++i)
  {
    out << m_died[i] << " ";
  }

  out << m_fitness_categories.size() << " ";
  for (unsigned long i = 0 ; i < m_fitness_categories.size() ; ++i)
  {
    out << m_fitness_categories[i] << " ";
  }

  out << m_n_mut.size() << " ";
  for (unsigned long i = 0 ; i < m_n_mut.size() ; ++i)
  {
    out << m_n_mut[i] << " ";
  }

  out << m_max_abundance.size() << " ";
  for (unsigned long i = 0 ; i < m_max_abundance.size() ; ++i)
  {
    out << m_max_abundance[i] << " ";
  }

  out << m_maxabtime1.size() << " ";
  for (unsigned long i = 0 ; i < m_maxabtime1.size() ; ++i)
  {
    out << m_maxabtime1[i] << " ";
  }

  out << m_maxabtime2.size() << " ";
  for (unsigned long i = 0 ; i < m_maxabtime2.size() ; ++i)
  {
    out << m_maxabtime2[i] << " ";
  }

  out << m_n_descend_spec.size() << " ";
  for (unsigned long i = 0 ; i < m_n_descend_spec.size() ; ++i)
  {
    out << m_n_descend_spec[i] << " ";
  }

  std::vector<long> NR_suspend = m_rng.suspend();
  out << NR_suspend.size() << " ";
  for (unsigned long i = 0 ; i < NR_suspend.size() ; ++i)
  {
    out << NR_suspend[i] << " ";
  }

  out.close();
}
