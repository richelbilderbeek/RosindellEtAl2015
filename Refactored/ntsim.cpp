#include "ntsim.h"

#include <cassert>

#include "helper.h"

NTsim::NTsim()
{
  total_clear();
  seeded = false;
  iter = 0;
  num_suspends = 0;
  num_file_outputs = 0;
}



void NTsim::output()
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
  std::vector<long> num_mut2 = num_mut;
  keep_data.clear();
  divergence_date.clear();
  intermediate_parents.clear();
  for (unsigned int i = 0 ; i < parents.size() ; ++i)
  {
    keep_data.push_back(true);
    divergence_date.push_back(-1.0);
    intermediate_parents.push_back(parents[i]);
  }
  for (unsigned int i = 0 ; i < parents.size() ; ++i)
  {
    if (abundances[i] < 0)
    {
      keep_data[i] = false;
    }
    if (keep_data[i])
    {
      long current_index = i;
      while(((current_index >= 0)&&(parents[current_index] >= 0))&&(num_descend_spec[parents[current_index]]==num_descend_spec[i]))
      {
        current_index = parents[current_index];
        if (current_index >= 0)
        {
          keep_data[current_index] = false;
          num_mut[i] += num_mut2[current_index]; // accumulate the mutations of all the intermediate link species that only have one descendent species %^&
        }
      }
      if (parents[current_index] >= 0)
      {
        divergence_date[i] = born[current_index];
        intermediate_parents[i] = parents[current_index];
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

  for (unsigned int i = 0 ; i < parents.size() ; ++i)
  {
    if (keep_data[i])
    {
      newi.push_back(iter);
      iter ++;
    }
    else
    {
      newi.push_back(-1);
    }
  }

  iter2 = 1; // note - I added iter 2 above to count the number of lines in this simulation on its own
  std::vector<long> newi2; // counts from 1 up in each simulation on its own
  newi2.clear();

  for (unsigned int i = 0 ; i < parents.size() ; ++i)
  {
    if (keep_data[i])
    {
      newi2.push_back(iter2);
      iter2 ++;
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

    for (unsigned int i = 0 ; i < parents.size() ; ++i)
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

    for (unsigned int i = 0 ; i < parents.size() ; ++i)
    {
      if (abundances[i] > 0) // automatically excludes any lineages that might not be kept
      {
        // we're at the end
        std::ostringstream ss;
        ss.precision(10);
        ss << alpha_format_auto(newi[i]);
        ss << "["; // %^&
        ss << newi2[i]; // %^&
        ss << "]"; // %^&
        newick_parts[i].push_back(ss.str());
        last_split_dates[i].push_back(generation);
        join_dates[i].push_back(generation);
        num_descend_done[i] += 1;
      }
    }

    bool keep_looping = true;
    while(keep_looping)
    {

      keep_looping = false;
      for (unsigned long i = 0 ; i < parents.size() ; ++i)
      {
        // check if we're done
        if ((num_descend_done[i] >= 0)&&(keep_data[i])) // num_descend done goes negative when we're done
        {
          keep_looping = true;

          if (num_descend_done[i] == num_descend_spec[i])
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
              thislastsplit = generation; // this can only really happen when were're at the tip of a tree
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
              num_descend_done[intermediate_parents[i]] += num_descend_spec[i];
              last_split_dates[intermediate_parents[i]].push_back(thislastsplit);
            }

          }
        }
      }
    }


    for (unsigned long i = 0 ; i < parents.size() ; ++i)
    {
      if (abundances[i] > 0)
      {
        // page needs updating
        if (num_descend_spec[i] > 1)
        {
          // it's the minimim of all the other species that join it
          page[i] = (generation - join_dates[i][((join_dates[i]).size()-2)]);
        }
        else
        {
          page[i] = generation - divergence_date[i];
          if (abundances[intermediate_parents[i]]<=0)
          {
            // it joins an exitinct species
            if (last_split_dates[i][((newick_parts[i]).size()-1)] == divergence_date[i])
            {
              // it was the last to join need to make a change
              page[i] = generation - (last_split_dates[i][((newick_parts[i]).size()-2)]);
            }
          }
        }
      }
    }

    FILE_final_newick.push_back(final_newick);

  }

  long topindex = -1;

  if (output_var)
  {

    // this code will produce the species variance data and store it

    // DATA TO BE PRODUCED

    // std::vector<double> FILE_generation_var; // will line up with the variance output rows for frequency of output here
    // std::vector<long> FILE_n; // value of n that applies for the rest of the readings (the first three above will probably repeat a lot)
    // std::vector<long> FILE_n_richness; // the species richness overall at this value of n
    // std::vector<long> FILE_num_subspec; // number of sub species - for this good species only
    //vector<long> FILE_tot_abund; // total abundance of all sub species - for this good species only
    // std::vector<double> FILE_meanc; // mean value of c - for this good species only
    // std::vector<double> FILE_varc; // variange in c - for this good species only
    // std::vector<long> FILE_maxc; // max value in c - for this good species only
    // std::vector<long> FILE_minc; // min value in c - for this good species only


    /*
           new data:

           std::vector<double> FILE_generation_var; // will line up with the variance output rows for frequency of output here
           std::vector<long> FILE_n; // value of n that applies for the rest of the readings (the first three above will probably repeat a lot)
           std::vector<long> FILE_n_richness; // the species richness overall at this value of n
           std::vector<long> FILE_num_subspec; // number of sub species - average over all species
           std::vector<long> FILE_tot_abund; // total abundance of all sub species - average over all species
           std::vector<double> FILE_meanc; // mean value of c - (within species) average over all species
           std::vector<double> FILE_varc; // variange in c - (within species) average over all species
           std::vector<long> FILE_rangec; // range in value of c - (within species) average over all species
           std::vector<double> FILE_meanc_w; // mean value of c - (within species) weighted average over all species
           std::vector<double> FILE_varc_w; // variange in c - (within species) weighted average over all species
           std::vector<long> FILE_rangec_w; // range in value of c - (within species) weighted average over all species


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
      for (unsigned long i = 0 ; i < abundances.size() ; ++i)
      {
        Species newspec;
        newspec.Setup(abundances[i],selection_level[i]);
        pruned_spec.push_back(newspec);
        num_done.push_back(0);
        num_child.push_back(0);
        if (abundances[i] >0)
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
      for (unsigned long i = 0 ; i < abundances.size() ; ++i)
      {
        if (parents[i] >= 0)
        {
          num_child[parents[i]] += 1;
        }
      }



      // to complete the initialisation we need to set min_num_mutations
      // an easy way to do this is to keep going until there's no change from a whole sweep of the tree
      bool looping = true;
      while(looping)
      {
        looping = false;
        for (unsigned long i = 0 ; i < abundances.size() ; ++i)
        {
          if (parents[i] >= 0)
          {
            // loop over whole tree and see if parent's min is set correctly
            long parent_min = min_num_mutations[parents[i]];
            long this_min = min_num_mutations[i];
            if (parent_min < 0)
            {
              looping = true;
              if (this_min >= 0)
              {
                min_num_mutations[parents[i]] = (this_min+num_mut[i]);
                min_num_dir[parents[i]] = (this_min+num_mut[i]);
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
                if ((this_min+num_mut[i]) < parent_min)
                {
                  min_num_mutations[parents[i]] = (this_min+num_mut[i]);
                  min_num_dir[parents[i]] = (this_min+num_mut[i]);
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
        for (unsigned long i = 0 ; i < abundances.size() ; ++i)
        {
          if (parents[i] >= 0)
          {
            if ((min_num_dir[parents[i]]+num_mut[i]) < min_num_dir[i])
            {
              min_num_dir[i] = (min_num_dir[parents[i]]+num_mut[i]);
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
        for (unsigned long i = 0 ; i < abundances.size() ; ++i)
        {
          if (parents[i] >= 0)
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
                if (((min_num_dir[parents[i]]+min_num_mutations[i]+num_mut[i])  <  this_n)||(0 > this_n )) // negative values mean infinity
                {
                  // lumping

                  pruned_spec[parents[i]].Add(pruned_spec[i]);
                  num_done[i] = -1;
                  num_done[parents[i]] ++;
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
                  num_done[parents[i]] ++;
                  pruned_spec[i].Clear();
                }
              }
            }
          }
          // else nothing to do for this loop around
        }
      }



      // could be that the final parent needs to be added still
      for (unsigned long i = 0 ; i < abundances.size() ; ++i)
      {
        if ((parents[i] == -1)&&((pruned_spec[i]).GetNumberOfSpecies()>0))
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

          temp_meanc += (double(final_spec[i].GetTraitMean()));
          temp_varc += (double(final_spec[i].GetTraitVariance()));
          temp_rangec += (double(final_spec[i].GetTraitMax()-final_spec[i].GetTraitMin()));

          temp_meanc_w += (double(final_spec[i].GetSumAbundances())*double(final_spec[i].GetTraitMean()));
          temp_varc_w += (double(final_spec[i].GetSumAbundances())*double(final_spec[i].GetTraitVariance()));
          temp_rangec_w += (double(final_spec[i].GetSumAbundances())*double(final_spec[i].GetTraitMax()-final_spec[i].GetTraitMin()));

        }

      }

      FILE_generation_var.push_back(true_generation); // will line up with the variance output rows for frequency of output here
      FILE_n.push_back(this_n); // value of n that applies for the rest of the readings (the first three above will probably repeat a lot)
      FILE_n_richness.push_back(num_spec_total); // the species richness overall at this value of n

      FILE_num_subspec.push_back(temp_num_subspec/double(num_spec_total));
      FILE_tot_abund.push_back(temp_tot_abund/double(num_spec_total));
      FILE_meanc.push_back(temp_meanc/double(num_spec_total));
      FILE_varc.push_back(temp_varc/double(num_spec_total));
      FILE_rangec.push_back(temp_rangec/double(num_spec_total));

      FILE_meanc_w.push_back(temp_meanc_w/double(num_ind_total));
      FILE_varc_w.push_back(temp_varc_w/double(num_ind_total));
      FILE_rangec_w.push_back(temp_rangec_w/double(num_ind_total));

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

  for (unsigned int i = 0 ; i < parents.size() ; ++i)
  {
    if (keep_data[i])
    {
      FILE_generation.push_back(true_generation); // when the data was taken - so replaces simnum
      FILE_iname.push_back(alpha_format_auto(newi[i])); // can remove this line later
      FILE_iname2.push_back(newi2[i]); // can remove this line later
      FILE_index.push_back(newi[i]);
      FILE_abundances.push_back(abundances[i]);

      if (intermediate_parents[i] >= 0)
      {
        FILE_parents.push_back(newi[intermediate_parents[i]]);
      }
      else
      {
        FILE_parents.push_back(-1);
        topindex = i;
      }
      FILE_born.push_back(born[i]);
      FILE_died.push_back(died[i]);
      FILE_selection_level.push_back(selection_level[i]);
      FILE_num_mut.push_back(num_mut[i]);
      FILE_max_abundance.push_back(max_abundance[i]);
      FILE_maxabtime1.push_back(maxabtime1[i]);
      FILE_maxabtime2.push_back(maxabtime2[i]);
      FILE_num_descend_spec.push_back(num_descend_spec[i]);
      FILE_divergence_date.push_back(divergence_date[i]);
      FILE_relgen.push_back(generation);
      if (abundances[i] > 0)
      {
        FILE_tage.push_back(generation-born[i]);
      }
      else
      {
        FILE_tage.push_back(-1);
      }
      FILE_page.push_back(page[i]) ;
    }
  }

  long tempindex = parents[topindex];
  // topindex is the index of the root node %^&

  double diffgen = born[topindex];
  parents[topindex] = -1;
  generation -= diffgen;
  for (unsigned int i = 0 ; i < parents.size() ; ++i)
  {
    if (abundances[i] >= 0)
    {
      born[i] -= diffgen;
      if (died[i] >= 0)
      {
        died[i] -= diffgen;
      }
      maxabtime1[i] -= diffgen;
      maxabtime2[i] -= diffgen;
    }
  }

  born[topindex] = 0.0;

  while (tempindex >= 0)
  {
    long   newtempindex = parents[tempindex];

    abundances[tempindex] = -2;
    max_abundance[tempindex] = -2;
    maxabtime1[tempindex] = -2.0;
    maxabtime2[tempindex] = -2.0;
    parents[tempindex] = -3;
    born[tempindex] = -2.0;
    died[tempindex] = -2.0;
    selection_level[tempindex] = -2000000;
    num_mut[tempindex]  = -2.0; // %^&
    num_descend_spec[tempindex] = -2;

    if (end_recycler < static_cast<int>(recycler.size()))
    {
      recycler[end_recycler] = tempindex;
    }
    else
    {
      recycler.push_back(tempindex);
    }
    end_recycler ++;

    tempindex = newtempindex;

  }

}

void NTsim::set_seed(long seedin)
{
  if (!seeded)
  {
    the_seed = seedin;
    seeded = true;
    NR.set_seed(seedin);
  }
}

void NTsim::setup(long seedin , long J_M_in , double mu_in , double s_in)
{
  // Prepare for a new simulation with the given parameters
  // note there is also a system restore funciton lower down in this file
  total_clear();
  set_seed(seedin);

  // Set simulation parmaeters correcty first
  J_M = J_M_in;
  mu = mu_in;
  s = s_in;

  // setup metacommunity
  for (long i = 0 ; i < J_M ; ++i)
  {
    metacommunity.push_back(upto);
  }
  upto ++;
  richness = 1;

  // setup burnin counting variables
  for (long i = 0 ; i < J_M ; ++i)
  {
    init_ancestor.push_back(i);
    num_descend.push_back(1);
  }
  num_remain = J_M;

  // species level data std::vectors to be set with initial conditions
  abundances.push_back(J_M);
  parents.push_back(-1);
  born.push_back(0.0);
  died.push_back(-1.0);
  selection_level.push_back(0);
  num_mut.push_back(0); // %^&

  max_abundance.push_back(J_M);
  maxabtime1.push_back(0.0);
  maxabtime2.push_back(0.0);

  num_descend_spec.push_back(1);

  // now setup the output file vars

  FILE_final_newick.clear();

  FILE_generation.clear();
  FILE_iname.clear();
  FILE_iname2.clear();
  FILE_index.clear();
  FILE_abundances.clear();
  FILE_parents.clear();
  FILE_born.clear();
  FILE_died.clear();
  FILE_selection_level.clear();
  FILE_num_mut.clear();
  FILE_max_abundance.clear();
  FILE_maxabtime1.clear();
  FILE_maxabtime2.clear();
  FILE_num_descend_spec.clear();
  FILE_divergence_date.clear();
  FILE_relgen.clear();
  FILE_tage.clear();
  FILE_page.clear();

  FILE_generation_var.clear();
  FILE_n.clear();
  FILE_n_richness.clear();
  FILE_num_subspec.clear();
  FILE_tot_abund.clear();
  FILE_meanc.clear();
  FILE_varc.clear();
  FILE_rangec.clear();
  FILE_meanc_w.clear();
  FILE_varc_w.clear();
  FILE_rangec_w.clear();

}

void NTsim::sim_all(long seedin , long J_M_in , double mu_in , double s_in , double simtime)
{
  // initialise timing points
  time_t start , end;
  time(&start);

  bool timeout = false;
  long steps = 0;

  setup(seedin , J_M_in , mu_in , s_in);

  while (((num_remain > 1)||(abundances[0] > 0))&&(!timeout))
  {
    sim_step();
    steps ++;

    if (steps%10000 == 0) //£$% - removed 00
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
    gen_burned_in = true_generation; // inside the if statement beacuse if we have run out of time then this is not the correct gen_burned_in

    while ((!timeout)&&(iter <= 1000000)) // checks in time and also that not too much HDD space has been used
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

      if (true_generation >= (gen_burned_in*double(simnum)/10.0 +4.0)) // 4 times the expected minimum burnin (hence the +4 here - up from 3) the /10 increases the frequency of outputting
      {

        if (simnum == 1)
        {
          std::cout << "simulation fully burned in\n";
        }
        output();
        simnum ++;
      }
    }
  }
  // do a final output just in case the system needs to be suspended here
  output();
}

void NTsim::sim_restore(char* fnamein,double simtime)
{
  std::cout << "restoring NNIM from " << fnamein << "\n";
  std::ifstream in;
  in.open(fnamein);

  // scalars

  in >> J_M ;
  in >> mu ;
  in >> s ;
  in >> num_suspends ;
  in >> num_file_outputs ;
  in >> generation ;
  in >> true_generation ;
  in >> gen_burned_in ;
  in >> simnum ;
  in >> max_selection_level ;
  in >> richness ;
  in >> upto ;
  in >> end_recycler ;
  in >> num_remain ;
  in >> iter ;
  in >> iter2 ;
  in >> seeded ;
  in >> the_seed ;

  // std::vectors

  long temp_size;
  long temp_long;
  double temp_double;

  in >> temp_size;
  recycler.clear();
  for (long i = 0 ; i < temp_size ; ++i)
  {
    in >> temp_long;
    recycler.push_back(temp_long);
  }

  in >> temp_size;
  init_ancestor.clear();
  for (long i = 0 ; i < temp_size ; ++i)
  {
    in >> temp_long;
    init_ancestor.push_back(temp_long);
  }

  in >> temp_size;
  num_descend.clear();
  for (long i = 0 ; i < temp_size ; ++i)
  {
    in >> temp_long;
    num_descend.push_back(temp_long);
  }

  in >> temp_size;
  metacommunity.clear();
  for (long i = 0 ; i < temp_size ; ++i)
  {
    in >> temp_long;
    metacommunity.push_back(temp_long);
  }

  in >> temp_size;
  abundances.clear();
  for (long i = 0 ; i < temp_size ; ++i)
  {
    in >> temp_long;
    abundances.push_back(temp_long);
  }

  in >> temp_size;
  parents.clear();
  for (long i = 0 ; i < temp_size ; ++i)
  {
    in >> temp_long;
    parents.push_back(temp_long);
  }

  in >> temp_size;
  born.clear();
  for (long i = 0 ; i < temp_size ; ++i)
  {
    in >> temp_double;
    born.push_back(temp_double);
  }

  in >> temp_size;
  died.clear();
  for (long i = 0 ; i < temp_size ; ++i)
  {
    in >> temp_double;
    died.push_back(temp_double);
  }

  in >> temp_size;
  selection_level.clear();
  for (long i = 0 ; i < temp_size ; ++i)
  {
    in >> temp_long;
    selection_level.push_back(temp_long);
  }

  in >> temp_size;
  num_mut.clear();
  for (long i = 0 ; i < temp_size ; ++i)
  {
    in >> temp_long;
    num_mut.push_back(temp_long);
  }

  in >> temp_size;
  max_abundance.clear();
  for (long i = 0 ; i < temp_size ; ++i)
  {
    in >> temp_long;
    max_abundance.push_back(temp_long);
  }

  in >> temp_size;
  maxabtime1.clear();
  for (long i = 0 ; i < temp_size ; ++i)
  {
    in >> temp_double;
    maxabtime1.push_back(temp_double);
  }

  in >> temp_size;
  maxabtime2.clear();
  for (long i = 0 ; i < temp_size ; ++i)
  {
    in >> temp_double;
    maxabtime2.push_back(temp_double);
  }

  in >> temp_size;
  num_descend_spec.clear();
  for (long i = 0 ; i < temp_size ; ++i)
  {
    in >> temp_long;
    num_descend_spec.push_back(temp_long);
  }

  std::vector<long> NR_suspend;

  in >> temp_size;
  NR_suspend.clear();
  for (long i = 0 ; i < temp_size ; ++i)
  {
    in >> temp_long;
    NR_suspend.push_back(temp_long);
  }

  NR.resume(NR_suspend);

  in.close();

  std::cout << "system restored \n";

  // initialise timing points
  time_t start , end;
  time(&start);

  bool timeout = false;
  long steps = 0;

  if (gen_burned_in < 0)
  {
    while (((num_remain > 1)||(abundances[0] > 0))&&(!timeout))
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
      gen_burned_in = true_generation;
    }
  }

  std::cout << "simulation part burned in\n";

  if (!timeout)
  {

    while ((!timeout)&&(iter <= 1000000)) // checks in time and also that not too much HDD space has been used
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

      if (true_generation >= (gen_burned_in*double(simnum)/10.0 +4.0)) // 4 times the expected minimum burnin (hence the +3 here) the /10 increase the frequency of outputting
      {

        if (simnum == 1)
        {
          std::cout << "simulation fully burned in\n";
        }

        output();
        simnum ++;

      }
    }
  }

}

void NTsim::sim_step()
{

  // increment generation counter
  generation += 2.0/J_M;
  true_generation += 2.0/J_M;

  // choose individual to die (and be replaced)
  long chosen2;
  chosen2 = NR.i0(J_M-1);

  // decide if speciation happened
  bool speciation = false;
  if (NR.d01() <= mu)
  {
    speciation = true;
  }

  // chooose replacement individual
  long chosen;
  bool repeat_selection = true;
  while(repeat_selection)
  {
    // keep going until selection criteria are satisfied
    chosen = NR.i0(J_M-1);
    if (chosen != chosen2)
    {
      if (s > 0)
      {
        // deals with selection
        long temp_selection = selection_level[metacommunity[chosen]];
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
          if (NR.d01()<=(pow((1.0+s),-1.0*double(max_selection_level - temp_selection -1))))
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
          if (NR.d01()<=(pow((1.0+s),-1.0*double(max_selection_level - temp_selection))))
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
  if (gen_burned_in<0)
  {
    num_descend[init_ancestor[chosen]]++;
    num_descend[init_ancestor[chosen2]]--;
    if (num_descend[init_ancestor[chosen2]] == 0)
    {
      num_remain --;
    }
    init_ancestor[chosen2] = init_ancestor[chosen];
  }

  // now update the reaminder of the variables based on
  // speciation, chosen and chosen2

  // update main data

  abundances[metacommunity[chosen2]]--;
  if (abundances[metacommunity[chosen2]] == 0)
  {

    died[metacommunity[chosen2]] = generation;
    richness --;

    long tempindex = metacommunity[chosen2]; // the species index of the now dead species
    bool keeplooping = true;
    while(keeplooping)
    {
      num_descend_spec[tempindex] --;

      long newtempindex = parents[tempindex]; // move up to see if more can be recycled now

      if (num_descend_spec[tempindex] == 0)
      {
        // recycle

        // set flags just to be sure it's recycled (good for redundancy in case of bugs so data that slipps through is obviously wrong)
        abundances[tempindex] = -2;
        max_abundance[tempindex] = -2;
        maxabtime1[tempindex] = -2.0;
        maxabtime2[tempindex] = -2.0;
        parents[tempindex] = -3;
        born[tempindex] = -2.0;
        died[tempindex] = -2.0;
        selection_level[tempindex] = -2000000;
        num_mut[tempindex] = -2;
        num_descend_spec[tempindex] = -2;

        if (end_recycler < static_cast<int>(recycler.size()))
        {
          recycler[end_recycler] = tempindex;
        }
        else
        {
          recycler.push_back(tempindex);
        }
        end_recycler ++;

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
    if (end_recycler == 0)
    {
      // no species numbers available to recycle
      metacommunity[chosen2] = upto;
      abundances.push_back(1);
      max_abundance.push_back(1);
      maxabtime1.push_back(generation);
      maxabtime2.push_back(generation);

      parents.push_back(metacommunity[chosen]);
      born.push_back(generation);
      died.push_back(-1.0);

      num_mut.push_back(1);

      if (NR.d01() < 0.5)
      {
        selection_level.push_back(selection_level[metacommunity[chosen]]-1);
      }
      else
      {
        selection_level.push_back(selection_level[metacommunity[chosen]]+1);
        if ((selection_level[metacommunity[chosen]]+1)>max_selection_level)
        {
          max_selection_level = (selection_level[metacommunity[chosen]]+1);
        }
      }

      num_descend_spec.push_back(1);

      upto ++;
    }
    else
    {
      // species are available to recycle
      long tempindex = recycler[end_recycler-1];

      metacommunity[chosen2] = tempindex;
      abundances[tempindex] = 1;
      max_abundance[tempindex] = 1;
      maxabtime1[tempindex] = generation;
      maxabtime2[tempindex] = generation;

      parents[tempindex] = metacommunity[chosen];
      born[tempindex] = generation;
      died[tempindex] = -1.0;

      num_mut[tempindex] = 1;

      if (NR.d01() < 0.5)
      {
        selection_level[tempindex] = (selection_level[metacommunity[chosen]]-1);
      }
      else
      {
        selection_level[tempindex] = (selection_level[metacommunity[chosen]]+1);

        if ((selection_level[tempindex])>max_selection_level)
        {
          max_selection_level = selection_level[tempindex];
        }

      }

      num_descend_spec[tempindex] = 1;

      end_recycler --;
    }

    long tempindex = metacommunity[chosen];
    while(tempindex >= 0)
    {
      num_descend_spec[tempindex] ++;
      tempindex = parents[tempindex]; // move up to see if more can be recycled now
    }

    richness ++;

  }
  else
  {
    // no speciation
    abundances[metacommunity[chosen]]++;
    if ((abundances[metacommunity[chosen]]) == max_abundance[metacommunity[chosen]])
    {
      maxabtime2[metacommunity[chosen]] = generation;
    }
    if ((abundances[metacommunity[chosen]]) > max_abundance[metacommunity[chosen]])
    {
      maxabtime2[metacommunity[chosen]] = generation;
      maxabtime1[metacommunity[chosen]] = generation;
      max_abundance[metacommunity[chosen]] = abundances[metacommunity[chosen]];
    }

    metacommunity[chosen2] = metacommunity[chosen];
  }
}

void NTsim::total_clear()
{
  generation = 0.0;
  true_generation = 0.0;
  gen_burned_in = -1.0;
  simnum = 1;
  max_selection_level = 0;
  richness = 0;

  upto = 0;
  recycler.clear();
  end_recycler = 0;

  init_ancestor.clear();
  num_descend.clear();
  num_remain = -1;

  metacommunity.clear();

  abundances.clear();
  parents.clear();
  born.clear();
  died.clear();
  selection_level.clear();
  num_mut.clear();

  max_abundance.clear();
  maxabtime1.clear();
  maxabtime2.clear();

  num_descend_spec.clear();
  // random number generator, simulation parameters are all set elsewhere
}



void NTsim::write_files()
{
  std::cout << "writing files \n";

  std::ofstream out;
  out.precision(10);

  if (FILE_generation.size() > 0) // check that there is anything to output
  {
    num_file_outputs ++;

    // now setup the output file vars
    // species level file
    char filename_ab[100];
    sprintf (filename_ab, "Data_%i_FV_Out_%i.csv", int(the_seed),int(num_file_outputs)); // edited here to make csv &&&&&&&
    out.open(filename_ab);
    // output just the header now
    out <<  "generation,iname,iname2,index,abundances,parents,born,died,selection_level,num_mut,max_abundance,maxabtime1,maxabtime2,num_descend_spec,divergence_date,relgen,tage,page\n";

    for (unsigned long i = 0 ; i < FILE_generation.size() ; ++i)
    {
      out <<  FILE_generation[i] <<",";
      out <<  FILE_iname[i] <<",";
      out <<  FILE_iname2[i] <<",";
      out <<  FILE_index[i] <<",";
      out <<  FILE_abundances[i] <<",";
      out <<  FILE_parents[i] <<",";
      out <<  FILE_born[i] <<",";
      out <<  FILE_died[i] <<",";
      out <<  FILE_selection_level[i] <<",";
      out <<  FILE_num_mut[i] <<",";
      out <<  FILE_max_abundance[i] <<",";
      out <<  FILE_maxabtime1[i] <<",";
      out <<  FILE_maxabtime2[i] <<",";
      out <<  FILE_num_descend_spec[i] <<",";
      out <<  FILE_divergence_date[i] <<",";
      out <<  FILE_relgen[i] <<",";
      out <<  FILE_tage[i] <<",";
      out <<  FILE_page[i] <<",\n";
    }
    out.close();

    // tree file
    char filename_newick[100];
    sprintf (filename_newick, "Data_%i_FV_Newick_%i.txt", int(the_seed),int(num_file_outputs));
    out.open(filename_newick);
    for (unsigned long i = 0 ; i < FILE_final_newick.size() ; ++i)
    {
      out << FILE_final_newick[i] << ";\n";
    }
    out.close();

    // variance file
    sprintf (filename_ab, "Data_%i_FV_FitVar_%i.csv", int(the_seed),int(num_file_outputs)); // edited here to make csv &&&&&&&
    out.open(filename_ab);
    // output just the header now
    out <<  "generation_var,n,n_richness,num_subspec,tot_abund,meanc,varc,rangec,meanc_w,varc_w,rangec_w\n";

    for (unsigned long i = 0 ; i < FILE_generation_var.size() ; ++i)
    {
      out <<  FILE_generation_var[i] <<",";
      out <<  FILE_n[i] <<",";
      out <<  FILE_n_richness[i] <<",";
      out <<  FILE_num_subspec[i] <<",";
      out <<  FILE_tot_abund[i] <<",";
      out <<  FILE_meanc[i] <<",";
      out <<  FILE_varc[i] <<",";
      out <<  FILE_rangec[i] <<",";
      out <<  FILE_meanc_w[i] <<",";
      out <<  FILE_varc_w[i] <<",";
      out <<  FILE_rangec_w[i] <<",\n";

    }
    out.close();
  }

  if (num_suspends == 0)
  {
    char filename_vars[100];
    sprintf (filename_vars, "Data_%i_FV_Vars.txt", int(the_seed));
    out.open(filename_vars);
    out << the_seed << " \n";
    out << J_M << " \n";
    out << mu << " \n";
    out << s << " \n";
    out.close();
  }

  num_suspends ++;

  char filename_vars[100];
  sprintf (filename_vars, "Data_%i_FV_Res_%i.txt", int(the_seed) ,int(num_suspends));
  out.open(filename_vars);
  // strore the simulation state so we can continue another time if necessary

  out.precision(15);

  // scalars

  out << J_M << " ";
  out << mu << " ";
  out << s << " ";
  out << num_suspends << " ";
  out << num_file_outputs << " ";
  out << generation << " ";
  out << true_generation << " ";
  out << gen_burned_in << " ";
  out << simnum << " ";
  out << max_selection_level << " ";
  out << richness << " ";
  out << upto << " ";
  out << end_recycler << " ";
  out << num_remain << " ";
  out << iter << " ";
  out << iter2 << " ";
  out << seeded << " ";
  out << the_seed << " ";

  // std::vectors

  out << recycler.size() << " ";
  for (unsigned long i = 0 ; i < recycler.size() ; ++i)
  {
    out << recycler[i] << " ";
  }

  out << init_ancestor.size() << " ";
  for (unsigned long i = 0 ; i < init_ancestor.size() ; ++i)
  {
    out << init_ancestor[i] << " ";
  }

  out << num_descend.size() << " ";
  for (unsigned long i = 0 ; i < num_descend.size() ; ++i)
  {
    out << num_descend[i] << " ";
  }

  out << metacommunity.size() << " ";
  for (unsigned long i = 0 ; i < metacommunity.size() ; ++i)
  {
    out << metacommunity[i] << " ";
  }

  out << abundances.size() << " ";
  for (unsigned long i = 0 ; i < abundances.size() ; ++i)
  {
    out << abundances[i] << " ";
  }

  out << parents.size() << " ";
  for (unsigned long i = 0 ; i < parents.size() ; ++i)
  {
    out << parents[i] << " ";
  }

  out << born.size() << " ";
  for (unsigned long i = 0 ; i < born.size() ; ++i)
  {
    out << born[i] << " ";
  }

  out << died.size() << " ";
  for (unsigned long i = 0 ; i < died.size() ; ++i)
  {
    out << died[i] << " ";
  }

  out << selection_level.size() << " ";
  for (unsigned long i = 0 ; i < selection_level.size() ; ++i)
  {
    out << selection_level[i] << " ";
  }

  out << num_mut.size() << " ";
  for (unsigned long i = 0 ; i < num_mut.size() ; ++i)
  {
    out << num_mut[i] << " ";
  }

  out << max_abundance.size() << " ";
  for (unsigned long i = 0 ; i < max_abundance.size() ; ++i)
  {
    out << max_abundance[i] << " ";
  }

  out << maxabtime1.size() << " ";
  for (unsigned long i = 0 ; i < maxabtime1.size() ; ++i)
  {
    out << maxabtime1[i] << " ";
  }

  out << maxabtime2.size() << " ";
  for (unsigned long i = 0 ; i < maxabtime2.size() ; ++i)
  {
    out << maxabtime2[i] << " ";
  }

  out << num_descend_spec.size() << " ";
  for (unsigned long i = 0 ; i < num_descend_spec.size() ; ++i)
  {
    out << num_descend_spec[i] << " ";
  }

  std::vector<long> NR_suspend = NR.suspend();
  out << NR_suspend.size() << " ";
  for (unsigned long i = 0 ; i < NR_suspend.size() ; ++i)
  {
    out << NR_suspend[i] << " ";
  }

  out.close();
}
