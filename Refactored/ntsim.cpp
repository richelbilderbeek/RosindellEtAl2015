#include "ntsim.h"

NTsim::NTsim()
{
  total_clear();
  seeded = false;
  iter = 0;
  num_suspends = 0;
  num_file_outputs = 0;
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

    for (long i = 0 ; i < FILE_generation.size() ; ++i)
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
    for (long i = 0 ; i < FILE_final_newick.size() ; ++i)
    {
      out << FILE_final_newick[i] << ";\n";
    }
    out.close();

    // variance file
    sprintf (filename_ab, "Data_%i_FV_FitVar_%i.csv", int(the_seed),int(num_file_outputs)); // edited here to make csv &&&&&&&
    out.open(filename_ab);
    // output just the header now
    out <<  "generation_var,n,n_richness,num_subspec,tot_abund,meanc,varc,rangec,meanc_w,varc_w,rangec_w\n";

    for (long i = 0 ; i < FILE_generation_var.size() ; ++i)
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
  for (long i = 0 ; i < recycler.size() ; ++i)
  {
    out << recycler[i] << " ";
  }

  out << init_ancestor.size() << " ";
  for (long i = 0 ; i < init_ancestor.size() ; ++i)
  {
    out << init_ancestor[i] << " ";
  }

  out << num_descend.size() << " ";
  for (long i = 0 ; i < num_descend.size() ; ++i)
  {
    out << num_descend[i] << " ";
  }

  out << metacommunity.size() << " ";
  for (long i = 0 ; i < metacommunity.size() ; ++i)
  {
    out << metacommunity[i] << " ";
  }

  out << abundances.size() << " ";
  for (long i = 0 ; i < abundances.size() ; ++i)
  {
    out << abundances[i] << " ";
  }

  out << parents.size() << " ";
  for (long i = 0 ; i < parents.size() ; ++i)
  {
    out << parents[i] << " ";
  }

  out << born.size() << " ";
  for (long i = 0 ; i < born.size() ; ++i)
  {
    out << born[i] << " ";
  }

  out << died.size() << " ";
  for (long i = 0 ; i < died.size() ; ++i)
  {
    out << died[i] << " ";
  }

  out << selection_level.size() << " ";
  for (long i = 0 ; i < selection_level.size() ; ++i)
  {
    out << selection_level[i] << " ";
  }

  out << num_mut.size() << " ";
  for (long i = 0 ; i < num_mut.size() ; ++i)
  {
    out << num_mut[i] << " ";
  }

  out << max_abundance.size() << " ";
  for (long i = 0 ; i < max_abundance.size() ; ++i)
  {
    out << max_abundance[i] << " ";
  }

  out << maxabtime1.size() << " ";
  for (long i = 0 ; i < maxabtime1.size() ; ++i)
  {
    out << maxabtime1[i] << " ";
  }

  out << maxabtime2.size() << " ";
  for (long i = 0 ; i < maxabtime2.size() ; ++i)
  {
    out << maxabtime2[i] << " ";
  }

  out << num_descend_spec.size() << " ";
  for (long i = 0 ; i < num_descend_spec.size() ; ++i)
  {
    out << num_descend_spec[i] << " ";
  }

  std::vector<long> NR_suspend = NR.suspend();
  out << NR_suspend.size() << " ";
  for (long i = 0 ; i < NR_suspend.size() ; ++i)
  {
    out << NR_suspend[i] << " ";
  }

  out.close();
}
