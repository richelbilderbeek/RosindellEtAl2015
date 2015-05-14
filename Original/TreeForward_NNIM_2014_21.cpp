// Code for the article
//   Rosindell, James, Luke J. Harmon, and Rampal S. Etienne.
//   "Unifying ecology and macroevolution with individual‐based theory."
//   Ecology letters 18.5 (2015): 472-482.

// Version 1.2 Jan 2014
// Author James Rosindell james@rosindell.org

// NNIM = nearly neutral individual mutations

// positive or negative mutations at individual level,
// speciation is the result of a number of mutations together in build up

// parameters:

// JM (simulations size of metacommunity - non spatial)
// u probability of mutation ±
// s size of mutation
// n number of mutations necessary to count as a new species (accounted for here for fitness variations as a spectrum)

bool debug = false;

/************************************************************
 INCLUDES
 ************************************************************/

// standard inludes
# include <stdio.h>
# include <fstream>
# include <vector>
# include <iostream>
# include <string>
# include <math.h>
# include <time.h>
# include <ctime>
# include <sstream>

using namespace std;

/************************************************************
 RANDOM NUMBER GENERATOR WRAPPER OBJECT
 FROM NUMERICAL RECIPES IN C
 ************************************************************/

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 5277
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-8
#define RNMX (1.0-EPS)

# include <stdio.h>
# include <iostream>
# include <string>
# include <math.h>

# include <vector>
# include <iostream>
# include <fstream>

class NRrand
{
    
private:
	long idum;
	long j;
	long k;
	long idum2;
	long iy;
	long iv[NTAB];
	bool seeded;
    
    double temp;
    
public:
    
	NRrand()
	{
		seeded = false;
	}
    
	void set_seed(long seed)
	{
		if (!seeded)
		{
			idum2 = 123456789;
			iy = 0;
			idum = seed;
			if (idum < 1) idum=1;
			//Be sure to prevent idum = 0.
			idum2=(idum);
			for (j=NTAB+7;j>=0;j--)
			{
				//Load the shuffle table (after 8 warm-ups).
				k=(idum)/IQ1;
				idum=IA1*(idum-k*IQ1)-k*IR1;
				if (idum < 0) idum += IM1;
				if (j < NTAB) iv[j] = idum;
			}
			iy=iv[0];
			seeded = true;
		}
	}
    
	double d01()
	{
		k=(idum)/IQ1;
		//Start here when not initializing.
		idum=IA1*(idum-k*IQ1)-k*IR1;
		//Compute idum=(IA1*idum) % IM1 without overflows by Schrage's method.
		if (idum < 0) idum += IM1;
		k=idum2/IQ2;
		idum2=IA2*(idum2-k*IQ2)-k*IR2;
		//Compute idum2=(IA2*idum) % IM2 likewise.
		if (idum2 < 0) idum2 += IM2;
		j=iy/NDIV;
		//Will be in the range 0..NTAB-1.
		iy=iv[j]-idum2;
		//Here idum is shuffled, idum and idum2 are combined to generate output.
		iv[j] = idum;
		if (iy < 1) iy += IMM1;
		if ((temp=AM*iy) > RNMX)
		{
			//cout << "random call = " << "RNMAX" << "\n";
			return RNMX; //Because users don't expect endpoint values.
		}
		else
		{
			return temp;
		}
	}
    
    bool event(double probin)
	{
		if (probin < 0.000001)
		{
			if (d01() <= 0.000001)
			{
				return (event(probin * 1000000.0));
			}
			else
			{
				return false;
			}
		}
		else
		{
			if (probin > 0.999999)
			{
				return (!(event(1.0-probin)));
			}
			else
			{
				return (d01() <= probin);
			}
		}
	}
    
    long i0(long max)
	// integer between 0 and max inclusive
	{
		return (long(floor(d01()*(max+1))));
	}
    
    vector<long> suspend()
    {
        vector<long> toret;
        toret.clear();
        if (seeded)
        {
            toret.push_back(idum);
            toret.push_back(j);
            toret.push_back(k);
            toret.push_back(idum2);
            toret.push_back(iy);
            
            for (int i = 0 ; i < NTAB ; i ++)
            {
                toret.push_back(iv[i]);
            }
            // we assume the system is seeded if we are suspending and resuming
        }
        return (toret);
    }
    
    void resume(vector<long> datain)
    {
        seeded = true;
        if (datain.size() == (5+NTAB))
        {
            idum = datain[0];
            j = datain[1];
            k = datain[2];
            idum2 = datain[3];
            iy = datain[4];
            for (int i = 5 ; i < NTAB+5 ; i ++)
            {
                iv[i-5] = datain[i];
            }
        }
        else
        {
            cout << "ERROR, cannot resume random number generator \n";
        }
    }
    
};

/************************************************************
 GOOD SPECIES OBJECT
 ************************************************************/

class species_obj {
    
    
private:
    //E.g. consider a good and incipition species, where
    // both trait and abundances are kept together
    //
    //   |
    //   |
    // +-+-+
    // |   |
    // +   +
    // A   B
    //
    //Each species has a collection of traits
    vector<long> c;
    vector<long> abund; //Number of individuals
    
    
public:
    
    species_obj()
    {
        clear();
    }
    
    void clear()
    {
        c.clear();
        abund.clear();
    }
    
    void setup(long abund_in , long c_in)
    {
        clear();
        abund.push_back(abund_in);
        c.push_back(c_in);
    }
    
    vector<long> get_c()
    {
        return(c);
    }
    
    vector<long> get_abund()
    {
        return(abund);
    }
    
    long get_tot_abund()
    {
        long toret = 0;
        if (abund.size() > 0)
        {
            
            for (long i = 0 ; i < abund.size() ; i ++)
            {
                toret += abund[i];
            }
            return toret;
        }
        else
        {
            return 0;
        }
    }
    
    void add(species_obj x)
    {
        vector<long> tempabund = x.get_abund();
        vector<long> tempc = x.get_c();
        if (tempabund.size() != tempc.size())
        {
            cout << "error in species object when combining two";
        }
        for (int i = 0 ; i < tempabund.size() ; i ++)
        {
            c.push_back(tempc[i]);
            abund.push_back(tempabund[i]);
        }
    }
    
    long get_subspec()
    {
        long num_spec_total = 0;
        for (long i = 0 ; i < abund.size() ; i ++)
        {
            if (abund[i] > 0)
            {
                num_spec_total++;
            }
        }
        return(num_spec_total);
    }
    
    double get_meanc()
    {
        if ((c.size() > 0)&&(get_tot_abund() > 0))
        {
            double toret = 0;
            for (long i = 0 ; i < c.size() ; i ++)
            {
                toret += double(c[i]*abund[i]);
            }
            return (toret/double(get_tot_abund()));
        }
        else
        {
            return 0;
        }
        
    }
    
    double get_varc()
    {
        if ((c.size() > 1)&&(get_tot_abund() > 1))
        {
            double tot = 0;
            double tot2 = 0;
            double numr = double(get_tot_abund());
            for (long i = 0 ; i < c.size() ; i ++)
            {
                tot += double(c[i]*abund[i]);
                tot2 += double(c[i])*double(c[i]*abund[i]);
            }
            return ((numr/(numr-1)) * ((tot2/numr)-((tot/numr)*(tot/numr))));
        }
        else
        {
            return 0;
        }
    }
    
    long get_maxc()
    {
        if (c.size() > 0)
        {
            long toret = c[0];
            for (long i = 1 ; i < c.size() ; i ++)
            {
                if (toret < c[i])
                {
                    toret = c[i];
                }
            }
            return toret;
        }
        else
        {
            return 0;
        }
    }
    
    long get_minc()
    {
        if (c.size() > 0)
        {
            long toret = c[0];
            for (long i = 1 ; i < c.size() ; i ++)
            {
                if (toret > c[i])
                {
                    toret = c[i];
                }
            }
            return toret;
        }
        else
        {
            return 0;
        }
    }
    
};

/************************************************************
 NEUTRAL SIMULATION OBJECT
 ************************************************************/

class NTsim {
	
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
    vector<long> recycler;
    // a list of indices in the per-species data that can now be used because the species are extinct
    long end_recycler;
    // the end of the recycler vector (in case it is required to contract) so it's equal to 0 when recyler is empty
    
	// *** burn in checker and per individual data - becomes irrelevant after first output ***
	
	vector<long> init_ancestor;
	// the ancestor of this individual in the initial state - index referres to currently living individuals in the metacommunity
	vector<long> num_descend;
	//the number of current descendents of this individual - index referrs to the index of individuals initialised in the metacommunity
  	long num_remain;
	// the number of values in num_descend that are non zero (when this decreases to one we're burned in)
    
    // *** per individual data that remains relevant after burn in
    
    vector<long> metacommunity;
	// metacommunity - contains the index in per-species data of the species of each individual (for ease of simulation) - the index of this vector is the position in the metacommunity
    
    // *** per-species simulation data ***
    
    vector<long> abundances;
	// the abundnance of each species
	vector<long> parents;
	// the parent species for each species
	vector<double> born;
	// the date of species birth for each species
	vector<double> died;
	// the date of species death for each species (-1 means not dead)
	vector<long> selection_level;
	// the selection level of each species (larger numbers mean greater advantage)
    vector<long> num_mut;
    // number of mutations along the branch that connects this to the parent node, these could be positive or negative mutations
    
    vector<long> max_abundance;
    // the maximum abundance ever experienced by this species
    vector<double> maxabtime1;
    // the earliest date when that maximum abundance was hit
    vector<double> maxabtime2;
    // the latest date when that maximum abundance was hit
    
    vector<long> num_descend_spec;
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
    vector<string> FILE_final_newick;
    
    // out (per species)
    
    // header of the filal file will read
    // generation,iname,iname2,index,abundances,parents,born,died,selection_level,num_mut,max_abundance,maxabtime1,maxabtime2,num_descend_spec,divergence_date,relgen,tage,page
    
    vector<double> FILE_generation;
    vector<string> FILE_iname;
    vector<long> FILE_iname2;
    vector<long> FILE_index;
    vector<long> FILE_abundances;
    vector<long> FILE_parents;
    vector<double> FILE_born;
    vector<double> FILE_died;
    vector<long> FILE_selection_level;
    vector<long> FILE_num_mut;
    vector<long> FILE_max_abundance;
    vector<double> FILE_maxabtime1;
    vector<double> FILE_maxabtime2;
    vector<long> FILE_num_descend_spec;
    vector<double> FILE_divergence_date;
    vector<double> FILE_relgen;
    vector<double> FILE_tage;
    vector<double> FILE_page;
    
    // res and vars are based on overall sim parameters at end of sim and don't need to be saved.
    
    // fitvar shows the fitness variation overall and within each species for a variety of differnet values of n
    // header of the final file will read: generation_var, meanc_all , varc_all , n , num_subspec , meanc , varc , maxc , minc
    
    vector<double> FILE_generation_var; // will line up with the variance output rows for frequency of output here
    vector<long> FILE_n; // value of n that applies for the rest of the readings (the first three above will probably repeat a lot)
    vector<long> FILE_n_richness; // the species richness overall at this value of n
    vector<double> FILE_num_subspec; // number of sub species - average over all species
    vector<double> FILE_tot_abund; // total abundance of all sub species - average over all species
    vector<double> FILE_meanc; // mean value of c - (within species) average over all species
    vector<double> FILE_varc; // variange in c - (within species) average over all species
  	vector<double> FILE_rangec; // range in value of c - (within species) average over all species
    vector<double> FILE_meanc_w; // mean value of c - (within species) weighted average over all species
    vector<double> FILE_varc_w; // variange in c - (within species) weighted average over all species
  	vector<double> FILE_rangec_w; // range in value of c - (within species) weighted average over all species
    
    
public:
    
	NTsim() // initialiser
	{
		total_clear();
		seeded = false;
        iter = 0;
        num_suspends = 0;
        num_file_outputs = 0;
	}
    
    void set_seed(long seedin) // set seed in random number generator
	{
		if (!seeded)
		{
			the_seed = seedin;
			seeded = true;
			NR.set_seed(seedin);
		}
	}
    
	void total_clear() // clears parameter values
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
    
	void setup(long seedin , long J_M_in , double mu_in , double s_in)
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
        for (long i = 0 ; i < J_M ; i ++)
		{
			metacommunity.push_back(upto);
		}
        upto ++;
        richness = 1;
        
        // setup burnin counting variables
        for (long i = 0 ; i < J_M ; i ++)
		{
			init_ancestor.push_back(i);
			num_descend.push_back(1);
		}
		num_remain = J_M;
        
        // species level data vectors to be set with initial conditions
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
    
 	void sim_step()
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
                    
                    if (end_recycler < recycler.size())
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
    
    // a set of tools for converting numbers to letters for species names in phylogenies
    string alpha_format(long x)
	{
    return 'A' + x;
    /*
		string alpha_res;
		
		switch (x) {
			case 0:
				alpha_res = "A";
				break;
			case 1:
				alpha_res = "B";
				break;
			case 2:
				alpha_res = "C";
				break;
			case 3:
				alpha_res = "D";
				break;
			case 4:
				alpha_res = "E";
				break;
			case 5:
				alpha_res = "F";
				break;
			case 6:
				alpha_res = "G";
				break;
			case 7:
				alpha_res = "H";
				break;
			case 8:
				alpha_res = "I";
				break;
			case 9:
				alpha_res = "J";
				break;
			case 10:
				alpha_res = "K";
				break;
			case 11:
				alpha_res = "L";
				break;
			case 12:
				alpha_res = "M";
				break;
			case 13:
				alpha_res = "N";
				break;
			case 14:
				alpha_res = "O";
				break;
			case 15:
				alpha_res = "P";
				break;
			case 16:
				alpha_res = "Q";
				break;
			case 17:
				alpha_res = "R";
				break;
			case 18:
				alpha_res = "S";
				break;
			case 19:
				alpha_res = "T";
				break;
			case 20:
				alpha_res = "U";
				break;
			case 21:
				alpha_res = "V";
				break;
			case 22:
				alpha_res = "W";
				break;
			case 23:
				alpha_res = "X";
				break;
			case 24:
				alpha_res = "Y";
				break;
			case 25:
				alpha_res = "Z";
				break;
			default:
				alpha_res = "_ERROR_";
		}
		return alpha_res;
    */
	}
	
	string alpha_format(long x_in, int num_digits)
	{
		string alpha_res;
		alpha_res = "";
		long x;
		x = x_in;
		if (pow(26.0,double(num_digits)) >= x)
		{
			for (int i = 0 ; i < num_digits ; i ++)
			{
				long remainder = (x % 26);
				alpha_res = alpha_format(remainder) + alpha_res;
				x -= remainder;
				x = x/26;
			}
		}
		else
		{
			alpha_res = "_ERROR_";
		}
		return alpha_res;
	}
	
	string alpha_format_auto(long x_in)
	{
		long upto = 25;
		long num_digits = 1;
		long x = x_in;
		while (upto < x_in)
		{
			x -= long(floor(pow(26.0,num_digits)));
			num_digits ++;
			upto += long(floor(pow(26.0,num_digits)));
		}
		return (alpha_format(x,num_digits));
	}
    
    // this function simulates a set of parameters for a set time and gets as many results as it can
    // all these reults are dropped in a common file named based on the seed and parameter compbination
    // the time is measured in minutes
    void sim_all(long seedin , long J_M_in , double mu_in , double s_in , double simtime)
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
                    // £$%
                    
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
                    
                    // £$%
                   // cout << difftime(end,start) << " , " << true_generation << " - burning in (phase 2)\n";
                    
                    //if (true_generation >= 154000) debug = true;
                    //if (debug) cout << "debug is on\n";
                    
                }
                
                if (true_generation >= (gen_burned_in*double(simnum)/10.0 +4.0)) // 4 times the expected minimum burnin (hence the +4 here - up from 3) the /10 increases the frequency of outputting
                {
                    
                    if (simnum == 1)
                    {
                        cout << "simulation fully burned in\n";
                    }
                    output();
                    simnum ++;
                }
            }
        }
        // do a final output just in case the system needs to be suspended here
        output();
	}
    
	// this function outputs the data including calculating the fitness variation data
    // it also sorts it so that no longer relevant parts are thrown out.
	void output()
	{
        // cout << "outputting \n";
        
        bool output_phy = true; // output the phylogeny?
        bool output_var = true; // output the variance?
        
        // first of all work out if each line should appear in the final output
        // extinct species that were not a most recent common ancestor between two species are not outputted
        // the date of divergence from the (grand)parent species is noted in another column of data now
        vector<bool> keep_data;
        vector<double> divergence_date;
        vector<long> intermediate_parents;
        vector<long> num_mut2 = num_mut;
        keep_data.clear();
        divergence_date.clear();
        intermediate_parents.clear();
        for (int i = 0 ; i < parents.size() ; i ++)
        {
            keep_data.push_back(true);
            divergence_date.push_back(-1.0);
            intermediate_parents.push_back(parents[i]);
        }
        for (int i = 0 ; i < parents.size() ; i ++)
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
        
        vector<long> newi;
        newi.clear();
        
        for (int i = 0 ; i < parents.size() ; i ++)
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
        vector<long> newi2; // counts from 1 up in each simulation on its own
        newi2.clear();
        
        for (int i = 0 ; i < parents.size() ; i ++)
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
        // each node stores a num descendents that are accounted for so far and two vectors
        // the vectors are of child nodes and store their newick strings and dates
        
        string final_newick; // the final answer
        vector<long> num_descend_done;
        vector< vector<string> > newick_parts;
        vector< vector<double> > join_dates;
        vector< vector<double> > last_split_dates; // the last spit of a lineage - needed so that joindates-lastsplit generally gives us a branch length
        
        vector<double> page; // the phylogenetic age of the species for outputting.
        
        num_descend_done.clear();
        newick_parts.clear();
        join_dates.clear();
        last_split_dates.clear();
        
        page.clear();
        
        
        if (output_phy)
        {
            
            for (int i = 0 ; i < parents.size() ; i ++)
            {
                string temps = "";
                vector<string> tempvs;
                tempvs.clear();
                vector<double> tempvd;
                tempvd.clear();
                if (keep_data[i])
                {
                    num_descend_done.push_back(0); // negative number indicates completed, positive number indicates number of descendents added to the newick_parts and join_dates vectors
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
            
            for (int i = 0 ; i < parents.size() ; i ++)
            {
                if (abundances[i] > 0) // automatically excludes any lineages that might not be kept
                {
                    // we're at the end
                    ostringstream ss;
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
                for (long i = 0 ; i < parents.size() ; i ++)
                {
                    // check if we're done
                    if ((num_descend_done[i] >= 0)&&(keep_data[i])) // num_descend done goes negative when we're done
                    {
                        keep_looping = true;
                        
                        if (num_descend_done[i] == num_descend_spec[i])
                        {
                            
                            // can sort out, otherwise can do nothing
                            // first order the two vectors
                            if ((newick_parts[i]).size() > 1) // only if ordering needed
                            {
                                bool unsorted = true;
                                while(unsorted)
                                {
                                    unsorted = false;
                                    for (int j = 1 ; j < (newick_parts[i]).size() ; j ++)
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
                                            
                                            string temps = newick_parts[i][j-1];
                                            newick_parts[i][j-1] = newick_parts[i][j];
                                            newick_parts[i][j] = temps;
                                        }
                                    }
                                }
                            }
                            
                            // this is the full string
                            ostringstream ss2;
                            ss2.precision(10);
                            double thislastsplit = -1.0;
                            
                            if ((newick_parts[i]).size() > 1) // only if ordering needed
                            {
                                // this is the interior node name
                                string int_node;
                                ostringstream ss;
                                ss.precision(10);
                                ss << alpha_format_auto(newi[i]);
                                ss << "["; // %^&
                                ss << newi2[i]; // %^&
                                ss << "]"; // %^&
                                int_node = ss.str();
                                
                                ss2 << "(";
                                
                                if((newick_parts[i]).size() >2)
                                {
                                    for (int j = 0 ; j < ((newick_parts[i]).size()-2) ; j ++)
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
            
            
            for (long i = 0 ; i < parents.size() ; i ++)
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
            
            // vector<double> FILE_generation_var; // will line up with the variance output rows for frequency of output here
            // vector<long> FILE_n; // value of n that applies for the rest of the readings (the first three above will probably repeat a lot)
            // vector<long> FILE_n_richness; // the species richness overall at this value of n
            // vector<long> FILE_num_subspec; // number of sub species - for this good species only
            //vector<long> FILE_tot_abund; // total abundance of all sub species - for this good species only
            // vector<double> FILE_meanc; // mean value of c - for this good species only
            // vector<double> FILE_varc; // variange in c - for this good species only
            // vector<long> FILE_maxc; // max value in c - for this good species only
            // vector<long> FILE_minc; // min value in c - for this good species only
            
            
            /*
             new data:
             
             vector<double> FILE_generation_var; // will line up with the variance output rows for frequency of output here
             vector<long> FILE_n; // value of n that applies for the rest of the readings (the first three above will probably repeat a lot)
             vector<long> FILE_n_richness; // the species richness overall at this value of n
             vector<long> FILE_num_subspec; // number of sub species - average over all species
             vector<long> FILE_tot_abund; // total abundance of all sub species - average over all species
             vector<double> FILE_meanc; // mean value of c - (within species) average over all species
             vector<double> FILE_varc; // variange in c - (within species) average over all species
             vector<long> FILE_rangec; // range in value of c - (within species) average over all species
             vector<double> FILE_meanc_w; // mean value of c - (within species) weighted average over all species
             vector<double> FILE_varc_w; // variange in c - (within species) weighted average over all species
             vector<long> FILE_rangec_w; // range in value of c - (within species) weighted average over all species

             
             */
            
            
            // DATA WE HAVE AVAILABLE TO US
            
            // vector<long> abundances;
            // the abundnance of each species
            // vector<long> parents;
            // the parent species for each species
            // vector<long> selection_level;
            // the selection level of each species (larger numbers mean greater advantage)
            // vector<long> num_mut;
            // number of mutations along the branch that connects this to the parent node, these could be positive or negative mutations
            
            // THE ALGORITHM
            
            // to do these calculations I'm going to use the species_obj to record details of good species
            // I will prune the tree and maintain a vector of good species to one side
            // I will repeat the process for different values of n including infinity
            
            // to perform this I'll store a number of useful vectors that will line up with the vectors of all sub species
            
            // choose a value of n to do this for
            //long this_n = -1; // loop over more values of n
            
            
            
            for (int this_n = -1 ; this_n < 70 ; )
            {
                
                // initialise a min_num_mutations (not including mutations up to parent) to living species vector for each position so that when working up the tree lumping can be done easily
                vector<long> min_num_mutations;
                min_num_mutations.clear();
                
                // initialise a min_num_mutations to living species vector for each position including the possibility that the fewest mutations is through the parent
                vector<long> min_num_dir;
                min_num_dir.clear();
                
                // initialise min_num_which which shows us which path gives the route to a living species with the minimum - this takes on the value of the node with the minimum number of mutations path to a living species.
                // Note that a draw doesn't matter either can be chosen as they'll be treated just the same later anyway - either both lumped or both split.
                // Note that the route could be through the parent
                //vector<long> min_num_which; // I don't think we actually need this after all commented out everywhere now
                //min_num_which.clear();
                
                // initialise a num_child vector (down to all direct chidren nodes) for all nodes as a useful tool
                vector<long> num_child;
                num_child.clear();
                
                // initilise a num_done vector (tells us how many of the children we have dealt with so far with pruning - a completely pruned branch should have this equal to num_child)
                // after a completely pruned branch has had it's details dealt with and pushed up to its parent the num_done tag changes to -1
                vector<long> num_done;
                num_done.clear();
                
                // initialise a vector pruned_spec of  species objects that gets passed up the tree and amalgamated as appropriate - this is stored on the parent of the species
                vector<species_obj> pruned_spec;
                pruned_spec.clear();
                
                // here is a vector that will store the final result
                vector<species_obj> final_spec;
                final_spec.clear();
                // species object vector contains the abundance and c values of all sub species and can produce the necessary outputs.
                
                
                
                // now we fill all these vectors with their initial conditions.
                for (long i = 0 ; i < abundances.size() ; i ++)
                {
                    species_obj newspec;
                    newspec.setup(abundances[i],selection_level[i]);
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
                for (long i = 0 ; i < abundances.size() ; i ++)
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
                    for (long i = 0 ; i < abundances.size() ; i ++)
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
                    for (long i = 0 ; i < abundances.size() ; i ++)
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
                
                
                // now do the actual pruning step after the vectors have been set up
                // we're going to keep looping until every node is dealt with
                looping = true;
                while(looping)
                {
                    looping = false;
                    for (long i = 0 ; i < abundances.size() ; i ++)
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
                                        
                                        pruned_spec[parents[i]].add(pruned_spec[i]);
                                        num_done[i] = -1;
                                        num_done[parents[i]] ++;
                                        pruned_spec[i].clear();
                                    }
                                    else
                                    {
                                        // splitting
                                        if ((pruned_spec[i]).get_tot_abund() > 0)
                                        {
                                            final_spec.push_back(pruned_spec[i]);
                                        }
                                        num_done[i] = -1;
                                        num_done[parents[i]] ++;
                                        pruned_spec[i].clear();
                                    }
                                }
                            }
                        }
                        // else nothing to do for this loop around
                    }
                }
                
                
                
                // could be that the final parent needs to be added still
                for (long i = 0 ; i < abundances.size() ; i ++)
                {
                    if ((parents[i] == -1)&&((pruned_spec[i]).get_subspec()>0))
                    {
                        final_spec.push_back(pruned_spec[i]);
                    }
                }
                
                
                long num_spec_total = 0;
                long num_ind_total = 0;
                for (long i = 0 ; i < final_spec.size() ; i ++)
                {
                    if ((final_spec[i].get_tot_abund())>0)
                    {
                        num_spec_total ++;
                        num_ind_total += (final_spec[i].get_tot_abund());
                    }
                }
                
                /*
                
                 vector<double> FILE_generation_var; // will line up with the variance output rows for frequency of output here
                 vector<long> FILE_n; // value of n that applies for the rest of the readings (the first three above will probably repeat a lot)
                 vector<long> FILE_n_richness; // the species richness overall at this value of n
                 vector<long> FILE_num_subspec; // number of sub species - average over all species
                 vector<long> FILE_tot_abund; // total abundance of all sub species - average over all species
                 vector<double> FILE_meanc; // mean value of c - (within species) average over all species
                 vector<double> FILE_varc; // variange in c - (within species) average over all species
                 vector<long> FILE_rangec; // range in value of c - (within species) average over all species
                 vector<double> FILE_meanc_w; // mean value of c - (within species) weighted average over all species
                 vector<double> FILE_varc_w; // variange in c - (within species) weighted average over all species
                 vector<long> FILE_rangec_w; // range in value of c - (within species) weighted average over all species
                 
                //*/
                
                
                double temp_num_subspec = 0;
                double temp_tot_abund = 0;
                
                double temp_meanc = 0;
                double temp_varc = 0;
                double temp_rangec = 0;
                
                double temp_meanc_w = 0;
                double temp_varc_w = 0;
                double temp_rangec_w = 0;
                
                for (long i = 0 ; i < final_spec.size() ; i ++)
                {
                    if ((final_spec[i].get_tot_abund())>0)
                    {
                        
                         temp_num_subspec += (double(final_spec[i].get_subspec()));
                         temp_tot_abund += (double(final_spec[i].get_tot_abund()));
                        
                         temp_meanc += (double(final_spec[i].get_meanc()));
                         temp_varc += (double(final_spec[i].get_varc()));
                         temp_rangec += (double(final_spec[i].get_maxc()-final_spec[i].get_minc()));
                        
                         temp_meanc_w += (double(final_spec[i].get_tot_abund())*double(final_spec[i].get_meanc()));
                         temp_varc_w += (double(final_spec[i].get_tot_abund())*double(final_spec[i].get_varc()));
                         temp_rangec_w += (double(final_spec[i].get_tot_abund())*double(final_spec[i].get_maxc()-final_spec[i].get_minc()));

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
        
        for (int i = 0 ; i < parents.size() ; i ++)
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
        for (int i = 0 ; i < parents.size() ; i ++)
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
            
            if (end_recycler < recycler.size())
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
    
    void sim_restore(char* fnamein,double simtime)
    {
        cout << "restoring NNIM from " << fnamein << "\n";
        ifstream in;
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
        
        // vectors
        
        long temp_size;
        long temp_long;
        double temp_double;
        
        in >> temp_size;
        recycler.clear();
        for (long i = 0 ; i < temp_size ; i ++)
        {
            in >> temp_long;
            recycler.push_back(temp_long);
        }
        
        in >> temp_size;
        init_ancestor.clear();
        for (long i = 0 ; i < temp_size ; i ++)
        {
            in >> temp_long;
            init_ancestor.push_back(temp_long);
        }
        
        in >> temp_size;
        num_descend.clear();
        for (long i = 0 ; i < temp_size ; i ++)
        {
            in >> temp_long;
            num_descend.push_back(temp_long);
        }
        
        in >> temp_size;
        metacommunity.clear();
        for (long i = 0 ; i < temp_size ; i ++)
        {
            in >> temp_long;
            metacommunity.push_back(temp_long);
        }
        
        in >> temp_size;
        abundances.clear();
        for (long i = 0 ; i < temp_size ; i ++)
        {
            in >> temp_long;
            abundances.push_back(temp_long);
        }
        
        in >> temp_size;
        parents.clear();
        for (long i = 0 ; i < temp_size ; i ++)
        {
            in >> temp_long;
            parents.push_back(temp_long);
        }
        
        in >> temp_size;
        born.clear();
        for (long i = 0 ; i < temp_size ; i ++)
        {
            in >> temp_double;
            born.push_back(temp_double);
        }
        
        in >> temp_size;
        died.clear();
        for (long i = 0 ; i < temp_size ; i ++)
        {
            in >> temp_double;
            died.push_back(temp_double);
        }
        
        in >> temp_size;
        selection_level.clear();
        for (long i = 0 ; i < temp_size ; i ++)
        {
            in >> temp_long;
            selection_level.push_back(temp_long);
        }
        
        in >> temp_size;
        num_mut.clear();
        for (long i = 0 ; i < temp_size ; i ++)
        {
            in >> temp_long;
            num_mut.push_back(temp_long);
        }
        
        in >> temp_size;
        max_abundance.clear();
        for (long i = 0 ; i < temp_size ; i ++)
        {
            in >> temp_long;
            max_abundance.push_back(temp_long);
        }
        
        in >> temp_size;
        maxabtime1.clear();
        for (long i = 0 ; i < temp_size ; i ++)
        {
            in >> temp_double;
            maxabtime1.push_back(temp_double);
        }
        
        in >> temp_size;
        maxabtime2.clear();
        for (long i = 0 ; i < temp_size ; i ++)
        {
            in >> temp_double;
            maxabtime2.push_back(temp_double);
        }
        
        in >> temp_size;
        num_descend_spec.clear();
        for (long i = 0 ; i < temp_size ; i ++)
        {
            in >> temp_long;
            num_descend_spec.push_back(temp_long);
        }
        
        vector<long> NR_suspend;
        
        in >> temp_size;
        NR_suspend.clear();
        for (long i = 0 ; i < temp_size ; i ++)
        {
            in >> temp_long;
            NR_suspend.push_back(temp_long);
        }
        
        NR.resume(NR_suspend);
        
        in.close();
        
        cout << "system restored \n";
        
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
        
        cout << "simulation part burned in\n";
        
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
                        cout << "simulation fully burned in\n";
                    }
                    
                    output();
                    simnum ++;
                    
                }
            }
        }
        
    }
    
    void write_files()
    {
        cout << "writing files \n";
        
        ofstream out;
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
            
            for (long i = 0 ; i < FILE_generation.size() ; i ++)
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
            for (long i = 0 ; i < FILE_final_newick.size() ; i ++)
            {
                out << FILE_final_newick[i] << ";\n";
            }
            out.close();
            
            // variance file
            sprintf (filename_ab, "Data_%i_FV_FitVar_%i.csv", int(the_seed),int(num_file_outputs)); // edited here to make csv &&&&&&&
            out.open(filename_ab);
            // output just the header now
            out <<  "generation_var,n,n_richness,num_subspec,tot_abund,meanc,varc,rangec,meanc_w,varc_w,rangec_w\n";
            
            for (long i = 0 ; i < FILE_generation_var.size() ; i ++)
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
        
        // vectors
        
        out << recycler.size() << " ";
        for (long i = 0 ; i < recycler.size() ; i ++)
        {
            out << recycler[i] << " ";
        }
        
        out << init_ancestor.size() << " ";
        for (long i = 0 ; i < init_ancestor.size() ; i ++)
        {
            out << init_ancestor[i] << " ";
        }
        
        out << num_descend.size() << " ";
        for (long i = 0 ; i < num_descend.size() ; i ++)
        {
            out << num_descend[i] << " ";
        }
        
        out << metacommunity.size() << " ";
        for (long i = 0 ; i < metacommunity.size() ; i ++)
        {
            out << metacommunity[i] << " ";
        }
        
        out << abundances.size() << " ";
        for (long i = 0 ; i < abundances.size() ; i ++)
        {
            out << abundances[i] << " ";
        }
        
        out << parents.size() << " ";
        for (long i = 0 ; i < parents.size() ; i ++)
        {
            out << parents[i] << " ";
        }
        
        out << born.size() << " ";
        for (long i = 0 ; i < born.size() ; i ++)
        {
            out << born[i] << " ";
        }
        
        out << died.size() << " ";
        for (long i = 0 ; i < died.size() ; i ++)
        {
            out << died[i] << " ";
        }
        
        out << selection_level.size() << " ";
        for (long i = 0 ; i < selection_level.size() ; i ++)
        {
            out << selection_level[i] << " ";
        }
        
        out << num_mut.size() << " ";
        for (long i = 0 ; i < num_mut.size() ; i ++)
        {
            out << num_mut[i] << " ";
        }
        
        out << max_abundance.size() << " ";
        for (long i = 0 ; i < max_abundance.size() ; i ++)
        {
            out << max_abundance[i] << " ";
        }
        
        out << maxabtime1.size() << " ";
        for (long i = 0 ; i < maxabtime1.size() ; i ++)
        {
            out << maxabtime1[i] << " ";
        }
        
        out << maxabtime2.size() << " ";
        for (long i = 0 ; i < maxabtime2.size() ; i ++)
        {
            out << maxabtime2[i] << " ";
        }
        
        out << num_descend_spec.size() << " ";
        for (long i = 0 ; i < num_descend_spec.size() ; i ++)
        {
            out << num_descend_spec[i] << " ";
        }
        
        vector<long> NR_suspend = NR.suspend();
        out << NR_suspend.size() << " ";
        for (long i = 0 ; i < NR_suspend.size() ; i ++)
        {
            out << NR_suspend[i] << " ";
        }
        
        out.close();
        
    }
};

/************************************************************
 MAIN
 ************************************************************/

int charconvertor(char charin)
{
	switch (charin)
	{
		case '0': return (0);
		case '1': return (1);
		case '2': return (2);
		case '3': return (3);
		case '4': return (4);
		case '5': return (5);
		case '6': return (6);
		case '7': return (7);
		case '8': return (8);
		case '9': return (9);
		default: return (-1);
	}
}

int jobconvertor(char* argin)
{
	int maxind = 0;
	while (charconvertor(argin[maxind]) != -1)
	{
		maxind ++;
	}
	int jobtoret = 0;
	int pow10 = 1;
	for (int i = maxind-1 ; i >=0 ; i --)
	{
		jobtoret += (pow10*charconvertor(argin[i]));
		pow10 = pow10*10;
	}
	return jobtoret;
}

int main(int argc, char* argv[])
{
    
    //NTsim test;
    //cout << "simulating M with fitness increase or decrease ... ";
    //test.sim_all(2 , 10000 , 0.05 , 0.001 , 60);
    //test.write_files();
    
    
    
    //*
     // first get the command line arg data in
     
     long task_iter = -1; // the task number
     double task_max = -1; // the max time (minutes)
     
     task_iter = jobconvertor(argv[1]) -1;
     task_max = jobconvertor(argv[2]);
     
     cout << task_iter << " = task \n";
     cout << task_max << " = max time \n";
     
     vector<long> JM_vect;
     vector<double> nu_vect;
     vector<double> S_vect;
     
     JM_vect.clear();
     nu_vect.clear();
     S_vect.clear();
     
     JM_vect.push_back(10000);
     JM_vect.push_back(100000);
     JM_vect.push_back(1000000);
    
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
    
     vector<long> JM_list;
     vector<double> nu_list;
     vector<double> S_list;
     
     for (int i = 0 ; i < JM_vect.size() ; i ++)
     {
     for (int j = 0 ; j < nu_vect.size() ; j ++)
     {
     for (int k = 0 ; k < S_vect.size() ; k ++)
     {
     JM_list.push_back(JM_vect[i]);
     nu_list.push_back(nu_vect[j]);
     S_list.push_back(S_vect[k]);
     }
     }
     }
     
     NTsim test;
     
     if (argc == 4)
     {
     // restore
     test.sim_restore(argv[3],60*task_max);
     }
     else
     {
     cout << "simulating M with fitness increase or decrease ... ";
     cout << JM_list[task_iter] << " , " << nu_list[task_iter] << " , " <<  S_list[task_iter] << " \n";
     
     test.sim_all(task_iter+10000 , JM_list[task_iter] , nu_list[task_iter] , S_list[task_iter] , 60*task_max);
     }
     test.write_files();
     // seed, JM, mu, S, simtime (60*seconds = mins)
     // long, long, double, double, double
     
     //*/
    
}

// notes: copy files to $WORK where there is 150 gb of storage available
// restrict output to 1M lines (0.22 GB storage and 1.1 hrs per sorting) - done in main simulation routine loop
// 0.22 Mb per 1K lines and 3.95 seconds sorting per 1K lines

