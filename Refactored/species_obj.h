#ifndef SPECIES_OBJ_H
#define SPECIES_OBJ_H

#include <iostream>
#include <vector>

/************************************************************
 GOOD SPECIES OBJECT
************************************************************/
class species_obj {

private:
  //E.g. consider a good and incipition species, where
  // both trait and abundances are` kept together
  //
  //   |
  //   |
  // +-+-+
  // |   |
  // +   +
  // A   B
  //
  //Each species has a collection of traits
  std::vector<long> c;
  std::vector<long> abund; //Number of individuals


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

  std::vector<long> get_c()
  {
    return(c);
  }

  std::vector<long> get_abund()
  {
    return(abund);
  }

  long get_tot_abund()
  {
    long toret = 0;
    if (abund.size() > 0)
    {

      for (long i = 0 ; i < static_cast<int>(abund.size()) ; ++i)
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
    std::vector<long> tempabund = x.get_abund();
    std::vector<long> tempc = x.get_c();
    if (tempabund.size() != tempc.size())
    {
      std::cout << "error in species object when combining two";
    }
    for (int i = 0 ; i < tempabund.size() ; ++i)
    {
      c.push_back(tempc[i]);
      abund.push_back(tempabund[i]);
    }
  }

  long get_subspec()
  {
    long num_spec_total = 0;
    for (long i = 0 ; i < abund.size() ; ++i)
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
      for (long i = 0 ; i < c.size() ; ++i)
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
      for (long i = 0 ; i < c.size() ; ++i)
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
      for (long i = 1 ; i < c.size() ; ++i)
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
      for (long i = 1 ; i < c.size() ; ++i)
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

#endif // SPECIES_OBJ_H
