#ifndef SPECIES_OBJ_H
#define SPECIES_OBJ_H

#include <iostream>
#include <vector>

/************************************************************
 GOOD SPECIES OBJECT
************************************************************/
struct Species
{
  Species();

  void Add(const Species& x);



  void Clear();

  long GetSumAbundances();

  void Setup(long abund_in, long c_in);

  private:


  const std::vector<long>& GetTraits() const noexcept { return m_traits; }
  const std::vector<long>& GetAbundances() const noexcept { return m_abundances; }

  public:

  long GetNumberOfSpecies();


  double GetMeanTrait()
  {
    if ((m_traits.size() > 0)&&(GetSumAbundances() > 0))
    {
      double toret = 0;
      for (long i = 0 ; i < m_traits.size() ; ++i)
      {
        toret += double(m_traits[i]*m_abundances[i]);
      }
      return (toret/double(GetSumAbundances()));
    }
    else
    {
      return 0;
    }

  }

  double get_varc()
  {
    if ((m_traits.size() > 1)&&(GetSumAbundances() > 1))
    {
      double tot = 0;
      double tot2 = 0;
      double numr = double(GetSumAbundances());
      for (long i = 0 ; i < m_traits.size() ; ++i)
      {
        tot += double(m_traits[i]*m_abundances[i]);
        tot2 += double(m_traits[i])*double(m_traits[i]*m_abundances[i]);
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
    if (m_traits.size() > 0)
    {
      long toret = m_traits[0];
      for (long i = 1 ; i < m_traits.size() ; ++i)
      {
        if (toret < m_traits[i])
        {
          toret = m_traits[i];
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
    if (m_traits.size() > 0)
    {
      long toret = m_traits[0];
      for (long i = 1 ; i < m_traits.size() ; ++i)
      {
        if (toret > m_traits[i])
        {
          toret = m_traits[i];
        }
      }
      return toret;
    }
    else
    {
      return 0;
    }
  }

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
  std::vector<long> m_traits;
  std::vector<long> m_abundances; //Number of individuals
};

#endif // SPECIES_OBJ_H
