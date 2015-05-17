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

#ifndef NRRAND_H
#define NRRAND_H

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

#include <cmath>
#include <iostream>
#include <vector>


/************************************************************
 RANDOM NUMBER GENERATOR WRAPPER OBJECT
 FROM NUMERICAL RECIPES IN C
************************************************************/

class NRrand
{

private:
  long m_idum;
  long m_j;
  long m_k;
  long m_idum2;
  long m_iy;
  long m_iv[NTAB];
  bool m_seeded;

  double m_temp;

public:

  NRrand();

  void set_seed(long seed)
  {
    if (!m_seeded)
    {
      m_idum2 = 123456789;
      m_iy = 0;
      m_idum = seed;
      if (m_idum < 1) m_idum=1;
      //Be sure to prevent idum = 0.
      m_idum2=(m_idum);
      for (m_j=NTAB+7;m_j>=0;m_j--)
      {
        //Load the shuffle table (after 8 warm-ups).
        m_k=(m_idum)/IQ1;
        m_idum=IA1*(m_idum-m_k*IQ1)-m_k*IR1;
        if (m_idum < 0) m_idum += IM1;
        if (m_j < NTAB) m_iv[m_j] = m_idum;
      }
      m_iy=m_iv[0];
      m_seeded = true;
    }
  }

  double d01()
  {
    m_k=(m_idum)/IQ1;
    //Start here when not initializing.
    m_idum=IA1*(m_idum-m_k*IQ1)-m_k*IR1;
    //Compute idum=(IA1*idum) % IM1 without overflows by Schrage's method.
    if (m_idum < 0) m_idum += IM1;
    m_k=m_idum2/IQ2;
    m_idum2=IA2*(m_idum2-m_k*IQ2)-m_k*IR2;
    //Compute idum2=(IA2*idum) % IM2 likewise.
    if (m_idum2 < 0) m_idum2 += IM2;
    m_j=m_iy/NDIV;
    //Will be in the range 0..NTAB-1.
    m_iy=m_iv[m_j]-m_idum2;
    //Here idum is shuffled, idum and idum2 are combined to generate output.
    m_iv[m_j] = m_idum;
    if (m_iy < 1) m_iy += IMM1;
    if ((m_temp=AM*m_iy) > RNMX)
    {
      //cout << "random call = " << "RNMAX" << "\n";
      return RNMX; //Because users don't expect endpoint values.
    }
    else
    {
      return m_temp;
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
    return (long(std::floor(d01()*(max+1))));
  }

  std::vector<long> suspend()
  {
    std::vector<long> toret;
    toret.clear();
    if (m_seeded)
    {
      toret.push_back(m_idum);
      toret.push_back(m_j);
      toret.push_back(m_k);
      toret.push_back(m_idum2);
      toret.push_back(m_iy);

      for (int i = 0 ; i < NTAB ; ++i)
      {
        toret.push_back(m_iv[i]);
      }
      // we assume the system is seeded if we are suspending and resuming
    }
    return (toret);
  }

  void resume(std::vector<long> datain)
  {
    m_seeded = true;
    if (datain.size() == (5+NTAB))
    {
      m_idum = datain[0];
      m_j = datain[1];
      m_k = datain[2];
      m_idum2 = datain[3];
      m_iy = datain[4];
      for (int i = 5 ; i < NTAB+5 ; ++i)
      {
        m_iv[i-5] = datain[i];
      }
    }
    else
    {
      std::cout << "ERROR, cannot resume random number generator \n";
    }
  }

};


#endif // NRRAND_H
