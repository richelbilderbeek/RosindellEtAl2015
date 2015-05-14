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
    return (long(std::floor(d01()*(max+1))));
  }

  std::vector<long> suspend()
  {
    std::vector<long> toret;
    toret.clear();
    if (seeded)
    {
      toret.push_back(idum);
      toret.push_back(j);
      toret.push_back(k);
      toret.push_back(idum2);
      toret.push_back(iy);

      for (int i = 0 ; i < NTAB ; ++i)
      {
        toret.push_back(iv[i]);
      }
      // we assume the system is seeded if we are suspending and resuming
    }
    return (toret);
  }

  void resume(std::vector<long> datain)
  {
    seeded = true;
    if (datain.size() == (5+NTAB))
    {
      idum = datain[0];
      j = datain[1];
      k = datain[2];
      idum2 = datain[3];
      iy = datain[4];
      for (int i = 5 ; i < NTAB+5 ; ++i)
      {
        iv[i-5] = datain[i];
      }
    }
    else
    {
      std::cout << "ERROR, cannot resume random number generator \n";
    }
  }

};


#endif // NRRAND_H
