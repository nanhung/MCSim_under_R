/* random.h

   Copyright (c) 1992-2008 Free Software Foundation, Inc.

   This file is part of GNU MCSim.

   GNU MCSim is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 3
   of the License, or (at your option) any later version.

   GNU MCSim is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with GNU MCSim; if not, see <http://www.gnu.org/licenses/>

   -- Revisions -----
     Logfile:  %F%
    Revision:  %I%
        Date:  %G%
     Modtime:  %U%
      Author:  @a
   -- SCCS  ---------

   Header for random number generator.  See random.c for extensive
   documentation.

   Gives prototypes for random functions.

*/

#ifndef RANDOM_H_DEFINED

#include <math.h>

/* ----------------------------------------------------------------------------
   Constants  */

#define SEED_MIN     1.0
#define SEED_MAX     2147483646.0
#define SEED_DEFAULT 314159265.3589793
#define PI           3.1415926535897932384626433
#define INV_SQRT_2PI 0.398942280401433
#define SQRT_2       1.4142135623731


/* ----------------------------------------------------------------------------
   Prototypes  */

/* Initialize the random generators, optional but recommended */
void InitRandom (double dSeed, int bWarmUp);

/* Two random generators */
/*     one that shuffles its output, */
double RandomShuffle (void);

/*     and one that doesn't */
double Randoms (void);


/* Several types of random variates */

double BetaRandom (double alpha, double beta, double a, double b);
double BinomialBetaRandom (double Expectation, double alpha, double beta);
double BinomialRandom (double p, long n);
double CauchyRandom (double dScale);
double Chi2Random (double dof);
double ExpRandom (double beta);
double InvGGammaRandom (double alpha, double beta);
double TruncInvGGammaRandom (double alpha, double beta, double a, double b);
double GammaRandom (double alpha);
double GGammaRandom (double alpha, double beta);
double LogNormalRandom (double dMean, double dStdDev);
double GenLogNormalRandom (double dMean, double dStdDevNorm, 
			   double dStdDevLogNorm);
double StudentTRandom (double dof, double dMean, double dStdDev);
double LogUniformRandom (double a, double b);
double NormalRandom (double dMean, double dStdDev);
double PiecewiseRandom (double min, double a, double b, double max);
double PiecewiseVariate (long cDim, double rg_x[], double rg_pdf[],
                         double rg_Cdf[], int iOrder, double *pVal_pdf);
long   PoissonRandom (double mu);
double TruncLogNormalRandom (double dMean, double dStdDev, double a, double b);
double TruncNormalRandom (double dMean, double dStdDev, double a, double b);
double UniformRandom (double a, double b);
void   Multinomial (long n, int dim, double *p, double *x);
void   WishartRandom (long n, long p, double *t, double *w, double *work);


/* ----------------------------------------------------------------------------
   Utility functions */

BOOL   and (BOOL A, BOOL B);
double CDFNormal (double z);
double InterpolateX (double rgX[], double rgY[], long lLower, double dY);
double erfc (double x);
double GetSeed (void);
double lnDFNormal (double x, double mu, double sd);
double lnGamma (double x);
double lnDFBeta (double x, double alpha, double beta, double min, double max);
void   CalcCumulative (long cDim, double *rg_x, double *rg_pdf,
                       double *rg_Cdf, int  iOrder);
void   SetSeed (double dSeed);
 
#define RANDOM_H_DEFINED
#endif

/* End */

