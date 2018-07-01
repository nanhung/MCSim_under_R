/* mh.h

   Originally written by Frederic Bois
   
   Copyright (c) 1996-2008 Free Software Foundation, Inc.

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

   Header file for mh.c
*/

#ifndef MH_H_DEFINED

/* ----------------------------------------------------------------------------
   Inclusions
*/

#include <float.h>  /* Floating point limits */
#include <stdarg.h>

#include "sim.h"
#include "config.h"

#define UPDATE_BASE 20  /* Change MH criterion */

/* ----------------------------------------------------------------------------
   Global/External variables
*/
#define NTEMP 5


/* ----------------------------------------------------------------------------
   Definitions
*/

#define MAX_ARGS              25

#define NULL_SUPPORT          -1.0E+100
#define MISSING_VALUE         (-DBL_MAX)
#define INPUT_MISSING_VALUE   -1

/* ----------------------------------------------------------------------------
   Typedefs
*/


/* ----------------------------------------------------------------------------
   Prototypes
*/

void   CalculateTotals (PLEVEL plevel, char **args);
void   CheckForFixed (PLEVEL plevel, char **args);
void   CheckPrintStatements (PLEVEL plevel, char **args);
void   CloneLikes (PLEVEL plevel, char **args);
void   CloneLikesL (PVOID pData, PVOID pUser1, PVOID pUser2, PVOID pUser3);
void   CloneMCVars (PLEVEL plevel, char **args);
void   CloneMCVarsL (PVOID pData, PVOID pUser1, PVOID pUser2, PVOID pUser3);
void   CloseMarkovFiles (PANALYSIS panal);
void   ConvertLists (PLEVEL plevel, char **args);
double dInvTemper( long indexT);
void   DoMarkov (PANALYSIS panal);
void   FindLikeParents (PLEVEL plevel, char **args);
void   FindMCDependents(PLEVEL plevel, char **args);
void   FindMCParents (PLEVEL plevel, char **args);
void   GetNumberOfMCVars (PLEVEL plevel, char **args);
void   InitArrays (long lDim, long nSubjs, double ***pdSum, 
                   double ****prgdSumProd);
void   InitMCVars (PLEVEL plevel, char **args);
void   ListToPMCArray (PANALYSIS panal, PLIST plist,
                       long *nMCVars, PMCVAR **rgpMCVars);
void   ListToPMCArrayL (PVOID pData, PVOID pUser1, PVOID pUser2, PVOID pUser3);
void   ListToPVArray (PANALYSIS panal, PLIST plist,
                      long *nFixedVars, PVARMOD **rgpFixedVars);
void   ListToPVArrayL (PVOID pData, PVOID pUser1, PVOID pUser2, PVOID pUser3);
double LnDensity (PMCVAR pMCVar, PANALYSIS panal);
double LnLike (PMCVAR pMCVar, PANALYSIS panal);
double LnLikeData (PLEVEL plevel, PANALYSIS panal);
double MaxMCVar (PMCVAR pMCVar);
double MinMCVar (PMCVAR pMCVar);
void   OpenMarkovFiles (PANALYSIS panal);
void   PrintAllExpts (PLEVEL plevel, PANALYSIS panal, PFILE pOutFile);
void   PrintDeps (PLEVEL plevel, char **args);
int    PrintExpt (PLEVEL plevel, char **args);
void   ReadData (PLEVEL plevel, char **args);
void   ReadRestart (FILE *pfileRestart, long nThetas, double *pdTheta, 
                    double *pdSum, double **prgdSumProd, long *pIter);
void   ReadRestartTemper (FILE *pfileRestart, long nThetas, 
                          int nInvTemperatures, double *pdTheta, 
                          double *pdSum, double **prgdSumProd,
                          long *pIter, long *pindexT, double *pdlnPi);
int    RestoreLikelihoods (PLEVEL plevel, char **args);
int    RunAllExpts(PANALYSIS panal, PDOUBLE pdLnData);
int    RunExpt (PLEVEL plevel, char **args);
double SampleTheta (PMCVAR pMCVar);
long   SampleTemperature (PGIBBSDATA pgd, double dLnPrior, double dLnData);
void   SampleThetas (PLEVEL plevel, char **args);
void   SampleThetasTempered (PLEVEL plevel, char **args);
void   SampleThetaVector (PLEVEL pLevel, PANALYSIS panal, long nThetas,
                          double *pdTheta, double *pdSum, double **prgdSumProd,
                          long iter, long nUpdateAt, long nTotal, 
                          PDOUBLE pdLnPrior, PDOUBLE pdLnData);
int    SaveLikelihoods (PLEVEL plevel, char **args);
void   SetFixedVars (PLEVEL plevel);
void   SetKernel (PLEVEL plevel, char **args);
void   ReadKernel (PLEVEL plevel, char **args);
void   WriteKernel (PLEVEL plevel, char **args);
void   SetModelVars (PLEVEL plevel);
int    SetMCVars (PLEVEL plevel, char **args);
int    SumAllExpts (PLEVEL plevel, char **args);
void   SetPointers (PLEVEL plevel, char **args);
double Temper_probabilities(int n, long i, long j);
BOOL   TestImpRatio (PGIBBSDATA pgd,  BOOL bExptIsDep,
                     double dLnKern,  double dLnKernNew,
                     double dLnPrior, double dLnPriorNew,
                     double dLnLike,  double dLnLikeNew,
                     double dLnData,  double dLnDataNew);
BOOL   TestImpRatioTemper (PGIBBSDATA pgd,  BOOL bExptIsDep, 
                           double dLnKern,  double dLnKernNew,
                           double dLnPrior, double dLnPriorNew,
                           double dLnLike,  double dLnLikeNew,
                           double dLnData,  double dLnDataNew, long indexT);
BOOL   TestTemper (PGIBBSDATA pgd, long indexT, long indexT_new, 
                   double dLnPrior, double dLnData, long attemptedTempe, 
                   double pseudo, double pseudonew);
void   TraverseLevels (PLEVEL plevel,
                       void (*routinePtr)(PLEVEL plevel, char **args), ...);
int    TraverseLevels1 (PLEVEL plevel,
                        int (*routinePtr)(PLEVEL plevel, char **args), ...);
void   WriteHeader (PLEVEL plevel, char **args);
void   WriteMCVars (PLEVEL plevel, PFILE pOutFile);

#define MH_H_DEFINED
#endif  /* _MH_H_ */

/* End */
