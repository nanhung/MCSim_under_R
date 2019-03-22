/* mh.c

   Written by Frederic Bois
   5 January 1996

   Copyright (c) 1996-2017 Free Software Foundation, Inc.

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

   This file calls at least one GNU Scientific Library function, if the
   symbol HAVE_LIBGSL is defined in config.h (or in the makefile).
   Otherwise, the corresponding features are disabled (and the program
   exits with an error message).
*/

#include <stdio.h>
#include <string.h>
#include <inttypes.h>

#include "mh.h"
#include "lsodes.h"
#include "lexerr.h"
#include "simmonte.h"

#ifdef HAVE_LIBGSL
#include <gsl/gsl_cdf.h>
#endif


/* Function -------------------------------------------------------------------
   AnnounceMarkov

   Print out the type of simulations to be performed.
*/
void AnnounceMarkov (int nSimTypeFlag, long nIter)
{

  switch (nSimTypeFlag) {

  case 0:
    printf ("\nDoing %ld Metropolis within Gibbs simulation",
            nIter);
    printf ((nIter != 1 ? "s\n" : "\n"));
    break;

  case 1:
    printf ("\nPrinting data and predictions for the last line of the "
            "restart file\n");
    break;

  case 2:
    printf ("\nDoing %ld Metropolis simulation",
            nIter);
    printf ((nIter != 1 ? "s\n" : "\n"));
    break;

  case 3:
    printf ("\nDoing %ld Metropolis within Gibbs posterior "
            "tempered simulation",
            nIter);
    printf ((nIter != 1 ? "s\n" : "\n"));

    break;

  case 4:
    printf ("\nDoing %ld Metropolis within Gibbs likelihood "
            "tempered simulation",
            nIter);
    printf ((nIter != 1 ? "s\n" : "\n"));

    break;

  case 5:
    printf("\nDoing Stochastic optimization\n");
    break;
  }

} /* AnnounceMarkov */


/* Function -------------------------------------------------------------------
   CalculateTotals

   Find total prior, likelihood for all MC vars
   Called from TraverseLevels
*/

void CalculateTotals (PLEVEL plevel, char **args)
{
  PANALYSIS panal = (PANALYSIS)args[0];
  double *pdLnPrior = (double*)args[1];

  long n;

  for (n = 0; n < plevel->nMCVars; n++) {
    *pdLnPrior += LnDensity(plevel->rgpMCVars[n], panal);
  }

} /* CalculateTotals */


/* Function -------------------------------------------------------------------
   CheckForFixed

   It is possible for an MC var to be fixed by a subsequent `=' statement in
   the level just below the level at which it is declared in the input file;
   This routine marks such vars as "Fixed" and set their dVal to the
   prescribed number.
   Called from TraverseLevels
*/

void CheckForFixed (PLEVEL plevel, char **args)
{
  long    n, m;
  PMCVAR  pMCVar;
  PVARMOD pFVar;

  for (n = 0; n < plevel->nMCVars; n++) {
    pMCVar = plevel->rgpMCVars[n];
    for (m = 0; m < plevel->nFixedVars; m++) {
      pFVar = plevel->rgpFixedVars[m];
      if (pMCVar->hvar == pFVar->hvar) {
        pMCVar->bIsFixed = TRUE;
        if (IsInput (pFVar->hvar)) {
          printf("Error: a sampled parameter cannot be assigned an input\n");
          exit(0);
        }
        else
          pMCVar->dVal = pFVar->uvar.dVal;
      }
    }
  }

} /* CheckForFixed */


/* Function -------------------------------------------------------------------
   CheckPrintStatements

   Tries to check that 'Print' and 'Data' statements are consistent.
   Note that experiments must have outputs; this is checked in
   PrepareOutSpec in siminit.c.
   Called from TraverseLevels
*/

void CheckPrintStatements (PLEVEL plevel, char **args)
{
  PANALYSIS panal = (PANALYSIS) args[0];
  POUTSPEC  pos;
  /* PMCVAR    pMCVar;
  BOOL      bFound, bOK; */
  long      i, j /* , k, l, dCount, dIndex */;

  if (plevel->pexpt == NULL)
    return;

  pos = &(plevel->pexpt->os);

  /* check that the same var does not appear in 2 or more print statements */
  for (j = 0; j < pos->nOutputs; j++)
    for (i = j+1; i < pos->nOutputs; i++)
      if (pos->phvar_out[j] == pos->phvar_out[i])
        ReportRunTimeError (panal, RE_DUPVARINEXPRT | RE_FATAL,
                            pos->pszOutputNames[i], "Print");

  /* check that the same var does not appear in 2 or more data statements */
  for (j = 0; j < pos->nData; j++)
    for (i = j+1; i < pos->nData; i++)
      if (pos->phvar_dat[j] == pos->phvar_dat[i])
        ReportRunTimeError (panal, RE_DUPVARINEXPRT | RE_FATAL,
                            pos->pszDataNames[i], "Data");

  /* check that print and matching data lists are of equal length */
  for (j = 0; j < pos->nOutputs; j++)
    for (i = 0; i < pos->nData; i++)
      if ((pos->phvar_out[j] == pos->phvar_dat[i]) &&
          (pos->pcOutputTimes[j] != pos->pcData[i])) {
        printf("\nError: unequal times in Print and Data statements for %s\n"
               "Exiting.\n\n", pos->pszOutputNames[j]);
        exit(0);
      }

} /* CheckPrintStatements */


/* Function -------------------------------------------------------------------
   CloneLikes

   Called from TraverseLevels
   First copy the likelihoods in the current plistLikes in the rgpLikes
   down to the next lower levels. Then eventually override them by definitions
   in rgpLikes at this level.
*/
void CloneLikes (PLEVEL plevel, char **args)
{
  long nLikes;
  long i, j, k;
  PLEVEL pLower;
  PMCVAR pClone;
  PMCVAR pMCVar;
  BOOL bFound;

  for (i = 0; i < plevel->iInstances; i++) {
    pLower = plevel->pLevels[i];
    /* The number of likelihoods is overestimated by this, but that will do
       for now: */
    pLower->nLikes = nLikes = plevel->nLikes + ListLength(plevel->plistLikes);

    if (pLower->nLikes != 0) {
      if (!(pLower->rgpLikes =
             (PMCVAR*) malloc (pLower->nLikes * sizeof(PMCVAR))))
        ReportError (NULL, RE_OUTOFMEM | RE_FATAL, "CloneLikes", NULL);
    }
  }

  /* Copy the likelihoods in the plistLikes at this level to the
     rgpLikes of the lower levels */
  nLikes = 0;
  ForAllList3 (plevel->plistLikes, &CloneLikesL, plevel, &nLikes, NULL);

  /* Note: nLikes now equals the number of likelihoods in the current
     plistLikes */

  /* Copy the likelihoods in rgpLikes at this level down to the lower
     levels, if they are not already defined */
  for (i = 0; i < plevel->iInstances; i++) {
    pLower = plevel->pLevels[i];

    for (j = 0; j < plevel->nLikes; j++) {
      pMCVar = plevel->rgpLikes[j];

      /* if that likelihood is already in pLower->rgpLikes forget it */
      bFound = FALSE;
      k = 0;
      while ((k < nLikes) && (!bFound)) {
        bFound = (pMCVar->hvar == pLower->rgpLikes[k]->hvar);
        if (!bFound) k++;
      }
      if (!bFound) {
        if (!(pClone = (PMCVAR) malloc (sizeof (MCVAR))))
          ReportError (NULL, RE_OUTOFMEM | RE_FATAL, "CloneLikes", NULL);

        memcpy (pClone, pMCVar, sizeof (MCVAR));
        pLower->rgpLikes[nLikes+j] = pClone;
      }
    }
  }

} /* CloneLikes */


/* Function -------------------------------------------------------------------
   CloneLikesL

   Clone and copy the members of a plistLikes from a level to the rgpLikes
   arrays of the levels below
   Called from ForAllList3
*/
void CloneLikesL (PVOID pData, PVOID pUser1, PVOID pUser2, PVOID pUser3)
{
  PMCVAR pMCVar = (PMCVAR) pData;
  PLEVEL plevel = (PLEVEL) pUser1;
  long   *pnLikes = (long *) pUser2;
  long   i;
  PLEVEL pLower;
  PMCVAR pClone;

  ++pMCVar->iDepth;
  for (i = 0; i < plevel->iInstances; i++) {
    pLower = plevel->pLevels[i];
    if (!(pClone = (PMCVAR) malloc (sizeof (MCVAR))))
      ReportError (NULL, RE_OUTOFMEM | RE_FATAL, "CloneLikeL", NULL);
    memcpy (pClone, pMCVar, sizeof (MCVAR));
    pLower->rgpLikes[*pnLikes] = pClone;
  }
  ++(*pnLikes);

} /* CloneLikesL */


/* Function -------------------------------------------------------------------
   CloneMCVars

   Called from TraverseLevels
   For all MC vars in list at given level, add to arrays of all instances of
   next (lower) level; the instances are created by the cloning routine
   CloneMCVarsL
*/
void CloneMCVars (PLEVEL plevel, char **args)
{
  long nMCVars = ListLength(plevel->plistMCVars);
  long n;
  PLEVEL pLower;

  for (n = 0; n < plevel->iInstances; n++) {
    pLower = plevel->pLevels[n];
    pLower->nMCVars = nMCVars;
    if (nMCVars != 0) {
      if (!(pLower->rgpMCVars = (PMCVAR*) malloc (nMCVars * sizeof(PMCVAR))))
        ReportError (NULL, RE_OUTOFMEM | RE_FATAL, "CloneMCVars", NULL);
    }
  }

  nMCVars = 0;
  ForAllList3 (plevel->plistMCVars, &CloneMCVarsL, plevel, &nMCVars, NULL);

} /* CloneMCVars */


/* Function -------------------------------------------------------------------
   CloneMCVarsL

   Clone and copy the members of a plistMCVars
   Called from ForAllList3
*/
void CloneMCVarsL (PVOID pData, PVOID pUser1, PVOID pUser2, PVOID pUser3)
{
  PMCVAR pMCVar = (PMCVAR) pData;
  PLEVEL plevel = (PLEVEL) pUser1;
  long   *pnMCVars = (long *) pUser2;
  long   i;
  PLEVEL pLower;
  PMCVAR pClone;

  ++pMCVar->iDepth;
  for (i = 0; i < plevel->iInstances; i++) {
    pLower = plevel->pLevels[i];
    if (!(pClone = (PMCVAR) malloc(sizeof (MCVAR))))
      ReportError(NULL, RE_OUTOFMEM | RE_FATAL, "CloneMCVarsL", NULL);
    memcpy (pClone, pMCVar, sizeof (MCVAR));
    pClone->plistDependents = InitList();
    pLower->rgpMCVars[*pnMCVars] = pClone;
  }
  ++(*pnMCVars);

} /* CloneMCVarsL */


/* Function -------------------------------------------------------------------
   CloseMarkovFiles

   Closes output files associated with the Markov sampler.
   The restart file has already been closed by the ReadRestart

*/
void CloseMarkovFiles (PANALYSIS panal)
{
  if (panal->gd.pfileOut) {
    fclose (panal->gd.pfileOut);
    printf ("\nWrote results to \"%s\"\n", panal->gd.szGout);
  }

} /* CloseMarkovFiles */


/* Function -------------------------------------------------------------------
   ConvertLists

   Converts lists to arrays
   Called from TraverseLevels
*/
void ConvertLists(PLEVEL plevel, char **args)
{
  PANALYSIS panal = (PANALYSIS)args[0];
  long m, n;
  PMCVAR pMCVar;

  if (plevel->pexpt == NULL)
    ListToPVArray (panal, plevel->plistVars, &plevel->nFixedVars,
                   &plevel->rgpFixedVars);
  else
    ListToPVArray (panal, plevel->pexpt->plistParmMods, &plevel->nFixedVars,
                   &plevel->rgpFixedVars);

  for (n = 0; n < plevel->nMCVars; n++) {
    pMCVar = plevel->rgpMCVars[n];
    ListToPMCArray (panal, pMCVar->plistDependents,
                    &pMCVar->nDependents, &pMCVar->rgpDependents);

    /* if there are no dependent, experiments are most likely to be dependent
       on that parameter (complete checking would require  verifying the
       model file - FB 2 May 1998 */
    if (pMCVar->nDependents == 0)
      pMCVar->bExptIsDep = TRUE;
    else {
      /* if there are dependent, experiments may be dependent on that
         parameter if they are not shielded by an intermediate instance of
         the same parameter */
      m = 0;
      while ((m < pMCVar->nDependents) &&
             !(pMCVar->bExptIsDep = (strcmp(pMCVar->rgpDependents[m]->pszName,
                                            pMCVar->pszName) ? TRUE : FALSE)))
        m++;
    }
  }

} /* ConvertLists */


/* Function -------------------------------------------------------------------
   SetdInvTemperatures

   Set default values for inverse temperatures of tempered MCMC algorithms, if
   the user has not specified values.
*/
void SetInvTemperatures (PGIBBSDATA pgd) {

  /* inverse temperatures not set by the user: use the default */
  if (((pgd->nSimTypeFlag == 3) || (pgd->nSimTypeFlag == 4)) &&
      (pgd->nInvTemperatures == 0)) {

    int i;

    pgd->nInvTemperatures = NTEMP;

    /* allocate inverse temperature array */
    if (!(pgd->rgInvTemperatures = InitdVector (NTEMP)))
      ReportError (NULL, RE_OUTOFMEM | RE_FATAL, "SetInvTemperatures", NULL);

    /* allocate working arrays */
    if (!(pgd->rgdlnPi = InitdVector (NTEMP)) ||
        !(pgd->rglTemp = InitlVector (NTEMP)))
      ReportError (NULL, RE_OUTOFMEM | RE_FATAL, "SetInvTemperatures", NULL);

    /* initialize rglTemp, counting perk visits */
    for (i = 0; i < NTEMP; i++)
      pgd->rglTemp[i] = 0;

    /* initialize the inverse temperatures */
    pgd->rgInvTemperatures[0] = 0.4;                         /* hot */

    for (i = 1; i < NTEMP - 1; i++)
      pgd->rgInvTemperatures[i] = pow(0.8, (NTEMP - 1 - i)); /* medium... */

    pgd->rgInvTemperatures[NTEMP - 1] = 1;                   /* cold */

  } /* end if */

} /* SetInvTemperatures */


/* Function -------------------------------------------------------------------
   DoMarkov

   Core routine of the MCMC sampler
*/
void DoMarkov (PANALYSIS panal)
{
  PGIBBSDATA pgd = &panal->gd;
  PLEVEL     pLevel0 = panal->pLevels[0];
  long       nThetas, nUpdateAt, nTotal;
  long       i, iter = 0;
  long       nIter = pgd->nMaxIter;  /* scheduled iterations of the sampler */
  double     *pdMCVarVals  = NULL;   /* read in values of thetas */
  double     *pdSum        = NULL;   /* column sums of read thetas */
  double     **prgdSumProd = NULL;   /* column sums of thetas cross-products */
  double     dTmp, dLnPrior = 0, dLnData = 0;

  AnnounceMarkov (pgd->nSimTypeFlag, nIter);

  SetInvTemperatures (pgd); /* (if needed) */

  OpenMarkovFiles (panal);

  ReadDataFile (panal);     /* (if needed) */

  /* MC variables must be placed in arrays at the next lower level */
  TraverseLevels (pLevel0, CloneMCVars, NULL);

  /* Likelihoods must percolate down to the experiment level */
  TraverseLevels (pLevel0, CloneLikes, NULL);

  /* Find the parents and dependents of the MC vars */
  TraverseLevels (pLevel0, FindMCParents,    panal, NULL);
  TraverseLevels (pLevel0, FindMCDependents, panal, NULL);

  /* Find the parents of the parameters in likelihood statements */
  TraverseLevels (pLevel0, FindLikeParents, panal, NULL);

  /* Now that we have the MC vars right, write the output file header
     unless it's just a run for fit checking */
  if (pgd->nSimTypeFlag != 1) {
    fprintf (pgd->pfileOut, "iter\t");
    TraverseLevels (pLevel0, WriteHeader, panal, pgd->pfileOut, NULL);
    if ((pgd->nSimTypeFlag == 3) || (pgd->nSimTypeFlag == 4)) { /* tempered */
      fprintf (pgd->pfileOut, "IndexT\t");
      for (i = 0; i < pgd->nInvTemperatures; i++)
        fprintf (pgd->pfileOut, "LnPseudoPrior(%ld)\t",i);
      fprintf (pgd->pfileOut, "LnPrior\tLnData\tLnPosterior\n");
    }
    else {
      fprintf (pgd->pfileOut, "LnPrior\tLnData\tLnPosterior\n");
    }
    fflush (pgd->pfileOut);
  }

  /* Convert the rest of the lists to arrays */
  TraverseLevels (pLevel0, ConvertLists, panal, NULL);

  /* Check for MC vars that have been fixed */
  TraverseLevels (pLevel0, CheckForFixed, NULL);

  /* Check variables in statements of type
     'Distrib(<var1>, <distrib>, Prediction, <var2>)'
     for identical 'Print' statements */
  TraverseLevels (pLevel0, CheckPrintStatements, panal, NULL);

  /* Print out the structure for debugging */
  if (panal->bDependents) {
    TraverseLevels (pLevel0, PrintDeps, NULL);
    CloseMarkovFiles (panal);
    return;
  }

  /* Change the MC vars hvar pointers from pointing to model parameters to
     pointing to the parents' dVal */
  TraverseLevels (pLevel0, SetPointers, panal, NULL);

  /* Get the initial values of the MC vars by reading th restart file if
     one is given */
  if (pgd->szGrestart) {

    /* Read them from the restart file */
    nThetas = 0;
    TraverseLevels (pLevel0, GetNumberOfMCVars, &nThetas, NULL);

    if ( !(pdMCVarVals = InitdVector (nThetas)) ||
         !(pdSum       = InitdVector (nThetas)) ||
         !(prgdSumProd = InitdMatrix (nThetas, nThetas))) {
      /* FB 6/4/98 */
      /* lack of memory can be severe here, free immediately some space */
      if (pdMCVarVals) free (pdMCVarVals);
      if (pdSum) free (pdSum);
      ReportRunTimeError (panal, RE_OUTOFMEM | RE_FATAL, "DoMarkov");
    }

    /* Read the starting values in the order they are printed and
       close the file when finished */

    if ((pgd->nSimTypeFlag == 3) || (pgd->nSimTypeFlag == 4))
      ReadRestartTemper (pgd->pfileRestart, nThetas, pgd->nInvTemperatures,
                         pdMCVarVals, pdSum, prgdSumProd, &iter,
                         &pgd->indexT, pgd->rgdlnPi);
    else
      ReadRestart (pgd->pfileRestart, nThetas, pdMCVarVals, pdSum, prgdSumProd,
                   &iter);

    /* Set the dVals of the variables to the values read in */
    nThetas = 0;
    if (!TraverseLevels1 (pLevel0, SetMCVars, pdMCVarVals, &nThetas, NULL)) {
      printf ("\nError: some read-in parameters are out of bounds - "
              "Exiting\n\n");
      exit(0);
    }

    /* If nSimTypeFlag is 1 print just the predictions and data and exit */
    if (pgd->nSimTypeFlag == 1) {
      PrintAllExpts (pLevel0, panal, pgd->pfileOut);
      CloseMarkovFiles (panal);
      return;
    }

    /* If component jumps, tempered or stochastic optimization */
    if ((pgd->nSimTypeFlag == 0) || (pgd->nSimTypeFlag >= 3)) {

      char szKernelFile[MAX_FILENAMESIZE+7];
      sprintf(szKernelFile, "%s%s", panal->gd.szGrestart, ".kernel");

      FILE *pfile = fopen (szKernelFile, "r");

      if (pfile) {
        printf("Reading kernel file %s\n", szKernelFile);
        TraverseLevels (pLevel0, ReadKernel, pfile, NULL);
      }
      else {
        /* Set the jumping kernel's SD */
        TraverseLevels (pLevel0, SetKernel, 1, pdMCVarVals, NULL);

        /* Free the now useless array */
        free (pdMCVarVals);
      }
    }
    else { /* vector jumping, initialize prior */
      TraverseLevels (pLevel0, CalculateTotals, panal, &dLnPrior, NULL);
    }

    /* Initialize the predictions arrays by running all experiments and
       save the likelihood */
    if (!RunAllExpts (panal, &dLnData)) {
      /* error, cannot compute */
      printf ("\nError: cannot compute at the starting point - Exiting\n\n");
      exit(0);
    }
  }

  else { /* no restart file, init by sampling */

    /* If nSimTypeFlag is 1 print an error message and exit */
    if (pgd->nSimTypeFlag == 1) {
      printf ("\nError: a restart file must be given to print data and"
              "         predictions - Exiting.\n\n");
      exit (0);
    }

    /* Set the jumping kernel's SD */
    TraverseLevels (pLevel0, SetKernel, 2, pdMCVarVals, NULL);

    /* Initialize the thetas by sampling from the prior,
       and write them out at the same time */
    fprintf (pgd->pfileOut, "0\t");
    TraverseLevels (pLevel0, InitMCVars, pgd->pfileOut, NULL);

    /* Calculate the total prior */
    TraverseLevels (pLevel0, CalculateTotals, panal, &dLnPrior, NULL);
    /* Initialize the predictions arrays by running all experiments */
    if (!RunAllExpts (panal, &dLnData)) {
      /* Error, cannot compute */
      printf ("\nError: cannot compute at the starting point - Exiting\n\n");
      exit(0);
    }

    if ((pgd->nSimTypeFlag == 3) || (pgd->nSimTypeFlag == 4)) {
      /* Output indexT and pseudo-priors */
      fprintf (pgd->pfileOut, "0\t");
      for (i = 0; i < pgd->nInvTemperatures; i++)
        fprintf (pgd->pfileOut, "%e\t", pgd->rgdlnPi[i]);
      fflush (pgd->pfileOut);
    }

    /* Output prior and likelihood */
    fprintf (pgd->pfileOut, "%e\t%e\t%e\n", dLnPrior, dLnData,
             dLnPrior + dLnData);
    fflush (pgd->pfileOut);
  }

  /* Save the data likelihoods */
  TraverseLevels1 (pLevel0, SaveLikelihoods, NULL);

  /* Initializations are finished, let's do the iterations */
  nTotal = UPDATE_BASE;
  nUpdateAt = iter + nTotal; /* kernel will be updated at that iteration */

  while (iter < nIter) {

    /* Output to screen, eventually */
    if (panal->bPrintIter && ((iter+1) % 100 == 0))
      printf("Iteration %ld\n", iter + 1);

    /* Start output to file, eventually */
    if (((iter + 1) % pgd->nPrintFreq == 0) &&
        (iter >= pgd->nMaxIter - pgd->nPrintIter))
      fprintf (pgd->pfileOut, "%ld\t", iter + 1);

    if ((pgd->nSimTypeFlag == 0) || (pgd->nSimTypeFlag == 5)) {
      /* component by component jumps or stochastic optimization */

      TraverseLevels (pLevel0, SampleThetas, panal, pgd, &iter, &nUpdateAt,
                      &nTotal, NULL);

      /* Output log-densities, eventually */
      if (((iter + 1) % pgd->nPrintFreq == 0) &&
          (iter >= pgd->nMaxIter - pgd->nPrintIter)) {
        dLnPrior = 0.0;
        TraverseLevels (pLevel0, CalculateTotals, panal, &dLnPrior, NULL);
        dLnData = 0.0;
        TraverseLevels1 (pLevel0, SumAllExpts, &dLnData, NULL);
        fprintf (pgd->pfileOut, "%e\t%e\t%e\n", dLnPrior, dLnData,
               dLnPrior + dLnData);
        fflush (pgd->pfileOut);
      }
    }
    else {

      if ((pgd->nSimTypeFlag == 3) || (pgd->nSimTypeFlag == 4)) {
        /* tempered */
        TraverseLevels (pLevel0, SampleThetasTempered, panal, pgd, &iter,
                        &nUpdateAt,  &nTotal, &pgd->indexT, NULL);

        dLnPrior = 0.0;
        TraverseLevels (pLevel0, CalculateTotals, panal, &dLnPrior, NULL);
        dLnData = 0.0;
        TraverseLevels1 (pLevel0, SumAllExpts, &dLnData, NULL);

        /* Robbins-Monro updating of the pseudo prior */
        for (i = 0; i < pgd->nInvTemperatures; i++) {
          dTmp = pgd->dCZero / (double) (iter + pgd->dNZero);
          if (i == pgd->indexT)
            pgd->rgdlnPi[i] -= dTmp;
          else
            pgd->rgdlnPi[i] += dTmp / (double) pgd->nInvTemperatures;
        }

        /* Output indexT and log-densities, eventually */
        if (((iter + 1) % pgd->nPrintFreq == 0) &&
            (iter >= pgd->nMaxIter - pgd->nPrintIter)) {
          fprintf (pgd->pfileOut, "%ld\t", pgd->indexT);
          for (i = 0; i < pgd->nInvTemperatures; i++)
            fprintf (pgd->pfileOut, "%e\t", pgd->rgdlnPi[i]);
          fprintf (pgd->pfileOut, "%e\t%e\t%e\n", dLnPrior, dLnData,
                   dLnPrior + dLnData);
          fflush (pgd->pfileOut);
        }

        /* update population count of current temperature */
        pgd->rglTemp[pgd->indexT] = pgd->rglTemp[pgd->indexT]+1;

        /* test the temperature and change indexT if necessary */
        pgd->indexT = SampleTemperature (pgd, dLnPrior, dLnData);
      }
      else { /* vector jumps */

        SampleThetaVector (pLevel0, panal, nThetas, pdMCVarVals, pdSum,
                           prgdSumProd, iter, nUpdateAt, nTotal, &dLnPrior,
                           &dLnData);

        /* Output, eventually */
        if (((iter + 1) % pgd->nPrintFreq == 0) &&
            (iter >= pgd->nMaxIter - pgd->nPrintIter)) {

          for (i = 0; i < nThetas; i++)
            fprintf (pgd->pfileOut, "%5g\t", pdMCVarVals[i]);

          fprintf (pgd->pfileOut, "%e\t%e\t%e\n", dLnPrior, dLnData,
                   dLnPrior + dLnData);
          fflush (pgd->pfileOut);
        }
      }
    }
    /* Adjust the update time, eventually */
    if (iter == nUpdateAt) {
      nTotal = nTotal * 3 / 2;
      nUpdateAt = iter + nTotal;
    }

    /* Increment the iteration counter */
    iter = iter + 1;

  } /* while iter */

  /* If component jumps, tempered or stochastic optimization: save kernel */
  if ((pgd->nSimTypeFlag == 0) || (pgd->nSimTypeFlag >= 3)) {

    char szKernelFile[MAX_FILENAMESIZE+7];
    sprintf(szKernelFile, "%s%s", panal->gd.szGout, ".kernel");

    FILE *pfile = fopen (szKernelFile, "w");
    if (!pfile) {
      printf ("Cannot create kernel saving file '%s'\n", panal->gd.szGdata);
      exit (0);
    }

    /* Write out the jumping kernel SDs */
    TraverseLevels (pLevel0, WriteKernel, pfile, NULL);

    fprintf(pfile, "\n");
    fclose(pfile);
    printf("\nWrote kernel SDs to \"%s\"", szKernelFile);
  }

  CloseMarkovFiles (panal);

  if ((pgd->nSimTypeFlag == 3) || (pgd->nSimTypeFlag == 4)) {
    printf ("\nInverse temperatures:\n");
    for (i = 0; i < pgd->nInvTemperatures - 1; i++)
      printf ("%f\t", pgd->rgInvTemperatures[i]);
    printf ("%f\n", pgd->rgInvTemperatures[pgd->nInvTemperatures - 1]);

    printf ("\nInverse temperatures' count:\n");
    for (i = 0; i < pgd->nInvTemperatures - 1; i++)
      printf ("%ld\t", pgd->rglTemp[i]);
    printf ("%ld\n\n", pgd->rglTemp[pgd->nInvTemperatures - 1]);
  }

} /* DoMarkov */


/* Function -------------------------------------------------------------------
   FindLikeParents

   Called from TraverseLevels
   Find the parents of the MC vars corresponding to likelihood at this level
   by looking at sampled parameters at this and all previous levels. Note that
   dependencies are not checked among likelihoods.
*/
void FindLikeParents (PLEVEL plevel, char **args)
{
  PANALYSIS panal = (PANALYSIS)args[0];
  long      k, l, m, n;
  PLEVEL    pPrevLev;
  PMCVAR    pMCVar1, pMCVar2;
  BOOL      bFound;

  /* Set the current level array as we pass through */
  panal->pCurrentLevel[plevel->iDepth] = plevel;

  for (k = 0; k < plevel->nLikes; k++) {
    pMCVar1 = plevel->rgpLikes[k]; /* current likelihood */

    for (l = 0; l < 4; l++) {
      if (pMCVar1->iParmType[l] == MCVP_PARM) {
        /* if there is a parent, find it, but parents must appear before
           current pMCVar */
        bFound = FALSE;

        /* First, this level */
        for (m = 0; m < plevel->nMCVars; m++) { /* for all previous distrib */
          pMCVar2 = plevel->rgpMCVars[m];
          if (pMCVar1->hParm[l] == pMCVar2->hvar) {
            pMCVar1->pMCVParent[l] = pMCVar2;
            bFound = TRUE;
          }
        } /* for m */

        /* Now, all previous levels */
        if (!bFound) {
          for (n = plevel->iDepth - 1; n >= 0; n--) {
            pPrevLev = panal->pCurrentLevel[n];

            for (m = 0; m < pPrevLev->nMCVars; m++) {
              pMCVar2 = pPrevLev->rgpMCVars[m];
              if (pMCVar1->hParm[l] == pMCVar2->hvar) {
                pMCVar1->pMCVParent[l] = pMCVar2;
                bFound = TRUE;
              }
            } /* for m */
          } /* for n */
        }

        if (!bFound) { /* oops, parent not found, error */
          printf ("\n"
                  "Error: parent in position %ld of %s must be\n"
                  "       declared before it when creating\n"
                  "       sampling dependencies - Exiting.\n\n",
                  l, pMCVar1->pszName);
          exit(0);
        }

      }
    } /* for l */
  } /* for k */

} /* FindLikeParents */


/* Function -------------------------------------------------------------------
   FindMCDependents

   Called from TraverseLevels
   Find the direct dependents of the MC vars at this level by looking at this
   and all lower levels
*/
void FindMCDependents (PLEVEL plevel, char **args)
{
  long   i, j;
  PMCVAR pMCVar;

  for (i = 0; i < plevel->nMCVars; i++) {
    pMCVar = plevel->rgpMCVars[i];
    for (j = 0; j < 4; j++)
      if ((pMCVar->pMCVParent[j] != NULL) &&
          (pMCVar->pMCVParent[j]->hvar == pMCVar->hParm[j]))
        QueueListItem(pMCVar->pMCVParent[j]->plistDependents, pMCVar);
  }

} /*FindMCDependents */


/* Function -------------------------------------------------------------------
   FindMCParents

   Called from TraverseLevels
   Find the parents of the MC vars at this level by looking at this and
   all previous levels.
*/
void FindMCParents (PLEVEL plevel, char **args)
{
  PANALYSIS panal = (PANALYSIS)args[0];
  long      k, l, m, n;
  PLEVEL    pPrevLev;
  PMCVAR    pMCVar1, pMCVar2;
  BOOL      bFound;

  /* Set the current level array as we pass through */
  panal->pCurrentLevel[plevel->iDepth] = plevel;

  for (k = 0; k < plevel->nMCVars; k++) {
    pMCVar1 = plevel->rgpMCVars[k]; /* current distrib */

    for (l = 0; l < 4; l++) {
      if (pMCVar1->iParmType[l] == MCVP_PARM) {
        /* if there is a parent, find it, but parents must appear before
           current pMCVar */
        bFound = FALSE;

        /* First, this level */
        for (m = 0; m < k; m++) { /* for all previous distrib */
          pMCVar2 = plevel->rgpMCVars[m];
          if (pMCVar1->hParm[l] == pMCVar2->hvar) {
            pMCVar1->pMCVParent[l] = pMCVar2;
            bFound = TRUE;
          }
        } /* for m */

        /* Now, all previous levels */
        if (!bFound) {
          for (n = plevel->iDepth - 1; n >= 0; n--) {
            pPrevLev = panal->pCurrentLevel[n];

            for (m = 0; m < pPrevLev->nMCVars; m++) {
              pMCVar2 = pPrevLev->rgpMCVars[m];
              if (pMCVar1->hParm[l] == pMCVar2->hvar) {
                pMCVar1->pMCVParent[l] = pMCVar2;
                bFound = TRUE;
              }
            } /* for m */
          } /* for n */
        }

        if (!bFound) { /* oops, parent not found, error */
          printf ("\n"
                  "Error: parent in position %ld of %s must be\n"
                  "       declared before it when creating\n"
                  "       sampling dependencies - Exiting.\n\n",
                  l, pMCVar1->pszName);
          exit(0);
        }

      }
    } /* for l */
  } /* for k */

} /* FindMCParents */


/* Function -------------------------------------------------------------------
   GetNumberOfMCVars

   Find the total number of MC vars
   Called from TraverseLevels
*/
void GetNumberOfMCVars (PLEVEL plevel, char **args)
{
  long *pnThetas = (long *) args[0];

  *pnThetas += plevel->nMCVars;

} /* GetNumberOfMCVars */


/* Function -------------------------------------------------------------------
   InitMCVars

   Sample initial values of thetas if not fixed
   Called from TraverseLevels
*/
void InitMCVars(PLEVEL plevel, char **args)
{
  FILE *pOutFile = (FILE*)args[0];
  long n;

  for (n = 0; n < plevel->nMCVars; n++)
    if ( !(plevel->rgpMCVars[n]->bIsFixed))
      CalculateOneMCParm (plevel->rgpMCVars[n]);

  /* Write out the sampled values */
  WriteMCVars (plevel, pOutFile);

} /* InitMCVars */


/* Function -------------------------------------------------------------------
   ListToPMCArray
*/
void ListToPMCArray (PANALYSIS panal, PLIST plist,
                     long *pnMCVars, PMCVAR **rgpMCVars)
{
  if ((*pnMCVars = ListLength(plist)) == 0)
    return;

  if (!(*rgpMCVars = (PMCVAR*) malloc (*pnMCVars * sizeof(PMCVAR))))
    ReportRunTimeError (panal, RE_OUTOFMEM | RE_FATAL, "ListToPMCArray");

  *pnMCVars = 0;
  ForAllList3 (plist, &ListToPMCArrayL, pnMCVars, *rgpMCVars, NULL);

} /*ListToPMCArray */


/* Function -------------------------------------------------------------------
   ListToPMCArrayL
*/
void ListToPMCArrayL (PVOID pData, PVOID pUser1, PVOID pUser2, PVOID pUser3)
{
  PMCVAR pMCVar = (PMCVAR)pData;
  long *pnMCVars = (long*)pUser1;
  PMCVAR *rgpMCVars = (PMCVAR*)pUser2;

  rgpMCVars[(*pnMCVars)++] = pMCVar;

} /* ListToPMCArrayL */


/* Function -------------------------------------------------------------------
   ListToPVArray
*/
void ListToPVArray (PANALYSIS panal, PLIST plist,
                    long *pnFixedVars, PVARMOD **rgpFixedVars)
{
  if ((*pnFixedVars = ListLength (plist)) == 0)
    return;

  if (!(*rgpFixedVars = (PVARMOD*) malloc (*pnFixedVars * sizeof(PVARMOD))))
    ReportRunTimeError (panal, RE_OUTOFMEM | RE_FATAL, "ListToPVArray");

  *pnFixedVars = 0;
  ForAllList3 (plist, &ListToPVArrayL, pnFixedVars, *rgpFixedVars, NULL);

} /*ListToPVArray */


/* Function -------------------------------------------------------------------
   ListToPVArrayL
*/
void ListToPVArrayL (PVOID pData, PVOID pUser1, PVOID pUser2, PVOID pUser3)
{
  PVARMOD pVar = (PVARMOD)pData;
  long    *pnFixedVars = (long*)pUser1;
  PVARMOD *rgpFixedVars = (PVARMOD*)pUser2;

  rgpFixedVars[(*pnFixedVars)++] = pVar;

} /* ListToPVArrayL */


/* Function -------------------------------------------------------------------
   LnDensity

   Returns the log of the (exact) density of variate under its distribution.
*/
#define LNSQRT2PI  0.918938533204672669541
#define LNINVERPI -1.144729885849400163877
#define LN2OVERPI -0.4515827052894548221396

double LnDensity (PMCVAR pMCVar, PANALYSIS panal)
{
  double dTmp = 0, density;
  double dTmp2, dTmp3, dTmp4;
  double dParm1 = *(pMCVar->pdParm[0]);
  double dParm2 = *(pMCVar->pdParm[1]);
  double dMin   = *(pMCVar->pdParm[2]);
  double dMax   = *(pMCVar->pdParm[3]);
  double dTheta = pMCVar->dVal;
  char str[10];

  /* This should take care of all dTheta checking */
  if (pMCVar->iType == MCV_BINOMIALBETA) {
    if (dTheta < 0) {
      printf ("Error: variate out of bounds in LnDensity\n");
      exit (0);
    }
  }
  else if (pMCVar->iType == MCV_GENLOGNORMAL ||
           pMCVar->iType == MCV_STUDENTT) {
    if (dParm1 < 0) {
      printf ("Error: parameter %g out of bounds in LnDensity\n",dParm1);
      exit (0);
    }
  }
  else {
    if (dTheta > dMax || dTheta < dMin) {
      return (NULL_SUPPORT);
    }
  }

  switch (pMCVar->iType) {

    case MCV_UNIFORM:
      if (dTheta > dParm2 || dTheta < dParm1)
        return (NULL_SUPPORT);
      if (dParm2 <= dParm1)
        ReportRunTimeError (panal, RE_BADUNIFORMDIST | RE_FATAL,
                            pMCVar->pszName, "LnDensity");
      return -log(dParm2 - dParm1);

    case MCV_LOGUNIFORM:
      if (dTheta > dParm2 || dTheta < dParm1)
        return (NULL_SUPPORT);
      if (dParm2 <= dParm1)
        ReportRunTimeError (panal, RE_BADUNIFORMDIST | RE_FATAL,
                            pMCVar->pszName, "LnDensity");
      return -log(dTheta * (dParm2 - dParm1));

    case MCV_NORMALV:  dParm2 = sqrt (dParm2);
      return lnDFNormal (dTheta, dParm1, dParm2);

    case MCV_NORMALCV: dParm2 = fabs(dParm1 * dParm2);
      return lnDFNormal (dTheta, dParm1, dParm2);
 
    case MCV_NORMAL:
    case MCV_HALFNORMAL:
      return lnDFNormal (dTheta, dParm1, dParm2);

    case MCV_LOGNORMALV: dParm2 = exp(sqrt(dParm2)); /* fall thru */
    case MCV_LOGNORMAL:
      if (dParm1 <= 0.0) {
        sprintf(str, "%5.2e", dParm1);
        ReportRunTimeError(panal, RE_BADLOGNORMALMEAN | RE_FATAL,
                           pMCVar->pszName, str, "LnDensity");
      }
      return (lnDFNormal(log(dTheta),log(dParm1),log(dParm2)) - log(dTheta));

    case MCV_TRUNCNORMALV: dParm2 = sqrt (dParm2);
      return lnDFNormal (dTheta, dParm1, dParm2) -
             log(CDFNormal ((dMax - dParm1) / dParm2) -
                 CDFNormal ((dMin - dParm1) / dParm2));

    case MCV_TRUNCNORMALCV: dParm2 = fabs(dParm1 * dParm2);
      return lnDFNormal (dTheta, dParm1, dParm2) -
             log(CDFNormal ((dMax - dParm1) / dParm2) -
                 CDFNormal ((dMin - dParm1) / dParm2));

    case MCV_TRUNCNORMAL:
      if (dParm2 <= 0.0) {
        sprintf(str, "%5.2e", dParm2);
        ReportRunTimeError(panal, RE_BADNORMALSD | RE_FATAL,
                           pMCVar->pszName, str, "LnDensity");
      }
      return lnDFNormal (dTheta, dParm1, dParm2) -
             log(CDFNormal ((dMax - dParm1) / dParm2) -
                 CDFNormal ((dMin - dParm1) / dParm2));

    case MCV_TRUNCLOGNORMALV: dParm2 = exp(sqrt(dParm2)); /* fall thru */
    case MCV_TRUNCLOGNORMAL:
      if (dParm1 <= 0.0 ) {
        sprintf(str, "%5.2e", dParm1);
        ReportRunTimeError(panal, RE_BADLOGNORMALMEAN | RE_FATAL,
                           pMCVar->pszName, str, "LnDensity");
      }
      if (dParm2 <= 1.0 ) {
        sprintf(str, "%5.2e", dParm2);
        ReportRunTimeError(panal, RE_BADLOGNORMALSD | RE_FATAL,
                           pMCVar->pszName, str, "LnDensity");
      }
      dTmp = log(dParm2);
      return lnDFNormal (log(dTheta), log(dParm1), dTmp) - log(dTheta) -
             log(CDFNormal (log(dMax / dParm1) / dTmp) -
                 CDFNormal (log(dMin / dParm1) / dTmp));

    case MCV_BETA:
      return lnDFBeta (dTheta, dParm1, dParm2, dMin, dMax);

    case MCV_CHI2:
      dTmp = 0.5 * dParm1;
      return (dTmp - 1) * log(dTheta) - 0.5 * dTheta +
             dTmp * (-6.9314718056E-01) - lnGamma (dTmp);

    case MCV_BINOMIAL:
      if ((dParm1 < 0) || (dParm1 > 1)) {
        printf ("Error: bad p for binomial variate in LnDensity\n");
        exit (0);
      }
      if (dTheta > dParm2) {
        return (NULL_SUPPORT);
      }
      /* log binomial coefficient n! / (x!(n-x)!) */
      dTmp = lnGamma (dParm2 + 1) - lnGamma (dTheta + 1) -
             lnGamma (dParm2 - dTheta + 1);

      if (dParm1 == 0) { /* revised FB 18/07/97 */
        if (dTheta != 0)
          return (NULL_SUPPORT); /* should be -INF */
      }
      else
        dTmp = dTmp + dTheta * log(dParm1);

      if (dParm1 == 1) { /* revised FB 18/07/97 */
        if ((dParm2 - dTheta) == 0)
          return dTmp; /* because log(0^0) is 0 */
        else
          return (NULL_SUPPORT); /* should be -INF */
      }
      else
        return dTmp + (dParm2 - dTheta) * log(1 - dParm1);

    case MCV_PIECEWISE:
      density = 2 / (dMax + dParm2 - dParm1 - dMin);

      if (dTheta <= dParm1)
        return log(density * (dTheta - dMin) / (dParm1 - dMin));

      else
        if (dTheta <= dParm2)
          return log(density);
        else
          return log(density * (dMax - dTheta) / (dMax - dParm2));

    case MCV_EXPONENTIAL:
      if (dParm1 <= 0) {
        printf ("Error: negative or null inverse scale (%g) for exponential "
                "variate in LnDensity\n", dParm1);
        exit (0);
      }
      return -dTheta * dParm1 + log(dParm1);

    case MCV_GGAMMA:
      if (dParm2 <= 0) {
        printf ("Error: bad inv. scale for gamma variate in LnDensity\n");
        exit (0);
      }
      return (dParm1 - 1) * log(dTheta) - dParm2 * dTheta +
             dParm1 * log(dParm2) - lnGamma (dParm1);

    case MCV_TRUNCINVGGAMMA:
#ifdef HAVE_LIBGSL
      dTmp = gsl_cdf_gamma_P (1/dMax, dParm1, dParm2) -
             gsl_cdf_gamma_P (1/dMin, dParm1, dParm2); /* fall thru */
#else
      printf ("Error: Truncated inverse gamma density cannot be evaluated\n");
      printf ("       if the GNU Scientific Library is not installed\n");
      exit (0);
#endif
    case MCV_INVGGAMMA:
      if (dParm2 <= 0) {
        printf ("Error: bad scale for inv. gamma variate in LnDensity\n");
        exit (0);
      }
      return (-dParm1 - 1) * log(dTheta) - dParm2 / dTheta +
             dParm1 * log(dParm2) - lnGamma (dParm1) - dTmp;

    case MCV_POISSON:
      if (dParm1 <= 0) {
        printf ("Error: bad rate for Poisson variate in LnDensity\n");
        exit (0);
      }
      return dTheta * log(dParm1) - dParm1 - lnGamma (dTheta + 1);

    case MCV_BINOMIALBETA:
      if (dParm1 < 0) {
        printf ("Error: bad expectation for BinomialBeta variate "
                "in LnDensity\n");
        exit (0);
      }
      if (dParm2 <= 0) {
        printf ("Error: bad alpha for BinomialBeta variate in LnDensity\n");
        exit (0);
      }
      if (dMin <= 0) {
        printf ("Error: bad beta for BinomialBeta variate in LnDensity\n");
        exit (0);
      }
      dTmp = floor (0.5 + dParm1 + dParm1 * dMin / dParm2); /* this is N */
      if (dTheta > dTmp)
        return (NULL_SUPPORT);
      else {
        if ((dParm2 == 0.5) && (dMin == 0.5))
          dTmp = lnGamma (0.5 + dTheta) +
                 lnGamma (0.5 + dTmp - dTheta) -
                 lnGamma (dTheta + 1) - lnGamma (dTmp - dTheta + 1);
        else
          dTmp = lnGamma (dParm2 + dMin) + lnGamma (dTmp + 1) +
                 lnGamma (dParm2 + dTheta) +
                 lnGamma (dMin + dTmp - dTheta) -
                 lnGamma (dTheta + 1) - lnGamma (dTmp - dTheta + 1) -
                 lnGamma (dParm2) - lnGamma(dMin) -
                 lnGamma (dParm2 + dMin + dTmp);
        return dTmp;
      }

    case MCV_GENLOGNORMAL: /* Generalized LogNormal */
      if (dParm1 < 0) {
        printf ("Error: bad expectation for GenLogNormal variate "
                "in LnDensity\n");
        exit (0);
      }
      /* This is relative stdev of lognormal part -- dMin is sigma for the
         lognormal part */
      dTmp = sqrt(exp(pow(dMin,2)) * (exp(pow(dMin,2)) - 1));
      /* This is the transformation parameter lambda */
      dTmp2 = pow(dParm2/dTmp,2);
      /* This is the transformed mean */
      dTmp3 = log(dParm1 + sqrt(pow(dParm1,2) + dTmp2));
      /* This is the transformed theta -- use Taylor series for
         negative dTheta when lambda/dTheta^2 < 0.01 */
      if ((dTheta < 0) && (dTmp2 < (0.01*dTheta*dTheta)))
        dTmp4 = log(dTmp2/(-2*dTheta)*(1+0.25*dTmp2/pow(dTheta,2)));
      else
        dTmp4 = log(dTheta + sqrt(pow(dTheta,2) + dTmp2));
      /* Normal log-density minus log-jacobian */
      return lnDFNormal( dTmp4, dTmp3, dTmp ) - 0.5*log(pow(dTmp4,2) + dTmp2);

  case MCV_STUDENTT: /* Student t, dParm1 is dof, dParm2 is m, dMin is sigma */
      if (dParm1 <= 0) {
        printf ("Error: bad dof for Student-T variate"
                "in LnDensity\n");
        exit(0);
      }
      dTmp = (dParm1 + 1)/ 2;
      return (lnGamma(dTmp) - lnGamma(dParm1 / 2)
              -0.5 * log(dParm1 * PI *dMin *dMin)
              -dTmp * log(1 + pow((dTheta - dParm2) / dMin, 2) / dParm1));

    case MCV_CAUCHY:
      return (LNINVERPI - log(dParm1 + dTheta * dTheta / dParm1));

    case MCV_HALFCAUCHY: /* dParm1 is scale */
      return (LN2OVERPI - log(dParm1 + dTheta * dTheta / dParm1));

    default:
      ReportRunTimeError(panal, RE_UNKNOWNDIST | RE_FATAL, "LnDensity");

  } /* switch */

  /* Not reached */
  return 0.0 ;

} /* LnDensity */


/* Function -------------------------------------------------------------------
   LnLike

   returns the log-likelihood of the dependents given the parameter specified
   by pMCVar
*/
double LnLike (PMCVAR pMCVar, PANALYSIS panal)
{
  long n;
  double dDensity, dLnLike = 0.0;

  for (n = 0; n < pMCVar->nDependents; n++) {
    dDensity = LnDensity(pMCVar->rgpDependents[n], panal);
    if (dDensity == NULL_SUPPORT)
      return NULL_SUPPORT;
    else
      dLnLike += dDensity;
  }
  return dLnLike;

} /* LnLike */


/* Function -------------------------------------------------------------------
   LnLikeData

   Likelihood of the data for one experiment
*/
double LnLikeData (PLEVEL plevel, PANALYSIS panal) {

  PMCVAR pMCVar;
  long   i, j, k;
  double dLnLike = 0.0;
  double dTmp;
  BOOL   bMissData, bMissOutp;
  static PDOUBLE pdBase[4];

  /* For all the likelihoods seen at the experiment level */
  for (i = 0; i < plevel->nLikes; i++) {
    pMCVar = plevel->rgpLikes[i];

    for (k = 0; k < 4; k++)
      pdBase[k] = pMCVar->pdParm[k]; /* store the base pointers of pdParm */

    /* Set dVal of pMCVar to the current data value */
    for (j = 0; j < pMCVar->lCount; j++) { /* for each data value */

      pMCVar->dVal = pMCVar->pdVal[j]; /* point to jth data value */

      /* If the main data is coded as missing do nothing, and otherwise: */
      if (pMCVar->dVal != INPUT_MISSING_VALUE) {
        /* Set the pdParms of pMCVar to point to data or output values */
        bMissData = FALSE; /* initialize */
        bMissOutp = FALSE; /* initialize */
        for (k = 0; k < 4; k++) {
          if (pMCVar->iParmType[k] == MCVP_PRED) {
            /* Advance in the prediction array */
            pMCVar->pdParm[k] = pdBase[k] + j;
            bMissOutp = bMissOutp + (*(pMCVar->pdParm[k]) == MISSING_VALUE);
          }
          else if (pMCVar->iParmType[k] == MCVP_DATA) {
            /* Advance in the data array */
            pMCVar->pdParm[k] = pdBase[k] + j;
            bMissData = bMissData +
                        (*(pMCVar->pdParm[k]) == INPUT_MISSING_VALUE);
          }
        } /* end for k */

        /* If no missing data among the k parameters */
        if (bMissData == FALSE) {
          /* If no missing model output among the k parameters */
          if (bMissOutp == FALSE) {
            dTmp = LnDensity (pMCVar, panal);
            if (dTmp == NULL_SUPPORT) {
              /* Reset the pdParms to the beginning of arrays  FB 9/6/99 */
              for (k = 0; k < 4; k++)
                pMCVar->pdParm[k] = pdBase[k];
              return (NULL_SUPPORT);
            }
            else
              dLnLike = dLnLike + dTmp;
          }
          else /* missing output, cannot compute */
            ReportRunTimeError (panal, RE_BADMODEL | RE_FATAL, "LnLikeData");
        } /* end if bMissData */

      } /* end if ! INPUT_MISSING_VALUE */
    } /* end for j */

    /* Reset the pdParms to the beginning of the output or data arrays */
    for (k = 0; k < 4; k++)
      pMCVar->pdParm[k] = pdBase[k];

  } /* end for i */

  return (dLnLike);

} /* LnLikeData */


/* Function -------------------------------------------------------------------
   MaxMCVar

   get the maximum of MC variables
*/
double MaxMCVar (PMCVAR pMCVar) {

  /* if the parameter has a discrete distribution, round it - FB 12/06/97 */
  if (pMCVar->iType == MCV_BINOMIAL || pMCVar->iType == MCV_POISSON ) {
   return( *(pMCVar->pdParm[3]));
  }
  else { /* FB fixed the uniform case - FB 30/06/97 */
    if (pMCVar->iType == MCV_UNIFORM || pMCVar->iType == MCV_LOGUNIFORM ) {
      return( *(pMCVar->pdParm[1]));
    }
    else {
      return( *(pMCVar->pdParm[3]));
    }
  }

} /* MaxMCVar*/


/* Function -------------------------------------------------------------------
   MinMCVar

    get the minimum of MC variables
*/
double MinMCVar (PMCVAR pMCVar) {

  /* if the parameter has a discrete distribution, round it - FB 12/06/97 */
  if (pMCVar->iType == MCV_BINOMIAL || pMCVar->iType == MCV_POISSON) {
   return( *(pMCVar->pdParm[2]));
  }
  else { /* FB fixed the uniform case - FB 30/06/97 */
    if (pMCVar->iType == MCV_UNIFORM || pMCVar->iType == MCV_LOGUNIFORM) {
      return( *(pMCVar->pdParm[0]));
    }
    else {
      return( *(pMCVar->pdParm[2]));
    }
  }

} /* MinMCVar*/


/* Function -------------------------------------------------------------------
   OpenMarkovFiles

   Opens the output file and the restart file.
*/
void OpenMarkovFiles (PANALYSIS panal)
{
  /* Take care of the output file first */

  /* Use command line spec if given */
  if (panal->bCommandLineSpec) {
    free (panal->gd.szGout);
    panal->bAllocatedFileName = FALSE;
    panal->gd.szGout = panal->szOutfilename;
  }

  /* Default if none given */
  else if (!(panal->gd.szGout))
    panal->gd.szGout = "MCMC.default.out";

  /* eventually open the restart file before crushing the output file */
  if (panal->gd.szGrestart)
    if (!(panal->gd.pfileRestart)
      && !(panal->gd.pfileRestart = fopen (panal->gd.szGrestart, "r")))
      ReportRunTimeError(panal, RE_FATAL | RE_CANNOTOPEN,
                         panal->gd.szGrestart, "OpenMarkovFiles");

  if (!(panal->gd.pfileOut)
      && !(panal->gd.pfileOut = fopen (panal->gd.szGout, "w")))
    ReportRunTimeError(panal, RE_FATAL | RE_CANNOTOPEN,
                       panal->gd.szGout, "OpenMarkovFiles");

} /* OpenMarkovFiles */


/* Function -------------------------------------------------------------------
   PrintAllExpts

   Run and print all experiments at level below `plevel'
*/
void PrintAllExpts (PLEVEL plevel, PANALYSIS panal, PFILE pOutFile)
{
  long n;

  for (n = 0; n < plevel->iInstances; n++)
    TraverseLevels1 (plevel->pLevels[n], PrintExpt, panal, pOutFile, NULL);

} /* PrintAllExpts */


/* Function -------------------------------------------------------------------
   PrintExpt

   Run all experiments and print the level code, experiment number, time, data,
   and predictions to the output file.

   Called from TraverseLevels1
*/
int PrintExpt (PLEVEL plevel, char **args)
{
  PANALYSIS   panal = (PANALYSIS)args[0];
  PFILE       pOutFile = (PFILE)args[1];
  long        k, l, m, n;
  PEXPERIMENT pExpt = plevel->pexpt;
  POUTSPEC    pos;
  static long printed_head = 0;

  if (!printed_head) {
    fprintf (pOutFile,
             "Level\tSimulation\tOutput_Var\tTime\tData\tPrediction\n");
    printed_head = 1;
  }

  /* Set level sequence */
  panal->pCurrentLevel[plevel->iDepth] = plevel;

  panal->iInstance[plevel->iDepth] = plevel->iSequence;

  if (pExpt != NULL) {
    InitModel ();

    /* Set the model vars that have been sampled in this experiment and
       above levels */
    for (n = 0; n <= plevel->iDepth; n++) {
      SetModelVars (panal->pCurrentLevel[n]);
      SetFixedVars (panal->pCurrentLevel[n]);
    }

    if (!DoOneExperiment (pExpt)) {
      /* Error */
      printf ("Warning: DoOneExperiment failed\n");
      return 0;
    }
    else {
      pos = &pExpt->os;
      for (m = 0; m < pos->nOutputs; m++) {
        /* find the corresponding data index in case Print and Data are not
           in the same order */
        for (k = 0; k < pos->nData; k++)
          if (!strcmp(pos->pszDataNames[k], pos->pszOutputNames[m]))
            break;

        for (l = 0; l < pos->pcOutputTimes[m]; l++) {
          for (n = 1; n < plevel->iDepth; n++)
            fprintf (pOutFile, "%d_", panal->iInstance[n]);
          fprintf (pOutFile, "%d\t", panal->iInstance[plevel->iDepth]);

          if (k != pos->nData) /* Data found */
            fprintf (pOutFile, "%d\t%s\t%g\t%g\t%g\n", pExpt->iExp,
                     pos->pszOutputNames[m], pos->prgdOutputTimes[m][l],
                     pos->prgdDataVals[k][l], pos->prgdOutputVals[m][l]);
          else /* data not found, empty field */
            fprintf (pOutFile, "%d\t%s\t%g\t\t%g\n", pExpt->iExp,
                     pos->pszOutputNames[m], pos->prgdOutputTimes[m][l],
                     pos->prgdOutputVals[m][l]);
        } /* end for l */
        fprintf (pOutFile, "\n");

      } /* end for m */
      fprintf (pOutFile, "\n");

    } /* end else */

  } /* end if pExpt */

  return (1);

} /* PrintExpt*/


/* Function -------------------------------------------------------------------
   ReadData

   Reads data for an experiment. There need to be one data item
   per output specified, and in the order given by the Print or PrintStep
   statements.
   Upgraded by FB 27/03/99 to deal with improved data handling.
   Called from TraverseLevels.
*/
void ReadData (PLEVEL plevel, char **args)
{
  FILE *pfileData = (FILE *) args[0];
  POUTSPEC pos;
  int cDat, i, j;

  if (plevel->pexpt == NULL)
    return;

  pos = &(plevel->pexpt->os);

  /* here the number of data must equal the number of outputs */
  cDat = pos->nOutputs;
  pos->prgdDataVals = InitpdVector (cDat);
  pos->pcData       = InitiVector (cDat);  /* FB 5/11/99 */
  pos->pszDataNames = (PSTR *) malloc (cDat * sizeof(PSTR));
  pos->phvar_dat    = (HVAR *) malloc (cDat * sizeof(HVAR));

  if (pos->prgdDataVals == NULL || pos->phvar_dat == NULL ||
      pos->pszDataNames == NULL || pos->pcData    == NULL)
    ReportError (NULL, RE_OUTOFMEM | RE_FATAL, "ReadData()", NULL);
  else {
    pos->nData = cDat;  /* FB 27/03/99: Set count of data */

    /* scan all data values for data file */
    for (i = 0; i < cDat; i++) {
      /* allocate space for pos->prgdDataVals[i] */
      if ( !(pos->prgdDataVals[i] = InitdVector (pos->pcOutputTimes[i])))
        ReportError (NULL, RE_OUTOFMEM | RE_FATAL, "ReadData()", NULL);;

      /* for each requested output time */
      for (j = 0; j < pos->pcOutputTimes[i]; j++)
        if (fscanf(pfileData, "%lg", &(pos->prgdDataVals[i][j])) == EOF) {
          printf ("Error: incorrect length for data file - Exiting\n");
          exit(0);
        }
      pos->pcData[i] = j; /* FB 5/11/99 */
      pos->phvar_dat[i] = pos->phvar_out[i];
      pos->pszDataNames[i] = pos->pszOutputNames[i];

    } /* for i */

  } /* else */

} /* ReadData */


/* Function -------------------------------------------------------------------
   ReadDataFile

   Reads in data from a specified data file, if specified.
*/
void ReadDataFile (PANALYSIS panal)
{
  if (panal->gd.szGdata) {

    char c;

    FILE *pfile = fopen (panal->gd.szGdata, "r");

    if (!pfile) {
      printf ("Cannot open data file '%s'\n", panal->gd.szGdata);
      exit (0);
    }

    /* skip the first line */
    do { c = getc(pfile); } while (c != '\n');

    TraverseLevels (panal->pLevels[0], ReadData, pfile, NULL);

    fclose (pfile);
  }

} /* ReadDataFile */


/* Function -------------------------------------------------------------------
   PrintDeps

   Called from TraverseLevels

   For debugging, print the variables, parents, and dependencies
*/
void PrintDeps (PLEVEL plevel, char **args)
{
  long n, m;
  PMCVAR pMCVar;

  printf ("Depth %d; Instance %d\n", plevel->iDepth, plevel->iSequence);

  for (n = 0; n < plevel->nMCVars; n++) {
    pMCVar = plevel->rgpMCVars[n];

    printf ("Variable %s (%d) [%" PRIxPTR "]\n",
            pMCVar->pszName, pMCVar->iDepth, (intptr_t) pMCVar);

    for (m = 0; m < 4; m++)
      if (pMCVar->pMCVParent[m] != NULL)
        printf ("  Parent %ld: %s (%d) [%" PRIxPTR "]\n", m,
                pMCVar->pMCVParent[m]->pszName, pMCVar->pMCVParent[m]->iDepth,
                (intptr_t) pMCVar->pMCVParent[m]);

    for (m = 0; m < pMCVar->nDependents; m++)
      printf ("  Dependent: %s (%d) [%" PRIxPTR "]\n",
              pMCVar->rgpDependents[m]->pszName,
              pMCVar->rgpDependents[m]->iDepth,
              (intptr_t) pMCVar->rgpDependents[m]);

    if (pMCVar->bExptIsDep)
      printf("  This variable influences experiments directly\n");
  }

} /* PrintDeps */


/* Function -------------------------------------------------------------------
   ReadRestart

   initialize the parameters by reading them in the restart file.
*/
void ReadRestart (FILE *pfileRestart, long nThetas,
                  double *pdTheta, double *pdSum, double **prgdSumProd,
                  long *pIter)
{
  register char c;
  register long i, j;

  *pIter = -1;

  for (i = 0; i < nThetas; i++) {
    pdSum[i] = 0.0;
    for (j = 0; j < nThetas; j++)
      prgdSumProd[i][j] = 0.0;
  }

  /* skip the first line. This allows a MC output file to be used
     directly as a restart file. */
  do { c = getc(pfileRestart); } while (c != '\n');

  /* as long as we have not reached the end of the file we keep reading lines
     and overwriting the thetas, they keep only their last value.
     We throw away first field, and we keep incrementing the global iteration
     counter iter:
   */
  while (!(feof (pfileRestart) ||
          (fscanf (pfileRestart, "%*s") == EOF))) {
    for (i = 0; i < nThetas; i++) {
      if (fscanf(pfileRestart, "%lg", &(pdTheta[i])) == EOF) {
        printf ("Error: incorrect length for restart file - Exiting\n");
        exit(0);
      }
      else { /* reading ok */
        /* update pdSum, the column sum of pdTheta */
        pdSum[i] = pdSum[i] + pdTheta[i];
      }
    }

    /* Throw away remainder of line. This allows a MC output file to be used
       directly as a restart file. */
    do { c = getc(pfileRestart); } while (c != '\n');

    /* update prgdSumProd */
    for (i = 0; i < nThetas; i++)
      for (j = 0; j < nThetas; j++)
        prgdSumProd[i][j] = prgdSumProd[i][j] + pdTheta[i] * pdTheta[j];


    /* increment pIter */
    *pIter = *pIter + 1;

  } /* end while */

  /* note that the theta returned is the last parameter set read */

  fclose (pfileRestart);

} /* ReadRestart */


/* Function -------------------------------------------------------------------
   ReadRestartTemper

   initialize parameters and temperature stuff (pseudoprior, indexT) by reading
   them in the restart file.
*/
void ReadRestartTemper (FILE *pfileRestart, long nThetas, int nInvTemperatures,
                        double *pdTheta, double *pdSum, double **prgdSumProd,
                        long *pIter, long *pindexT, double *pdlnPi)
{
  register char c;
  register long i, j;
  double *pdAux;
  long iT;
  *pIter = -1;

  pdAux = InitdVector (nThetas);

  for (i = 0; i < nThetas; i++) {
    pdSum[i] = 0.0;
    for (j = 0; j < nThetas; j++)
      prgdSumProd[i][j] = 0.0;
  }

  /* skip the first line. This allows a MC output file to be used
     directly as a restart file. */
  do { c = getc(pfileRestart); } while (c != '\n');

  /* as long as we have not reached the end of the file we keep reading lines
     and overwriting the thetas, they keep only their last value.
     We throw away first field, and we keep incrementing the global iteration
     counter iter: */
  while (!(feof (pfileRestart) ||
           (fscanf (pfileRestart, "%*s") == EOF))) {
    for (i = 0; i < nThetas; i++) {
      if (fscanf(pfileRestart, "%lg", &(pdTheta[i])) == EOF) {
        printf ("Error: incorrect length for restart file - Exiting\n");
        exit(0);
      }
      else { /* reading ok */
        /* update pdSum, the column sum of pdTheta */
        pdSum[i] = pdSum[i] + pdTheta[i];
      }
    }

    if (fscanf(pfileRestart,"%ld", &(iT)) == EOF) {
      printf ("Error: incorrect length for restart file - Exiting\n");
      exit(0);
    }

    for (i = 0; i < nInvTemperatures; i++) {
      if (fscanf(pfileRestart,"%lg", &(pdAux[i])) == EOF) {
        printf ("Error: incorrect length for restart file - Exiting\n");
        exit(0);
      }
    }

    /* Throw away remainder of line. This allows a MC output file to be used
       directly as a restart file. */
    do { c = getc(pfileRestart); } while (c != '\n');

    *pindexT = iT;
    for (i = 0; i < nInvTemperatures; i++)
      pdlnPi[i] = pdAux[i];

    /* update prgdSumProd */
    for (i = 0; i < nThetas; i++)
      for (j = 0; j < nThetas; j++)
        prgdSumProd[i][j] = prgdSumProd[i][j] + pdTheta[i] * pdTheta[j];

    /* increment pIter */
    *pIter = *pIter + 1;

  } /* end while */

  /* note that the theta returned is the last parameter set read */

  fclose (pfileRestart);
  free (pdAux);

} /* ReadRestartTemper */


/* Function -------------------------------------------------------------------
   RestoreLikelihoods

   Called from TraverseLevels1
*/
int RestoreLikelihoods (PLEVEL plevel, char **args)
{
  PEXPERIMENT pExpt = plevel->pexpt;

  if (pExpt != NULL) {
    pExpt->dLnLike = pExpt->dLnLikeSave;
  }

  return (1);

} /* RestoreLikelihoods */


/* Function -------------------------------------------------------------------
   RunAllExpts

   Run all experiments (i.e. experiments at all levels below
   panal->pLevels[0]).
*/
int RunAllExpts (PANALYSIS panal, PDOUBLE pdLnData)
{
  PLEVEL plevel0 = panal->pLevels[0];
  long n;

  for (n = 0; n < plevel0->iInstances; n++) {
    if (!TraverseLevels1 (plevel0->pLevels[n], RunExpt, panal,
                          pdLnData, NULL)) {
      /* error */
      return (0);
    }
  }

  /* success */
  return (1);

} /* RunAllExpts */


/* Function -------------------------------------------------------------------
   RunExpt

   If `plevel' has experiments, modify the variables that have been sampled
   at this level and above using its two lists, and run the experiments.
   Return 1 on success and 0 if failure.
   Computes the data likelihood in case of success.

   Called from TraverseLevels1
*/
int RunExpt (PLEVEL plevel, char **args)
{
  PANALYSIS   panal = (PANALYSIS) args[0];
  double      *pdLnData = (double *) args[1];
  long        i;
  PEXPERIMENT pExpt = plevel->pexpt;

  /* Set level sequence */
  panal->pCurrentLevel[plevel->iDepth] = plevel;

  if (pExpt != NULL) {
    InitModel ();

    /* Set the model vars that have been sampled in this level and above */
    for (i = 0; i <= plevel->iDepth; i++) {
      SetModelVars (panal->pCurrentLevel[i]);
      SetFixedVars (panal->pCurrentLevel[i]);
    }

    if (!DoOneExperiment (pExpt)) {
      /* Error */
      printf ("Warning: DoOneExperiment failed\n");
      return 0;
    }
    else {
      pExpt->dLnLike = LnLikeData (plevel, panal);
      *pdLnData = (*pdLnData) + pExpt->dLnLike;
    }
  } /* if */

  return (1);

} /* RunExpt */


/* Function -------------------------------------------------------------------
   SampleTemperature
*/
long SampleTemperature (PGIBBSDATA pgd, double dLnPrior, double dLnData)
{
  int    indexT = pgd->indexT;
  int    indexT_new;

  /* Propose a new inverse temperature */
  if (indexT == 0) indexT_new = 1;
  else {
    if (indexT == pgd->nInvTemperatures - 1) indexT_new = indexT - 1;
    else {
      if (Randoms() > 0.5) indexT_new = indexT + 1;
        else indexT_new = indexT - 1;
    }
  }

  /* Test the temperature */
  if (TestTemper (pgd, indexT, indexT_new, dLnPrior, dLnData,
                  pgd->rgdlnPi[indexT], pgd->rgdlnPi[indexT_new])) {
    /* jump */
    return (indexT_new);
  }
  else
    return (indexT);

}  /* SampleTemperature */


/* Function -------------------------------------------------------------------
   SampleTheta, samples from a normal kernel
*/
double SampleTheta (PMCVAR pMCVar) {

  /* if the parameter has a discrete distribution, round it - FB 12/06/97 */
  if (pMCVar->iType == MCV_BINOMIAL || pMCVar->iType == MCV_POISSON) {
   return floor (0.5 + TruncNormalRandom (pMCVar->dVal, pMCVar->dKernelSD,
                                          MinMCVar(pMCVar),
                                          MaxMCVar(pMCVar)));
  }
  else { /* FB fixed the uniform case - FB 30/06/97 */
    return TruncNormalRandom (pMCVar->dVal, pMCVar->dKernelSD,
                              MinMCVar(pMCVar), MaxMCVar(pMCVar));
  }

} /* SampleTheta */


/* Function -------------------------------------------------------------------
   SampleThetaUnif, samples from a uniform kernel
*/
double SampleThetaUnif (PMCVAR pMCVar) {

  /* if the parameter has a discrete distribution, round it */
  if (pMCVar->iType == MCV_BINOMIAL || pMCVar->iType == MCV_POISSON)
    return floor (0.5 + UniformRandom (MinMCVar(pMCVar), MaxMCVar(pMCVar)));
  else
    return UniformRandom (MinMCVar(pMCVar), MaxMCVar(pMCVar));

} /* SampleThetaUnif */


/* Function -------------------------------------------------------------------
   SampleThetas

   Perform a Metropolis step on all the random (MC) variables at the
   level passed as argument 1.
   Sample thetas in sequence - test using prior and likelihood -
   restore old values if necessary - write new values to output file.

   Called from TraverseLevels
*/
void SampleThetas (PLEVEL plevel, char **args)
{
  PANALYSIS  panal = (PANALYSIS) args[0];
  PGIBBSDATA pgd = (PGIBBSDATA) args[1];
  long *pnIter = (long *) args[2];
  long *pnUpdateAt = (long *) args[3];
  long *pnTotal = (long *) args[4];

  double dLnPrior, dLnLike, dLnData, dLnKern;
  double dLnPriorNew, dLnLikeNew, dLnDataNew, dLnKernNew;
  double dTheta, dJumps;
  PMCVAR pMCVar;
  long   n;

  /* Set level sequence. This is needed to call RunExpt from the middle
     of the tree with the right initialization of the nodes above */
  panal->pCurrentLevel[plevel->iDepth] = plevel;

  /* For all MC vars at this level */
  for (n = 0; n < plevel->nMCVars; n++) {

    pMCVar = plevel->rgpMCVars[n];

    /* If the MC var is fixed, no sampling is made, just write it out */
    if (pMCVar->bIsFixed)
      goto WriteIt;

    /* Compute prior and likelihood */
    dLnPrior = LnDensity (pMCVar, panal);
    dLnLike  = LnLike (pMCVar, panal);

    dLnData = 0.0;

    /* If data are dependent compute the data likelihood */
    if (pMCVar->bExptIsDep) {
      /* Form the likelihood of all experiments at this level or beneath.*/
      TraverseLevels1 (plevel, SumAllExpts, &dLnData, NULL);
    }

    /* Save current value */
    dTheta = pMCVar->dVal;

    /* Adjust the jumping kernel SD, depending on acceptance rates,
       make sure it does not exceed DBL_MAX or a third of the range */
    if (*pnIter == *pnUpdateAt) {

      dJumps = (double) pMCVar->lJumps / (double) (*pnTotal);

      if (dJumps > 0.3) { /* we will increase the kernel spread */
        if (dJumps == 1) {
          /* pathological case: chances are the parameter has such a
             small kernel that it does not change values, this leads
             to a "jump" each time. Drastic increase of the kernel is
             attempted */
          if (pMCVar->dKernelSD < sqrt(DBL_MAX)) {
            if (pMCVar->dKernelSD > 2)
              pMCVar->dKernelSD = pMCVar->dKernelSD * pMCVar->dKernelSD;
            else /* FB 28/03/1999 */
              pMCVar->dKernelSD = pMCVar->dKernelSD * 20;
          }
          else
            pMCVar->dKernelSD = DBL_MAX;
        }
        else { /* more normal case */
          if (pMCVar->dKernelSD < DBL_MAX / 2)
            pMCVar->dKernelSD = pMCVar->dKernelSD * 2;
          else
            pMCVar->dKernelSD = DBL_MAX;
        }

        /* check that kernel SD does not increase wildly */
        if (pMCVar->dKernelSD > pMCVar->dMaxKernelSD)
          pMCVar->dKernelSD = pMCVar->dMaxKernelSD;
      }
      else { /* we will decrease the kernel spread */
        if (dJumps == 0) {
          /* pathological case: chances are the parameter has such a
             big kernel that it will never jump. Drastic decrease of
             the kernel is attempted, dissymetric from the increase */
          if (pMCVar->dKernelSD > pow(DBL_MIN, 0.45)) {
            if (pMCVar->dKernelSD > 2) { /* FB 28/03/1999 */
              pMCVar->dKernelSD = pow(pMCVar->dKernelSD, 0.45);
            }
            else {
              pMCVar->dKernelSD = pMCVar->dKernelSD * 0.04;
            }
          }
          else
            pMCVar->dKernelSD = DBL_MIN;
        }
        else { /* more normal case */
          /* the kernel should not be infinitely decreased; decrease
             is otherwise slighly different from the increase to avoid
             oscillations */
          if (pMCVar->dKernelSD > DBL_MIN / 0.4)
            pMCVar->dKernelSD = pMCVar->dKernelSD * 0.4;
          else
            pMCVar->dKernelSD = DBL_MIN;
        }
      }

      pMCVar->lJumps = 0; /* reset the jumps counter */

    } /* end of kernel stuff */

    /* compute the integral of the jump kernel */
    dLnKern  = log(CDFNormal ((MaxMCVar(pMCVar) - pMCVar->dVal) /
                              pMCVar->dKernelSD) -
                   CDFNormal ((MinMCVar(pMCVar) - pMCVar->dVal) /
                              pMCVar->dKernelSD));

    /* sample a new value */
    pMCVar->dVal = SampleTheta (pMCVar);

    /* update the integral of the jump kernel */
    dLnKernNew  =  log(CDFNormal ((MaxMCVar(pMCVar) - pMCVar->dVal) /
                                  pMCVar->dKernelSD) -
                       CDFNormal ((MinMCVar(pMCVar) - pMCVar->dVal) /
                                  pMCVar->dKernelSD));

    /* update prior and likelihood */
    dLnPriorNew = LnDensity (pMCVar, panal);
    dLnLikeNew  = LnLike (pMCVar, panal);

    dLnDataNew  = 0.0;

    /* If data are dependent compute the data likelihood */
    if (pMCVar->bExptIsDep) {
      /* Run all experiments at this level or beneath.
         We should in fact run only the dependent experiments ! */

      if (!TraverseLevels1 (plevel, RunExpt, panal, &dLnDataNew, NULL)) {
        /* If running experiments fails, do not jump */
        pMCVar->dVal = dTheta;
        TraverseLevels1 (plevel, RestoreLikelihoods, NULL);
        goto WriteIt;
      }
    }

    /* Test the results and act accordingly */
    if (!TestImpRatio (pgd, pMCVar->bExptIsDep, dLnKern, dLnKernNew,
                       dLnPrior, dLnPriorNew,
                       dLnLike, dLnLikeNew, dLnData, dLnDataNew)) {
      pMCVar->dVal = dTheta;

      if(pMCVar->bExptIsDep)
        TraverseLevels1 (plevel, RestoreLikelihoods, NULL);
    }
    else {
      pMCVar->lJumps = pMCVar->lJumps + 1;

      if(pMCVar->bExptIsDep)
        TraverseLevels1 (plevel, SaveLikelihoods, NULL);
    }

    WriteIt: /* Write the MC var value to output file */

    if (((*pnIter+1) % pgd->nPrintFreq == 0) &&
        (*pnIter >= pgd->nMaxIter - pgd->nPrintIter)) {
      fprintf(pgd->pfileOut, "%5g\t", pMCVar->dVal);
    }
  } /* end for all pMCVar */

} /* SampleThetas */


/* Function -------------------------------------------------------------------
   SampleThetasTempered

   Use annealing MCMC
   Sample thetas in sequence - test using prior and likelihood -
   restore old values if necessary - write new values to output file
   Called from TraverseLevels

*/
void SampleThetasTempered (PLEVEL plevel, char **args)
{
  PANALYSIS  panal = (PANALYSIS) args[0];
  PGIBBSDATA pgd = (PGIBBSDATA) args[1];
  long *pnIter = (long *) args[2];
  long *pnUpdateAt = (long *) args[3];
  long *pnTotal = (long *) args[4];
  long *pindexT = (long *) args[5];

  double dLnPrior, dLnLike, dLnData, dLnKern;
  double dLnPriorNew, dLnLikeNew, dLnDataNew, dLnKernNew;
  double dTheta, dJumps, old_dKernelSD;
  PMCVAR pMCVar;
  long   n;

  /* Set level sequence. This is needed to call RunExpt from the middle
     of the tree with the right initialization of the nodes above */
  panal->pCurrentLevel[plevel->iDepth] = plevel;

  /* For all MC vars at this level */
  for (n = 0; n < plevel->nMCVars; n++) {

    pMCVar = plevel->rgpMCVars[n];

    /* If the MC var is fixed, no sampling is made, just write it out */
    if (pMCVar->bIsFixed)
      goto WriteIt;

    /* Compute prior and likelihood */
    dLnPrior = LnDensity (pMCVar, panal);
    dLnLike = LnLike (pMCVar, panal);

    dLnData = 0.0;

    /* If data are dependent compute the data likelihood */
    if (pMCVar->bExptIsDep) {
      /* Form the likelihood of all experiments at this level or beneath.*/
      TraverseLevels1 (plevel, SumAllExpts, &dLnData, NULL);
    }

    /* Save current value */
    dTheta = pMCVar->dVal;

    /* Adjust the jumping kernel SD, depending on acceptance rates,
       make sure it does not exceed DBL_MAX or a third of the range */
    if (*pnIter == *pnUpdateAt) {

      dJumps = (double) pMCVar->lJumps / (double) (*pnTotal);

      if (dJumps > 0.3) { /* we will increase the kernel spread */
        if (dJumps == 1) {
          /* pathological case: chances are the parameter has such a
             small kernel that it does not change values, this leads
             to a "jump" each time. Drastic increase of the kernel is
             attempted */
          if (pMCVar->dKernelSD < sqrt(DBL_MAX)) {
            if (pMCVar->dKernelSD > 2)
              pMCVar->dKernelSD = pMCVar->dKernelSD * pMCVar->dKernelSD;
            else /* FB 28/03/1999 */
              pMCVar->dKernelSD = pMCVar->dKernelSD * 20;
          }
          else
            pMCVar->dKernelSD = DBL_MAX;
        }
        else { /* more normal case */
          if (pMCVar->dKernelSD < DBL_MAX / 2)
            pMCVar->dKernelSD = pMCVar->dKernelSD * 2;
          else
            pMCVar->dKernelSD = DBL_MAX;
        }

        /* check that kernel SD does not increase wildly */
        if (pMCVar->dKernelSD > pMCVar->dMaxKernelSD)
          pMCVar->dKernelSD = pMCVar->dMaxKernelSD;
      }
      else { /* we will decrease the kernel spread */
        if (dJumps == 0) {
          /* pathological case: chances are the parameter has such a
             big kernel that it will never jump. Drastic decrease of
             the kernel is attempted, dissymetric from the increase */
          if (pMCVar->dKernelSD > pow(DBL_MIN, 0.45)) {
            if (pMCVar->dKernelSD > 2) { /* FB 28/03/1999 */
              pMCVar->dKernelSD = pow(pMCVar->dKernelSD, 0.45);
            }
            else {
              pMCVar->dKernelSD = pMCVar->dKernelSD * 0.04;
            }
          }
          else
            pMCVar->dKernelSD = DBL_MIN;
        }
        else { /* more normal case */
          /* the kernel should not be infinitely decreased; decrease
             is otherwise slighly different from the increase to avoid
             oscillations */
          if (pMCVar->dKernelSD > DBL_MIN / 0.4)
            pMCVar->dKernelSD = pMCVar->dKernelSD * 0.4;
          else
            pMCVar->dKernelSD = DBL_MIN;
        }
      }

      pMCVar->lJumps = 0; /* reset the jumps counter */

    } /* end kernel updating stuff */

    /* sample a new value */

    if (pgd->rgInvTemperatures[*pindexT] > 0) { /* usual case */

      /* first scale temporarily the kernelSD according to the temperature */

      /* save the current kernelSD */
      old_dKernelSD = pMCVar->dKernelSD;

      /* update the kernelSD according to the temperature, the formula
         comes from that of the variance of the powered normal */
      pMCVar->dKernelSD = pow(2 * PI,
                              0.5 * (1 - pgd->rgInvTemperatures[*pindexT])) *
                          pow(pMCVar->dKernelSD,
                              3 - pgd->rgInvTemperatures[*pindexT]) /
                          pow(pgd->rgInvTemperatures[*pindexT], 1.5);

      /* check that kernel SD does not increase wildly */
      if (pMCVar->dKernelSD > pMCVar->dMaxKernelSD)
        pMCVar->dKernelSD = pMCVar->dMaxKernelSD;

      /* compute the integral of the jump kernel */
      dLnKern  = log(CDFNormal ((MaxMCVar(pMCVar) - pMCVar->dVal) /
                                pMCVar->dKernelSD) -
                     CDFNormal ((MinMCVar(pMCVar) - pMCVar->dVal) /
                                pMCVar->dKernelSD));

      /* sample a new value */
      pMCVar->dVal = SampleTheta (pMCVar);

      /* update the integral of the jump kernel */
      dLnKernNew  =  log(CDFNormal ((MaxMCVar(pMCVar) - pMCVar->dVal) /
                                    pMCVar->dKernelSD) -
                         CDFNormal ((MinMCVar(pMCVar) - pMCVar->dVal) /
                                    pMCVar->dKernelSD));
    }
    else { /* inverse temperature is zero, sample uniformly or from the prior */
      dLnKern = dLnKernNew  = 1;
      if (pgd->nSimTypeFlag == 3) { /* the posterior is tempered */
        /* sample a new value uniformly in its range */
        pMCVar->dVal = SampleThetaUnif (pMCVar);
      }
      else { /* the likelihood is tempered, sample from the prior */
        CalculateOneMCParm (pMCVar);
      }
    }

    /* recompute prior and likelihood */
    dLnPriorNew = LnDensity (pMCVar, panal);
    dLnLikeNew  = LnLike (pMCVar, panal);

    dLnDataNew = 0.0;

    /* If data are dependent compute the data likelihood */
    if (pMCVar->bExptIsDep) {
      /* Run all experiments beneath this level.
         We should in fact run only the dependent experiments ! */
      if (!TraverseLevels1 (plevel, RunExpt, panal, &dLnDataNew, NULL)) {
        /* If running experiments fails, do not jump */
        pMCVar->dVal = dTheta;
        TraverseLevels1 (plevel, RestoreLikelihoods, NULL);
        goto WriteIt;
      }
    }

    /* Test the results and act accordingly */
    if (!TestImpRatioTemper (pgd, pMCVar->bExptIsDep, dLnKern, dLnKernNew,
                             dLnPrior, dLnPriorNew,
                             dLnLike, dLnLikeNew, dLnData, dLnDataNew,
                             *pindexT)) {
      pMCVar->dVal = dTheta;
       if(pMCVar->bExptIsDep)
        TraverseLevels1 (plevel, RestoreLikelihoods, NULL);
    }
    else {
      pMCVar->lJumps = pMCVar->lJumps + 1;
      if(pMCVar->bExptIsDep)
        TraverseLevels1 (plevel, SaveLikelihoods, NULL);
    }

    if (pgd->rgInvTemperatures[*pindexT] > 0) /* restore kernelSD */
      pMCVar->dKernelSD = old_dKernelSD;

WriteIt: /* Write the MC var value to output file  */

    if (((*pnIter+1) % pgd->nPrintFreq == 0) &&
        (*pnIter >= pgd->nMaxIter - pgd->nPrintIter)) {
      fprintf(pgd->pfileOut, "%5g\t", pMCVar->dVal);
    }
  }  /*  (Metrolopolis-Hastings */

} /* SampleThetasTempered */


/* Function -------------------------------------------------------------------
   SampleThetaVector

   Sample thetas in block, do Metropolis test, write the sampled values to
   the output file
*/
void SampleThetaVector (PLEVEL pLevel, PANALYSIS panal, long nThetas,
                        double *pdTheta, double *pdSum, double **prgdSumProd,
                        long iter, long nUpdateAt, long nTotal,
                        PDOUBLE pdLnPrior, PDOUBLE pdLnData)
{
  register long i, j;
  double dTmp, dAccept, dLnPrior_old, dLnData_old;
  BOOL   bInBounds;

  static long lAccepted = 0;
  static double dJumpSpread;
  static PDOUBLE pdTheta_old = NULL; /* previous model parameters values */
  static PDOUBLE *prgdComponent;
  static PDOUBLE *prgdVariance;
  static PDOUBLE dNormVar; /* storage for nParms normal deviates */

  if ((pdTheta_old == NULL) || (iter == nUpdateAt)) {

    if (pdTheta_old == NULL) { /* initialize */
      if ( !(pdTheta_old = InitdVector (nThetas)) ||
           !(dNormVar    = InitdVector (nThetas)) ||
           !(prgdVariance  = InitdMatrix (nThetas, nThetas)) ||
           !(prgdComponent = InitdMatrix (nThetas, nThetas)))
        ReportRunTimeError (panal, RE_OUTOFMEM | RE_FATAL,
                            "SampleThetaVector");

      /* initialize dJumpSpread */
      dJumpSpread = 2.4 / sqrt(nThetas); /* Gelman's Normal theory result */
    }
    else {
      /* done if iter == nUpdateAt, but not at start:
         if some vector samplings have been made, check that the current
         dJumpSpread leads to an acceptation rate of 15 to 30% over the
         last batch of simulations. Adjust eventually. */
      dAccept = ((double) lAccepted) / (double) (nTotal);

      if ( dAccept > 0.3) dJumpSpread = dJumpSpread * 1.5;
      else if (dAccept < 0.15) dJumpSpread = dJumpSpread * 0.7;

      printf ("Monitoring: iter\t%ld\t", iter);
      printf ("success rate\t%g\tspread\t%g\n", dAccept, dJumpSpread);
      lAccepted = 0; /* reset the counter */
    }

    /* other updates: */

    /* generate the covariance matrix */
    for (i = 0; i < nThetas; i++)
      for (j = 0; j < nThetas; j++)
        prgdVariance[i][j] = (prgdSumProd[i][j] -
                              pdSum[i] * pdSum[j] / (double) (iter+1)) /
                             (double) iter;

    /* do the Cholesky decomposition of the covariance matrix */
    if (!Cholesky (prgdVariance, prgdComponent, nThetas)) {
      /* try to save the computation by zeroing all non-diagonal elements */
      for (i = 0; i < nThetas; i++)
        for (j = 0; j < nThetas; j++) {
          if (i == j)
            prgdVariance[i][j] = prgdSumProd[i][j] / (double) (iter);
          else
            prgdVariance[i][j] = 0.0;
        }

        /* if it still does not work, exit */
      if (!Cholesky (prgdVariance, prgdComponent, nThetas)) {
        printf ("Error: impossible to compute a jumping kernel - Exiting."
                "(You should check or change the restart file).\n\n");
        exit (0);
      }
    }
  }

  /* keep the value of all thetas */
  for (i = 0; i < nThetas; i++)
    pdTheta_old[i] = pdTheta[i];

  /* keep old prior and old likelihood */
  dLnPrior_old = *pdLnPrior;
  dLnData_old  = *pdLnData;

  /* generate new pdTheta vector */
  for (i = 0; i < nThetas; i++)
    dNormVar[i] = NormalRandom(0.0, 1.0);

  for (i = 0; i < nThetas; i++) {
    dTmp = 0;
    for (j = 0; j <= i; j++) /* only the non-zero part of prgdComponent */
      dTmp = dTmp + dNormVar[j] * prgdComponent[i][j];

    pdTheta[i] = pdTheta_old[i] + dJumpSpread * dTmp;
  }

  /* Set the dVals of the variables to the values sampled and check that
     we are within bounds */
  long iTmp = 0; /* dummy variable to avoid resetting nThetas FB 16/07/2016 */
  bInBounds = TraverseLevels1 (pLevel, SetMCVars, pdTheta, &iTmp, NULL);

  if (!bInBounds) { /* reject */
    for (i = 0; i < nThetas; i++)
      pdTheta[i] = pdTheta_old[i]; /* restore */
    goto DontJump;
  }

  /* Calculate the new prior */
  *pdLnPrior = 0.0;
  TraverseLevels (pLevel, CalculateTotals, panal, pdLnPrior, NULL);

  /* compute the model at the newly drawn point and the likelihood */
  *pdLnData = 0.0;
  if (!RunAllExpts (panal, pdLnData)) {
    /* computation failed, don't jump, restore */
    for (i = 0; i < nThetas; i++)
      pdTheta[i] = pdTheta_old[i];
    *pdLnPrior = dLnPrior_old;
    *pdLnData  = dLnData_old;
  }
  else {
    /* Test */
    if (log(Randoms()) >
        ((*pdLnPrior) + (*pdLnData) - dLnPrior_old - dLnData_old)) {
      /* don't jump, restore */
      for (i = 0; i < nThetas; i++)
        pdTheta[i] = pdTheta_old[i];
      *pdLnPrior = dLnPrior_old;
      *pdLnData  = dLnData_old;
    }
    else { /* jump */
      lAccepted++; /* this is used above to adjust the acceptation rate */
    }
  }

  DontJump:

  /* update arrays */
  for (i = 0; i < nThetas; i++) {
    pdSum[i] = pdSum[i] + pdTheta[i];
    for (j = 0; j < nThetas; j++)
      prgdSumProd[i][j] = prgdSumProd[i][j] + pdTheta[i] * pdTheta[j];
  }

} /* SampleThetaVector */


/* Function -------------------------------------------------------------------
   SaveLikelihoods

   Called from TraverseLevels1
*/
int SaveLikelihoods (PLEVEL plevel, char **args)
{
  PEXPERIMENT pExpt = plevel->pexpt;

  if (pExpt != NULL) {
    pExpt->dLnLikeSave = pExpt->dLnLike;
  }

  return (1);

} /* SaveLikelihoods */


/* Function -------------------------------------------------------------------
   SetFixedVars

   Set the array of fixed variables
*/
void SetFixedVars (PLEVEL plevel)
{
  long n;
  PVARMOD pFVar;

  for (n = 0; n < plevel->nFixedVars; n++) {
    pFVar = plevel->rgpFixedVars[n];
    if (IsInput (pFVar->hvar))
      SetInput (pFVar->hvar, pFVar->uvar.pifn);
    else
      SetVar (pFVar->hvar, pFVar->uvar.dVal);
  }

} /* SetFixedVars */


/* Function -------------------------------------------------------------------
   SetKernel

   Set initial values of the MCMC jumping kernel
*/
void SetKernel (PLEVEL plevel, char **args)
{
  intptr_t useMCVarVals = (intptr_t) args[0]; /* 1 to restore dVal, else 2 */
  double *pdMCVarVals = (double *) args[1];
  double dMin, dMax, dTmp;
  long   n, m;
  static long nThetas;
  PMCVAR pMCVar;

  /* set the jumping kernel's SD: sample 4 variates and take the range */
  for (n = 0; n < plevel->nMCVars; n++) {

    if (!(plevel->rgpMCVars[n]->bIsFixed)) {

      pMCVar = plevel->rgpMCVars[n];
      CalculateOneMCParm (pMCVar);

      /* set the maximum kernel SD value */
      if (pMCVar->iType == MCV_UNIFORM || pMCVar->iType == MCV_LOGUNIFORM )
        pMCVar->dMaxKernelSD = (*(pMCVar->pdParm[1]) - *(pMCVar->pdParm[0])) /
                               6.0;
      else
        pMCVar->dMaxKernelSD = (*(pMCVar->pdParm[3]) - *(pMCVar->pdParm[2])) /
                               6.0;

      /* sample 4 variates */
      dMin = dMax = pMCVar->dVal;
      for (m = 0; m < 3; m++) {
        CalculateOneMCParm (pMCVar);
        dTmp = pMCVar->dVal;
        if (dMin >= dTmp) dMin = dTmp;
        else if (dMax < dTmp) dMax = dTmp;
      }

      /* set the range safely */
      if ((*(pMCVar->pdParm[2]) == -DBL_MAX) ||
          (*(pMCVar->pdParm[3]) == DBL_MAX))
        pMCVar->dKernelSD = (0.5 * dMax) - (0.5 * dMin);
      else
        pMCVar->dKernelSD = dMax - dMin;

      /* take care of the case in which the range is zero
         (can happens for discrete variables) */
      if (pMCVar->dKernelSD == 0) {
        pMCVar->dKernelSD = pMCVar->dMaxKernelSD;
      }
    }

    /* restore the value of the variable - FB 02/07/97 */
    if (useMCVarVals == 1)
      plevel->rgpMCVars[n]->dVal = pdMCVarVals[nThetas++];

  }

} /* SetKernel */


/* Function -------------------------------------------------------------------
   WriteKernel

   Write values of the MCMC jumping kernel
*/
void WriteKernel (PLEVEL plevel, char **args)
{
  FILE   *pfile = (FILE *) args[0];
  long   n;
  PMCVAR pMCVar;

  for (n = 0; n < plevel->nMCVars; n++) {
    if ( !(plevel->rgpMCVars[n]->bIsFixed)) {
      pMCVar = plevel->rgpMCVars[n];
      fprintf(pfile, "%lg\t", pMCVar->dKernelSD);
    }
  }

} /* WriteKernel */


/* Function -------------------------------------------------------------------
   ReadKernel

   Read values of the MCMC jumping kernel
*/
void ReadKernel (PLEVEL plevel, char **args)
{
  FILE   *pfile = (FILE *) args[0];
  long   n;
  PMCVAR pMCVar;

  for (n = 0; n < plevel->nMCVars; n++) {
    if ( !(plevel->rgpMCVars[n]->bIsFixed)) {
      pMCVar = plevel->rgpMCVars[n];

      /* set the maximum kernel SD value */
      pMCVar->dMaxKernelSD = (MaxMCVar(pMCVar) - MinMCVar(pMCVar)) / 6.0;

      if (!fscanf(pfile, "%lg", &(pMCVar->dKernelSD))) {
        ReportError (NULL, RE_READERROR | RE_FATAL, "kernel file", NULL);
      }
    }
  }

} /* ReadKernel */


/* Function -------------------------------------------------------------------
   SetModelVars

   Sets the array of model variables to the sampled values. Does not set fixed
   variables. That has to be done by SetFixedVars.
*/
void SetModelVars(PLEVEL plevel)
{
  long n;
  PMCVAR  pMCVar;

  for (n = 0; n < plevel->nMCVars; n++) {
    pMCVar = plevel->rgpMCVars[n];
    if (!(pMCVar->bIsFixed) && (IsParm (pMCVar->hvar)))
      SetVar (pMCVar->hvar, pMCVar->dVal);
  }

} /* SetModelVars */


/* Function -------------------------------------------------------------------
   SetMCVars

   Set initial values of thetas after reading input file or sampling them.
   Values are assumed to be in proper order.
   It also checks that the ranges are respected.
   Called from TraverseLevels.
*/
int SetMCVars (PLEVEL plevel, char **args)
{
  double *pdMCVarVals = (double *) args[0];
  long   *nThetas = (long *) args[1];
  PMCVAR pMCVar;
  double dVar;
  long   n;

  for (n = 0; n < plevel->nMCVars; n++) {
    dVar = pdMCVarVals[(*nThetas)++];
    pMCVar = plevel->rgpMCVars[n];
    if ((pMCVar->iType == MCV_UNIFORM) || (pMCVar->iType == MCV_LOGUNIFORM)) {
      if ((dVar < *(pMCVar->pdParm[0])) || (dVar > *(pMCVar->pdParm[1]))) {
        /* error */
        return (0);
      }
    }
    else {
      if ((dVar < *(pMCVar->pdParm[2])) || (dVar > *(pMCVar->pdParm[3]))) {
        /* error */
        return (0);
      }
    }

    pMCVar->dVal = dVar;
  }

  /* success */
  return (1);

} /* SetMCVars */


/* Function -------------------------------------------------------------------
   SetPointers

   Called from TraverseLevels

   FB 26 nov 96: For each Monte Carlo variable, pdParms are set to point to the
   parent's dVals rather than to model parameters. If there is no parent,
   pdParms point to their own dParms.
   The pdVals and pdParms of the data and output arrays are also set.
*/
void SetPointers (PLEVEL plevel, char **args)
{
  long i, j, k;
  PMCVAR pMCVar;
  POUTSPEC pos;
  BOOL bFound;

  for (i = 0; i < plevel->nMCVars; i++) { /* for parameters */
    pMCVar = plevel->rgpMCVars[i];

    /* For each distribution parameter */
    for (j = 0; j < 4; j++) {
      if (pMCVar->pMCVParent[j] == NULL) /* Point to its own values */
        pMCVar->pdParm[j] = &(pMCVar->dParm[j]);
      else /* Point to the parent dVal */
        pMCVar->pdParm[j] = &(pMCVar->pMCVParent[j]->dVal);
    }
  }

  if (plevel->pexpt != NULL) /* if the level has experiments */
  for (i = 0; i < plevel->nLikes; i++) { /* for likelihoods */
    pMCVar = plevel->rgpLikes[i];
    pos = &(plevel->pexpt->os);

    /* Set pdVal of pMCVar to a data array, first find it */
    bFound = FALSE;
    j = 0;
    while ((j < pos->nData) && (!bFound)) {
      bFound = (pMCVar->hvar == pos->phvar_dat[j]);
      if (!bFound) j++;
    }
    if (bFound) {
      pMCVar->pdVal = pos->prgdDataVals[j];
      pMCVar->lCount = pos->pcData[j]; /* number of data points */
    }
    else { /* no corresponding Data found: error */
      printf ("Error: no Data for %s in Simulation %d - Exiting.\n\n",
              pMCVar->pszName, plevel->pexpt->iExp);
      exit (0);
    }

    /* we also have to set the pdParms of pMCVar to either data or output
       arrays */
    for (j = 0; j < 4; j++) {

      if (pMCVar->iParmType[j] == MCVP_PRED) { /* it's an output */
        /* set pdParm[j] to an output, first find it */
        bFound = FALSE;
        k = 0;
        while ((k < pos->nOutputs) && (!bFound)) {
          bFound = (pMCVar->hParm[j] == pos->phvar_out[k]);
          if (!bFound) k++;
        }
        if (bFound) {
          pMCVar->pdParm[j] = &(pos->prgdOutputVals[k][0]);
        }
        else { /* no corresponding Print statement found: error */
          printf ("Error: missing Print statement for parameter number %ld\n"
                  "of %s distribution - Exiting.\n\n", j, pMCVar->pszName);
          exit (0);
        }
      } /* end if MCVP_PRED */

      else if (pMCVar->iParmType[j] == MCVP_DATA) { /* it's data */
        /* set pdParm[j] to a data array, first find it */
        bFound = FALSE;
        k = 0;
        while ((k < pos->nData) && (!bFound)) {
          bFound = (pMCVar->hParm[j] == pos->phvar_dat[k]);
          if (!bFound) k++;
        }
        if (bFound) {
          pMCVar->pdParm[j] = &(pos->prgdDataVals[k][0]);
        }
        else { /* no Data found: error */
          printf ("Error: no Data for %s in Simulation %d - Exiting.\n\n",
                  pMCVar->pszName, plevel->pexpt->iExp);
          exit (0);
        }
      } /* end if MCVP_DATA */

      else { /* it's either a parameter or a numeric */
        if (pMCVar->pMCVParent[j] == NULL) /* Point to its own values */
          pMCVar->pdParm[j] = &(pMCVar->dParm[j]);
        else /* Point to the parent dVal */
          pMCVar->pdParm[j] = &(pMCVar->pMCVParent[j]->dVal);
      }
    } /* end for j */
  }

} /* SetPointers */


/* Function -------------------------------------------------------------------
   SumAllExpts

   If `plevel' has experiments, add the current Likelihood to the total
   passed as argument 2.

   Called from TraverseLevels1
*/
int SumAllExpts (PLEVEL plevel, char **args)
{
  double      *pdLnData = (double*)args[0];
  PEXPERIMENT pExpt = plevel->pexpt;

  if (pExpt != NULL) {
    *pdLnData += pExpt->dLnLike;
  }
  return (1);

} /* SumAllExpts */


/* Function -------------------------------------------------------------------
   TestImpRatio

   Test prior, likelihood against random number between 0 and 1
*/
BOOL TestImpRatio (PGIBBSDATA pgd,  BOOL bExptIsDep,
                   double dLnKern,  double dLnKernNew,
                   double dLnPrior, double dLnPriorNew,
                   double dLnLike,  double dLnLikeNew,
                   double dLnData,  double dLnDataNew)
{
  double dPjump;

  if (dLnKernNew  == NULL_SUPPORT ||
      dLnPriorNew == NULL_SUPPORT ||
      dLnLikeNew  == NULL_SUPPORT ||
      dLnDataNew  == NULL_SUPPORT)
    return FALSE;

  dPjump = dLnPriorNew - dLnPrior + dLnLikeNew - dLnLike +
           dLnKern - dLnKernNew;

  if (bExptIsDep)
    dPjump += dLnDataNew - dLnData;

  if (pgd->nSimTypeFlag == 0)
    return ((BOOL) (log(Randoms()) <= dPjump));
  else {
    if (pgd->nSimTypeFlag == 5)
      return ((BOOL) (0 <= dPjump));
    else {
      printf ("Error: simTypeFlag should be 0 or 5 in TestImpRatio "
              "- Exiting.\n\n");
      exit (0);
    }
  }

} /* TestImpRatio */


/* Function -------------------------------------------------------------------
   TestImpRatioTemper

   Test prior, likelihood against a random number between 0 and 1 according to
   the temperature
*/
BOOL TestImpRatioTemper (PGIBBSDATA pgd,  BOOL bExptIsDep,
                         double dLnKern,  double dLnKernNew,
                         double dLnPrior, double dLnPriorNew,
                         double dLnLike,  double dLnLikeNew,
                         double dLnData,  double dLnDataNew, long indexT)
{
  double dPjump;

  if (dLnPriorNew == NULL_SUPPORT ||
      dLnLikeNew  == NULL_SUPPORT ||
      dLnDataNew  == NULL_SUPPORT)
    return FALSE;

  if (pgd->nSimTypeFlag == 3) { /* posterior is tempered */
    dPjump = pgd->rgInvTemperatures[indexT] *
             (dLnPriorNew - dLnPrior + dLnLikeNew - dLnLike +
              dLnKern - dLnKernNew);
  }
  else { /* only the likelihood is tempered */
    dPjump = dLnPriorNew - dLnPrior +
             pgd->rgInvTemperatures[indexT] *
             (dLnLikeNew - dLnLike + dLnKern - dLnKernNew);
  }

  if (bExptIsDep)
    dPjump += pgd->rgInvTemperatures[indexT] * (dLnDataNew - dLnData);

  return ((BOOL) (log(Randoms()) <= dPjump));

} /* TestImpRatioTemper */


/* Function -------------------------------------------------------------------
   TestTemper

   Test temperature against a random number between 0 and 1
*/
BOOL TestTemper (PGIBBSDATA pgd, long indexT, long indexT_new, double dLnPrior,
                 double dLnData, double pseudo, double pseudonew)
{
  double dPjump = 0;
  #define MINUSLN2  -0.6931471805599452862268

  if (dLnPrior + dLnData == NULL_SUPPORT)
    return FALSE;

  if (pgd->nSimTypeFlag == 3) { /* the posterior is tempered */
    dPjump = (pgd->rgInvTemperatures[indexT_new] -
              pgd->rgInvTemperatures[indexT]) * (dLnPrior + dLnData) +
             pseudonew - pseudo +
             ((indexT_new == 0) || (indexT_new == pgd->nInvTemperatures - 1) ?
              0 : MINUSLN2) -
             ((indexT     == 0) || (indexT     == pgd->nInvTemperatures - 1) ?
              0 : MINUSLN2);
  }
  else { /* only the likelihood is tempered */
    dPjump = (pgd->rgInvTemperatures[indexT_new] -
              pgd->rgInvTemperatures[indexT]) * dLnData +
             pseudonew - pseudo +
             ((indexT_new == 0) || (indexT_new == pgd->nInvTemperatures - 1) ?
              0 : MINUSLN2) -
             ((indexT     == 0) || (indexT     == pgd->nInvTemperatures - 1) ?
              0 : MINUSLN2);
  }

  return ((BOOL) (log(Randoms()) <= dPjump));

} /* TestTemper */


/* Function -------------------------------------------------------------------
   TraverseLevels (recursive)

   Called with variable argument list ending in NULL;
   arguments should be pointers only; if you call this with a value
   that can be zero, you will be very sorry

   Find all allocated levels, execute `routinePtr' for each, starting at the
   top, passing the argument list as char**

   The argument list is saved from the initial call; on recursive calls the
   list is NULL
*/
void TraverseLevels (PLEVEL plevel,
                     void (*routinePtr)(PLEVEL plevel, char **args), ...)
{
  va_list ap;
  static char *arg[MAX_ARGS], **args = arg;
  char *arg1;
  long n, nargs = 0;

  va_start(ap, routinePtr);
  if ((arg1 = va_arg (ap, char*)) != NULL) {
    arg[0] = arg1;
    while ((arg[++nargs] = va_arg(ap, char*)) != NULL) {};
  }
  va_end (ap);

  routinePtr (plevel, args);

  for (n = 0; n < plevel->iInstances; n++)
    TraverseLevels (plevel->pLevels[n], routinePtr, NULL);

} /* TraverseLevels */


int TraverseLevels1 (PLEVEL plevel,
                     int (*routinePtr)(PLEVEL plevel, char **args), ...)
{
  va_list ap;
  static char *arg[MAX_ARGS], **args = arg;
  char *arg1;
  long n, nargs = 0;

  va_start (ap, routinePtr);
  if ((arg1 = va_arg (ap, char*)) != NULL) {
    arg[0] = arg1;
    while ((arg[++nargs] = va_arg(ap, char*)) != NULL) {};
  }
  va_end (ap);

  if (routinePtr (plevel, args)) {

    for (n = 0; n < plevel->iInstances; n++) {
      if (!TraverseLevels1(plevel->pLevels[n], routinePtr, NULL)) {
        /* error */
        return (0);
      }
    }
  }
  else /* error */
    return (0);

  /* success */
  return (1);

} /* TraverseLevels1 */


/* Function -------------------------------------------------------------------
   WriteHeader

   Called from Traverse Levels
   Write the names of the sampled parameters to output file header
*/
void WriteHeader (PLEVEL plevel, char **args)
{
  PANALYSIS panal = (PANALYSIS)args[0];
  FILE *outFile = (FILE*)args[1];
  long n, m;

  panal->iInstance[plevel->iDepth] = plevel->iSequence;

  for (n = 0; n < plevel->nMCVars; n++) {
    fprintf (outFile, "%s(", plevel->rgpMCVars[n]->pszName);
    for (m = 1; m < plevel->iDepth; m++)
      fprintf (outFile, "%d.", panal->iInstance[m]);
    fprintf (outFile, "%d)\t", panal->iInstance[plevel->iDepth]);
  }

} /* WriteHeader */


/* Function -------------------------------------------------------------------

   WriteMCVars

   Write the values of MC vars for one level to output file
*/
void WriteMCVars (PLEVEL plevel, PFILE pOutFile)
{
  long n;
  PMCVAR pMCVar;

  for (n = 0; n < plevel->nMCVars; n++) {
    pMCVar = plevel->rgpMCVars[n];
    fprintf(pOutFile, "%5g\t", pMCVar->dVal);
  }

} /* WriteMCVars */

/* End */
