/* sim.c

   Copyright (c) 1993-2017 Free Software Foundation, Inc.

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

   Entry point and main simulation routines for 'sim' program.

   This file can use the SUNDIALS CVODES libraries if they are installed.
   If so, the symbols HAVE_LIBSUNDIALS_CVODES and HAVE_LIBSUNDIALS_NVECSERIAL
   should be defined in config.h (or in the makefile).
   Otherwise, the corresponding features are disabled.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>

#include "delays.h"
#include "yourcode.h"
#include "getopt.h"
#include "mh.h"
#include "optdsign.h"
#include "lexerr.h"
#include "lsodes.h"
#include "sim.h"
#include "simi.h"
#include "siminit.h"
#include "simo.h"
#include "simmonte.h"
#include "strutil.h"
#include "config.h"

// CVODES specific includes and routines

#ifdef HAVE_LIBSUNDIALS_CVODES

#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */

#include <cvodes/cvodes.h>           /* prototypes CVODE fcts. and consts. */
#include <cvodes/cvodes_band.h>      /* prototype for CVBand */
#include <cvodes/cvodes_dense.h>     /* prototype for CVDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */
#include <sundials/sundials_math.h>  /* definition of ABS and EXP */


typedef struct {
  int nVars; // number of state and output variables
  int npes, my_pe;
} UserData;

/* Check function return value...
     opt == 0 means SUNDIALS function allocates memory so check if
              returned NULL pointer
     opt == 1 means SUNDIALS function returns a flag so check if
              flag >= 0
     opt == 2 means function allocates memory so check if returned
              NULL pointer */
static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */

  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  /* Check if flag < 0 */

  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */

  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  return(0);
}


/* ----------------------------------------------------------------------------
   f_for_cvodes

   stupid routine interfacing cvodes derivative function call and CalcDeriv.
*/
static int f_for_cvodes(realtype t, N_Vector u, N_Vector udot, void *user_data)
{
  static realtype *rgMVars = NULL;
  static int nStates, nVars;
  int i;

  /* problem with the output variables, derivatives have the proper length
     and do not need to be translated */
  if (rgMVars == NULL) { // initialize

    nStates = NV_LENGTH_S(udot);

    /* Extract number of outputs from user_data */
    nVars = ((UserData *) user_data)->nVars;

    // state and output vector
    rgMVars = (realtype *) malloc(nVars * sizeof(realtype));

    if (/*!dudata ||*/ !rgMVars)
      ReportError (NULL, RE_OUTOFMEM | RE_FATAL, "f_for_cvodes", NULL);
  }

  // copy u to rgMVars start
  for (i = 0; i < nStates; i++) {
    rgMVars[i] = NV_Ith_S(u, i);
  }

  CalcDeriv ((PDOUBLE) rgMVars, (PDOUBLE) NV_DATA_S(udot), (PDOUBLE) &t);

  return(0);

} /* f_for_cvodes */

#endif


/* ----------------------------------------------------------------------------
   CorrectInputToTransition

   resets the integrator and inputs when an input transition occurs.

   returns the simulation time pexp->dTime and input values to
   the input discontinuity, or transition point *pdTtrans.

   The inputs are updated to reflect their state just after the
   transition.  The integrator is initialized for a new segment.

   This does NOT affect state and output definitions.
*/
void CorrectInputToTransition (PEXPERIMENT pexp, PDOUBLE pdTtrans)
{
  pexp->dTime = *pdTtrans;
  UpdateInputs (&pexp->dTime, pdTtrans);

} /* CorrectInputToTransition */


/* ----------------------------------------------------------------------------
   Euler

   Simple Euler integrator.
*/
int Euler (long neq, double *y, double *t, double tout, double dTStep)
{
  static PDOUBLE rgdDeriv;
  double dTmp_step;
  long   i;

  if (!(rgdDeriv))
    if ( !(rgdDeriv = InitdVector (neq)))
      ReportError (NULL, RE_OUTOFMEM | RE_FATAL, "Euler", NULL);

  /* Iterate through time out to prescrebed output */
  while (*t < tout) {

    /* Compute derivatives at current time */
    CalcDeriv (y, rgdDeriv, t);

    /* Update time */
    *t = *t + dTStep;

    /* But do not exceed prescribed tout */
    if (*t > tout) {
      dTmp_step = tout - (*t - dTStep);
      *t = tout;
    }
    else
      dTmp_step = dTStep;

    /* Update the state variables */
    for (i = 0; i < neq; i++)
      y[i] = y[i] + dTmp_step * rgdDeriv[i];

    if (bDelays) 
      StoreDelayed(*t);

    DoStep_by_Step();
  }

  /* Calculate the derivatives at t = tout for cleanliness */
  CalcDeriv (y, rgdDeriv, t);

  return (0);

} /* Euler */


/* ----------------------------------------------------------------------------
   FreeVarMod

   Callback for FreeList().

   Frees the memory for one parameter modification.  If the parameter
   is an input, the memory allocated for the input function is also
   free'd.  Note that FreeList() will pass the data as a PVOID which
   needs to be re-cast.
*/
void FreeVarMod (PVOID pData)
{
  PVARMOD pvarmod = (PVARMOD) pData;

  if (IsInput (pvarmod->hvar))
    if (pvarmod->uvar.pifn) free (pvarmod->uvar.pifn);

  free (pvarmod);

} /* FreeVarMod */


/* ----------------------------------------------------------------------------
   ModifyOneParm

   Callback function for ModifyParms.
*/

int ModifyOneParm (PVOID pData, PVOID pNullInfo)
{
  PVARMOD pvarmod = (PVARMOD) pData;

  if (IsInput(pvarmod->hvar))
    SetInput (pvarmod->hvar, pvarmod->uvar.pifn);
  else
    SetVar (pvarmod->hvar, pvarmod->uvar.dVal);

  return 0;

} /* ModifyOneParm */


/* ----------------------------------------------------------------------------
   ModifyParms

   Modifies the parameters in the plistParmMods LIST of the experiment
   spec by call ForAllList to increment through the list.
*/
void ModifyParms (PLIST plistParmMods)
{

  assert (plistParmMods);
  ForAllList (plistParmMods, &ModifyOneParm, NULL);

} /* ModifyParms */


/* ----------------------------------------------------------------------------
   DoOneExperiment

   Runs one experiment - return 1 on success and 0 in case of errors
*/
int DoOneExperiment (PEXPERIMENT pexp)
{

  double dTout;     /* next output time */
  double dTtrans;   /* next exposure transition time */
  double dTup;      /* the smaller one of dTout or dTtrans*/
  int    iOut;      /* index to next output time */
  PMODELINFO pmod;  /* pointer to the current model info */
  PINTSPEC   pis;   /* pointer to the integrator specs */

#ifdef HAVE_LIBSUNDIALS_CVODES
  // CVODES specific variables
  static N_Vector u = NULL;
  static UserData user_data;
  static void *cvode_mem = NULL;
  int flag, i;
#endif

  if (!pexp) return 0;

  pmod = pexp->pmodelinfo;
  pis  = &(pexp->is);

  if (!InitOutputs (pexp, &iOut, &dTout)) return 0;

  /* Resolve dependency for initial time, eventually */
  if (pexp->hT0)
    pexp->dT0 = GetVarValue (pexp->hT0);

   /* Resolve dependent inputs, which calls ScaleModel */
  UpdateInputs (&pexp->dT0, &dTtrans);

  if (bDelays) 
    InitDelays(pexp->hT0);

  if (pexp->dT0 > dTtrans) {
    printf ("\nError: starting time is greater than first discontinuity,"
            "       check your inputs - Exiting.\n");
    exit (0);
  }

  if (pexp->dT0 > dTout) {
    printf ("\nError: starting time is greater than first output time,"
            "       check your outputs - Exiting.\n");
    exit (0);
  }

  pexp->dTime = pexp->dT0;

  // integrator initializations
  if (pis->iAlgo == IAL_LSODES) { /* Lsodes algorithm */
    /* set lsodes return flag to 1 for first call */
    pis->iDSFlag = 1;
  }
  else {
    if (pis->iAlgo == IAL_CVODES) { /* Sundials CVODE algorithm */

#ifdef HAVE_LIBSUNDIALS_CVODES

      if (1 || u == NULL) { /* always done for now */

        /* Create a serial vector for state variables */
        u = N_VNew_Serial(pmod->nStates);  /* Allocate u vector */
        if (check_flag((void*)u, "N_VNew_Serial", 0)) return(1);

        /* Call CVodeCreate to create the solver memory and specify the 
         * Backward Differentiation Formula and the use of a Newton iteration */
        cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
        if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

        /* set initial state values */
        for (i = 0; i < pmod->nStates; i++)
          NV_Ith_S(u, i) = pmod->pdModelVars[i];
    
        /* Call CVodeInit to initialize the integrator memory and specify the
         * user's right hand side function in u'=f(t,u), the inital time T0, and
         * the initial dependent variable vector u. */
        flag = CVodeInit(cvode_mem, f_for_cvodes, pexp->hT0, u);
        if (check_flag(&flag, "CVodeInit", 1)) return(1);

        /* Call CVodeSStolerances to specify the scalar relative tolerance
         * and scalar absolute tolerance */
        flag = CVodeSStolerances(cvode_mem, 
                                 RCONST(pis->dRtol), RCONST(pis->dAtol));
        if (check_flag(&flag, "CVodeSStolerances", 1)) return(1);

        /* Set the pointer to user-defined data used to store the total 
           number of state and output variables */
        user_data.nVars = pmod->nModelVars;
        flag = CVodeSetUserData(cvode_mem, &user_data);
        if (check_flag(&flag, "CVodeSetUserData", 1)) return(1);

        /* Call CVDense to specify the CVDENSE dense linear solver */
        flag = CVDense(cvode_mem, pmod->nStates);
        if (check_flag(&flag, "CVDense", 1)) return(1);
    
        /* Set the user-supplied Jacobian routine Jac, not used now
        flag = CVDlsSetBandJacFn(cvode_mem, Jac);
        if(check_flag(&flag, "CVDlsSetBandJacFn", 1)) return(1); */

      }
      else { /* disabled for now */
        /* reset initial state values */
        for (i = 0; i < pmod->nStates; i++)
          NV_Ith_S(u, i) = pmod->pdModelVars[i];

        flag = CVodeReInit(cvode_mem, pexp->dT0, u);
      }
#endif
    } /* end if IAL_CVODES */
  }

  /* Iterate to final time */
  while (pexp->dTime < pexp->dTfinal) {

    /* If dynamics equations are defined */
    if (pmod->nStates > 0) {

      /* the upper limit of integration dTup should be either dTout
         or dTtrans, whichever is smaller */
      /* F. Bois, 08 April 2007: before that, if the difference between dTout 
         and dTtrans is too small make dTtrans = dTout to avoid problems with
         the integration */
      if (fabs(dTout - dTtrans) < DBL_EPSILON * 2.0 * 
                                  mymax(fabs(dTout), fabs(dTtrans)))
        dTtrans = dTout;

      dTup = (dTout < dTtrans) ? dTout : dTtrans;

      /* F. Bois, 07 October 2008: same rounding problem fix: */
      if (fabs(dTup - pexp->dTime) < DBL_EPSILON * 2.0 * 
                                     mymax(fabs(dTup), fabs(pexp->dTime)))
        pexp->dTime = dTup;

      if (pis->iAlgo == IAL_LSODES) { /* Lsodes algorithm */

        pis->rwork[0] = dTup; /* do not overshoot dTup - FB 01/07/97 */

        lsodes_(&pmod->nStates, pmod->pdModelVars, &(pexp)->dTime,
                &dTup, &pis->itol, &pis->dRtol, &pis->dAtol,
                &pis->itask, &pis->iDSFlag, &pis->iopt, pis->rwork,
                &pis->lrw, pis->iwork, &pis->liw, &pis->iMf);

        /* Handle error returns : FB 25/11/96 : */
        if (pis->iDSFlag < 0) {
          /* We cannot guarantee the accuracy of the results, exit the routine
             with an error flag */
          return (0);
        }
      }
      else {
        if (pis->iAlgo == IAL_CVODES) { /* Sundials CVODE algorithm */

#ifdef HAVE_LIBSUNDIALS_CVODES
	  if (dTup > (pexp)->dTime) {
	    /* do not overshoot */
            flag = CVodeSetStopTime(cvode_mem, (realtype) dTup);
            flag = CVode(cvode_mem, dTup, u, &(pexp)->dTime, CV_NORMAL);
            if (check_flag(&flag, "CVode", 1))
              break;
            /* copy back state values */
            for (i = 0; i < pmod->nStates; i++) {
	      pmod->pdModelVars[i] = NV_Ith_S(u, i);
	    }
          }
#endif
        }
        else if (pis->iAlgo == IAL_EULER) { /* Euler algorithm */
          Euler(pmod->nStates, pmod->pdModelVars, &(pexp)->dTime, dTup,
                pis->dTStep);
        }
      }
    }
    else {
      /* We still need to advance the time */
      pexp->dTime = (dTout < dTtrans) ? dTout : dTtrans;
    }

    if (dTtrans <= dTout) {
      /* dTime == dTtrans <= dTout: we are at a discontinuity.
         This point belongs to the NEW period UNLESS we are at
         the final time */
      if (dTtrans < dTout) {
        if (dTtrans < pexp->dTfinal) {
          CorrectInputToTransition (pexp, &dTtrans);
          pis->iDSFlag = 1;
        }
      }
      else {
        /* dTtrans == dTout */
        if (dTtrans < pexp->dTfinal) {
          CorrectInputToTransition (pexp, &dTtrans);
          pis->iDSFlag = 1;
        }
        SaveOutputs (pexp, &dTout);
        NextOutputTime (pexp, &dTout, &iOut);
      }
    }
    else {
      /* dTime == dTout < dTtrans: */
      SaveOutputs (pexp, &dTout);
      NextOutputTime (pexp, &dTout, &iOut);
    }

  } /* while dTime < final time */

  if (pis->iAlgo == IAL_CVODES) { /* cleanup Sundials CVODE algorithm */
#ifdef HAVE_LIBSUNDIALS_CVODES
    /* Free vector u */
    N_VDestroy_Serial(u);    
    CVodeFree(&cvode_mem);  /* Free the integrator memory */
#endif
  }

  /* success */
  return 1;

} /* DoOneExperiment */


/* ----------------------------------------------------------------------------
   DoOneNormalExp

   Does one AT_DEFAULTSIM simulation.

   Return 1 on success and 0 in case of failure
*/
int DoOneNormalExp (PANALYSIS panal, PEXPERIMENT pexp)
{
  printf (" %d", pexp->iExp); /* Show what experiment it is */

  InitModel ();
  ModifyParms (panal->expGlobal.plistParmMods); /* Global modifications */
  ModifyParms (pexp->plistParmMods); /* Mods for this experiment */
  if (!DoOneExperiment (pexp)) {
    /* Error */
    return 0;
  }

  printf ("\n");

  return (1);

} /* DoOneNormalExp */


/* ----------------------------------------------------------------------------
   DoOneMCExp

   Does one AT_MONTECARLO simulation.

   Can maybe merge this with DoOneNormalExp() in the future.

   The major issue is the order of setting parameters.  For each
   experiment in a Monte Carlo run of an analysis, the order must be
   as follows:

   Each Run
    calc mc mods

     Each Simulation (old "Experiment")
     1)  Init the model
     2)  Global parm mods
     3)  Monte Carlo mods
     4)  Local mods override everything

   The problem becomes that for the simulation to be started over
   again, the inputs have to be told to initialize and parm mods for
   the current experiment must be made (body weight, etc).  This
   currently won't happen unless the model is init'd.  Maybe create a
   ResetInputs() with starting time which will do the funky stuff
   done by the global variable right now.

   Return 1 on success and 0 in case of failure
*/
int DoOneMCExp (PANALYSIS panal, PEXPERIMENT pexp)
{
  register MONTECARLO *pmc = &panal->mc;

  InitModel ();
  ModifyParms (panal->expGlobal.plistParmMods); /* Global modifications */
  SetParms (pmc->nParms, pmc->rghvar, pmc->rgdParms); /* MC mods */
  ModifyParms (pexp->plistParmMods); /* Mods for this experiment */
  if (!DoOneExperiment (pexp)) {
    /* Error */
    return 0;
  }

  return (1);

} /* DoOneMCExp */


/* ----------------------------------------------------------------------------
   DoNormal

   Does a normal analysis
*/
void DoNormal (PANALYSIS panal)
{
  int nExps = panal->expGlobal.iExp;
  int i;

  printf ("\nDoing analysis - %d normal experiment%c\n", nExps,
       (nExps > 1 ? 's' : ' '));

  for (i = 0; i < nExps; i++) {
    if (DoOneNormalExp (panal, panal->rgpExps[i])) {
      /* if successfull write out the results */
      WriteNormalOutput (panal, panal->rgpExps[i]);
    }
    else
      printf ("Warning: Integration failed - No output generated\n");
  }

} /* DoNormal */


/* ----------------------------------------------------------------------------
   DoMonteCarlo

   Does a Monte Carlo analysis or a Set Points analysis.  The latter is
   handled here because the looping is basically the same, with one
   difference.

   If the number of runs (nRuns) for SetPoints() analysis is
   specified to be zero, set points are read from the set points file
   until end of file is reached.  Otherwise, the number of runs
   explicity stated are read.  Not having enough points in the file
   in this latter case yields an error.

   If nRuns == 0, the test at the end of the while{} loop is not
   used, and the decision to continue is made by the return value of
   GetMCMods().  Since for MonteCarlo analyses, this return value is
   always TRUE (i.e. you can always pick another random number),
   nRuns is locally modified to 1, if it has been spec'd to zero,
   thus preventing the the Monte Carlo's from accidentaly running
   forever.
*/
void DoMonteCarlo (PANALYSIS panal)
{
  int nExps = panal->expGlobal.iExp;
  long nRuns = panal->mc.nRuns;
  MCPREDOUT mcpredout;
  BOOL bOK = FALSE, bNotDone; /* Not finished with analysis */
  int i;

  mcpredout.pred = NULL;

  if (panal->iType == AT_MONTECARLO && nRuns <= 0)
    nRuns = 1; /* Don't let MonteCarlo run forever */

  /* if cannot open files, Abort */
  if (OpenMCFiles (panal)) exit(0);

  printf ("\nDoing analysis - %ld %s run%c... %d experiment%c%s\n",
          nRuns,
          (panal->iType == AT_MONTECARLO ? "Monte Carlo" : "Set point"),
          (nRuns != 1 ? 's' : ' '),
          nExps, (nExps > 1 ? 's' : ' '),
          (nRuns != 1 ? " each" : " "));

  if (!nRuns)
    printf ("0 runs specified for SetPoint().  Reading entire file.\n\n");

  /* FB 21/07/97: dependencies between parameters are handled here. Pointers 
     to the parent parameters are stored in hparm of the pmcvar structure for
     each parameter. We now have to make pdParms point to the parent's dVals. */
  if (panal->iType == AT_MONTECARLO) {
    SetParents (&panal->mc, 0); /* start at 0, do them all */
  }
  else { /* panal->iType == AT_SETPOINTS */
    SetParents (&panal->mc, panal->mc.nSetParms);
  }

  panal->mc.lRun = 0; /* First run */
  bNotDone = TRUE;

  while (bNotDone) {

    bNotDone = GetMCMods (panal, NULL); /* Mods for this run */

    if (bNotDone) {
      /* Do analysis if not finished */
      for (i = 0; i < nExps; i++) { /* Do all experiments */
        bOK = DoOneMCExp (panal, panal->rgpExps[i]);
        if (!bOK) break;
      }

      if (bOK) {
        /* If successful write results */
        TransformPred (panal, &mcpredout); /* transform output run */
        WriteMCOutput (panal, &mcpredout);
      }
      else
        printf ("Warning: Integration failed on iteration %ld, experiment %d:\n"
                "         No output generated\n", panal->mc.lRun+1, i+1);
    } /* if bNotDone */

    panal->mc.lRun++; /* Next run */
    if (nRuns) /* If a number of runs spec'd... */
      bNotDone = (panal->mc.lRun < nRuns);

  } /* while */

  CloseMCFiles (panal);

  if (mcpredout.pred) free(mcpredout.pred);

} /* DoMonteCarlo */


/* ----------------------------------------------------------------------------
   DoAnalysis

   Does the analysis in the given specification.
*/
void DoAnalysis (PANALYSIS panal)
{

  InitRandom (panal->dSeed, TRUE);

  switch (panal->iType) {

    default:
    case AT_DEFAULTSIM:
      DoNormal (panal);
      break;

    case AT_SETPOINTS: /* Not really Monte Carlo */
    case AT_MONTECARLO:
      DoMonteCarlo (panal);
      break;

    case AT_MCMC:
      DoMarkov (panal);
      break;

    case AT_OPTDESIGN:
      DoOptimalDesign (panal);
      break;

  } /* switch */

  if (panal->pfileOut) {
    fclose (panal->pfileOut);
    printf ("Wrote output file \"%s\"\n", panal->szOutfilename);
  }

} /* DoAnalysis */


/* ----------------------------------------------------------------------------
   FreeMemory

   To use in the case of simulations without Levels.
*/
void FreeMemory (PANALYSIS panal)
{
  int i, j;

  free(panal->modelinfo.pStateHvar);

  FreeList (&panal->mc.plistMCVars, NULL, TRUE);
  if (panal->mc.rgdParms) {
    free (panal->mc.rgdParms);
    free (panal->mc.rghvar);
  }

  PINTSPEC pis = &panal->rgpExps[0]->is;
  free (pis->iwork);
  free (pis->rwork);

  for (i = 0; i < panal->expGlobal.iExp; i++) {
    if (panal->rgpExps[i] != NULL) {
      POUTSPEC pos = &panal->rgpExps[i]->os;

      FreeList (&panal->rgpExps[i]->plistParmMods, NULL, TRUE);  
      free (pos->pszOutputNames);
      free (pos->phvar_out);
      free (pos->pcOutputTimes);
      free (pos->piCurrentOut);
      free (pos->prgdOutputTimes);
      for (j = 0; j < pos->nOutputs; j++)
        free(pos->prgdOutputVals[j]);
      free (pos->prgdOutputVals);
      free (pos->rgdDistinctTimes);
      ForAllList (pos->plistPrintRecs, &FreePrintRec, NULL);
      FreeList (&pos->plistPrintRecs, NULL, FALSE);  
      free (pos->plistPrintRecs);
      ForAllList (pos->plistDataRecs, &FreeDataRec, NULL);
      FreeList (&pos->plistDataRecs, NULL, FALSE);  
      free (pos->plistDataRecs);
      free (panal->rgpExps[i]);
    }
  }
  if (panal->bAllocatedFileName) {
    if (panal->szOutfilename)           free (panal->szOutfilename);
    if (panal->mc.szMCOutfilename)      free (panal->mc.szMCOutfilename);
    if (panal->gd.szGout)               free (panal->gd.szGout);
  }
  
  if (panal->mc.szSetPointsFilename)  free (panal->mc.szSetPointsFilename);
  if (panal->gd.szGrestart)           free (panal->gd.szGrestart);
  if (panal->gd.szGdata)              free (panal->gd.szGdata);

  FreeList (&panal->expGlobal.plistParmMods, NULL, TRUE);
  free (panal);

} /* FreeMemory */


/* ----------------------------------------------------------------------------
   MCVarListToArray

   converts a list of MCVAR to an array.  This must be a callback for
   ForAllList() since we are making the change here that will let us
   not to be forced to use list traversal in the future.
*/

MCVAR **vrgpMCVar; /* Avoid hairy pointers in here */
int   viMCVar;     /* Index to the array */

int MCVarListToArray (PVOID pv_pMCVar, PVOID pv_Null)
{

  vrgpMCVar[viMCVar] = (MCVAR *) pv_pMCVar; /* Copy the pointer and.. */
  viMCVar++; /* Advance to next element of array */
  return 1;

} /* MCVarListToArray */


/* ----------------------------------------------------------------------------
   PrepAnalysis

   makes the ANALYSIS structure easier to work with in the simulation
   code. Specifically, changes lists to arrays.
*/
void PrepAnalysis (PANALYSIS panal)
{
  register MONTECARLO *pmc = &panal->mc;
  register int l;

  pmc->nParms = ListLength (pmc->plistMCVars);
  /* avoid zero pmc->nParms which can confuse some implementations of
     malloc. If pmc->nParms is zero  no use is going to be made of these
     arrays anyway */
  if (pmc->nParms == 0) return;
  
  pmc->rgdParms = InitdVector (pmc->nParms);
  pmc->rgpMCVar = (MCVAR **) malloc((pmc->nParms)*sizeof(MCVAR *));
  if (!(pmc->rgdParms && pmc->rgpMCVar))
    ReportError (NULL, RE_OUTOFMEM | RE_FATAL, "PrepAnalysis", NULL);

  /* Address of the first pointer */
  vrgpMCVar = &pmc->rgpMCVar[0];

  /* Initialize global array index */
  viMCVar = 0;
  ForAllList (pmc->plistMCVars, MCVarListToArray, (PVOID) NULL);
  FreeList (&pmc->plistMCVars, NULL, FALSE);

  /* Make a handle vector */
  pmc->rghvar = (HVAR *) malloc((pmc->nParms)*sizeof(HVAR));
  if (pmc->rghvar) {
    for (l = 0; l < pmc->nParms; l++)
      pmc->rghvar[l] = pmc->rgpMCVar[l]->hvar;
  }
  else
    ReportError (NULL, RE_OUTOFMEM | RE_FATAL, "PrepAnalysis", NULL);

} /* PrepAnalysis */


/* Get the command line argument stuff */

/* ----------------------------------------------------------------------------
   SansPath

   returns a pointer to just the filename of a full path.
*/
/*
char *SansPath (char *szFullPathname)
{
  register char *szFile;

  if ((szFile = szFullPathname))
    while (*szFullPathname) {
      if (*szFullPathname == '/')
        szFile = szFullPathname+1;
      szFullPathname++;
    }

  return szFile;

} */ /* SansPath */


/* ----------------------------------------------------------------------------
   PromptFilenames

   prompts for both input and output file names.  The space allocated
   for inputting the files is reallocated to their actual size.
*/
void PromptFilenames (PSTR *pszFileIn, PSTR *pszFileOut)
{
  *pszFileIn  = (PSTR) calloc (1, MAX_FILENAMESIZE);
  *pszFileOut = (PSTR) calloc (1, MAX_FILENAMESIZE);

  printf ("Input filename? ");

  if (!fgets (*pszFileIn, MAX_FILENAMESIZE, stdin)) {
    ReportError (NULL, RE_READERROR | RE_FATAL, "stdin", NULL);
  }
  else
    *pszFileIn = strtok (*pszFileIn, " \t\n");

  if (!(*pszFileIn)) /* Nothing entered, quit */
    return;

  if ((*pszFileIn)[0]) { /* Input file specified */
    printf ("Output filename? ");

    if (!fgets (*pszFileOut, MAX_FILENAMESIZE, stdin)) {
      ReportError (NULL, RE_READERROR | RE_FATAL, "stdin", NULL);
    }
    else
      *pszFileOut = strtok (*pszFileOut, " \t\n");
  }

  if (!(*pszFileOut) || !(*pszFileOut)[0]) { /* If no output specified */
    free (*pszFileOut);                      /* .. use default later */
    *pszFileOut = NULL;
  }
  else {
    *pszFileIn = (PSTR) realloc (*pszFileIn, MyStrlen(*pszFileIn) + 1);
    *pszFileOut = (PSTR) realloc (*pszFileOut, MyStrlen(*pszFileOut) + 1);
  }

} /* PromptFilenames */


/* ----------------------------------------------------------------------------
   Command Line Parsing

*/

/*
#define WarnArgumentReqd(szOption, szArg) \
  printf ("* Command-line option \"%s\" requires an argument \"%s\"\n",\
          szOption, szArg);

#define WarnUnknownArg(szOption, szArg) \
  printf ( "* Unknown argument \"%s\" to option \"%s\"\n", szArg, szOption);


void GetOutputFlagOption (PANALYSIS panal, char *optarg)
{
  WarnUnknownArg ("-O", optarg);

} */ /* GetOutputFlagOption */


/* ----------------------------------------------------------------------------
   GetCmdLineArgs

   retrieves options and filenames from the command line arguments passed to
   the program.

   The command line syntax is:

     mcsim [-options | --options] [input-file [output-file]]

   If the output filename is not given a default is used.
   If neither the input, nor output filenames are given, the
   program prompts for them both.

   The options can appear anywhere in the line and in any order.

   The options are parsed with _getopt(). After _getopt() is called,
   the args in rgszArg have been permuted so that non-option args are
   first, which in this case means the filenames.

   Uses the following globals:

     char *optarg;    -- Contains the string argument to each option in turn
     int   optind;    -- Index in ARGV of the next elem to be scanned
     char *nextchar;  -- The next char to be scanned in the option-element
     int   opterr;    -- 0 value flags to inhibit GNU error messages

*/

static char vszOptions[] = "h:H:D:";

void GetCmdLineArgs (int cArg, char *const *rgszArg, PSTR *pszFileIn, 
                     PSTR *pszFileOut, PANALYSIS panal)
{
  int c;

  *pszFileIn = *pszFileOut = (PSTR) NULL;

  while (1) {

    c = _getopt (cArg, rgszArg, vszOptions);
    if (c == EOF) /* Finish with option args */
      break;

    switch (c) {
      case '?':
        optarg = 0;
        /* Fall through! */

      case 'H':
      case 'h':
        /* if (optarg && *optarg) disabled for now, not up to date
          ShowHelp (optarg);
        else
          ShowHelpMessage (SansPath (rgszArg[0])); */
        exit (-1);
        break;

      case 'D':
        printf (">> Debug mode: Using option '%s'\n", optarg);
        /* Could setup to run with certain debug flags, not enabled  */
        break;

      default:
        printf ("Unknown option in command-line, %c = code 0%o ?\n", c, c);
        break;

    } /* switch */

  } /* while */

  switch (cArg - optind) { /* Remaining args are filenames */
    case 2: /* Output and input file specificed */
      *pszFileOut = rgszArg[optind + 1];

      /* Fall through! */

    case 1: /* Input file specificed */
      *pszFileIn = rgszArg[optind];
      break;

    case 0: /* No file names specified */
      PromptFilenames (pszFileIn, pszFileOut);
      break;

    default:
      /* ShowHelp ("Usage"); */ /* disabled for now, not updated */
      exit (-1);
      break;
  } /* switch */

  while (*pszFileIn && (*pszFileIn)[0] &&      /* Files specified   */
         !MyStrcmp(*pszFileIn, *pszFileOut)) { /* and not different */

    printf ("\n** Input and output filename must be different.\n");
    PromptFilenames (pszFileIn, pszFileOut);

  } /* while */

  if (!(*pszFileIn && (*pszFileIn)[0])) { /* no input name given is an error */
    printf ("Error: an input file name must be specified - Exiting\n\n");
    exit (-1);
  }

} /* GetCmdLineArgs */


/* ----------------------------------------------------------------------------
*/
void AnnounceProgram (void)
{
  printf ("\n________________________________________\n");
  printf ("\nMCSim " VSZ_VERSION "\n\n");
  printf (VSZ_COPYRIGHT "\n\n");

  printf ("MCSim comes with ABSOLUTELY NO WARRANTY;\n"
          "This is free software, and you are welcome to redistribute it\n"
          "under certain conditions; see the GNU General Public License.\n\n");

  printf ("* Using `%s' model in file \"%s\" created by %s\n\n",
          szModelDescFilename, szModelSourceFilename, szModelGenAndVersion);

} /* AnnounceProgram */


/* ----------------------------------------------------------------------------
   main

   Entry point for simulation and analysis program.
*/

int main (int nArg, char **rgszArg)
{
  PSTR szFileIn, szFileOut;
  INPUTBUF ibIn;
  PANALYSIS panal = (PANALYSIS) malloc (sizeof(ANALYSIS));

  AnnounceProgram ();

  if (!panal)
    ReportError (NULL, RE_OUTOFMEM | RE_FATAL,
                 "ANALYSIS specification too large", NULL);

  InitAnalysis (panal);
  GetCmdLineArgs (nArg, rgszArg, &szFileIn, &szFileOut, panal);

  /* Define the output file as the global experiment default  */
  panal->szOutfilename = szFileOut;
  szFileOut == NULL ? (panal->bCommandLineSpec = FALSE) :
                      (panal->bCommandLineSpec = TRUE);

  if (!InitBuffer (&ibIn, szFileIn))
    ReportError (&ibIn, RE_INIT | RE_FATAL, "ReadInput", NULL);

  ibIn.pInfo = (PVOID) panal; /* Attach analysis specification to input */

  if (ReadAnalysis (&ibIn)) {
    PrepAnalysis (panal);
    DoAnalysis (panal);
  }

  if (panal->iType == AT_MCMC)
    FreeLevels (panal);
  else {
    FreeMemory (panal);
    free (ibIn.pbufOrg);
  }

  return 0;

} /* main */
