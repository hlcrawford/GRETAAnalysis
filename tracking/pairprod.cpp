
#include <stdlib.h>
#include <fcntl.h>
#include <stddef.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/param.h>
#include <sys/resource.h>
#include <sys/times.h>
#include <assert.h>
#include <limits.h>

#include "ctk.h"

volatile extern TRACKINGPAR Pars;

/* npp / ppe / ... via GretaTrackingWorkspace macros in greta_tracking_std.hpp */

/*----------------------------------------------------------------------------*/

int
findpairprod2 (SHELLHIT *shellhit, STAT *ctkStat)
{

  /* second implemention where we search  */
  /* for flashes and look for 511s or 1022s  */
  /* in the vicinity */

  /*-------------------------------------------------*/
  /* NOTE:: used[] does not work as shellhit->used[] */
  /*-------------------------------------------------*/

  /* declarations */

  static int firsttime = 1;
  int i, j, k, l, m, n;

  int ii, jj, kk, ll, ninterest, ncol;
  float etot, etest, d1, d2;
  static FILE *fp;
  int nflash, nsto, ok;
  static int nprint;
  int pair_can[200], used[200];
  float flashx[200], flashy[200], flashz[200], flashe[200];
  long long int flashTS[200];
  float eexpect, efinal;
  unsigned long long int TStrace;

//  float emin = 1.022;           /* in MeV */
//  float de = 0.010;             /* in MeV */
//  float eflash = 3.0;           /* in MeV */
//  float maxdist = 5;            /* in cm */
//  float dflash = 2.0;           /* in cm */

  /* open log file? */

  if (firsttime)
    {
      fp = fopen ("pairprod.log", "w");
      firsttime = 0;
      nprint = 500;
      fprintf (fp, "\n");
      fprintf (fp, "processing pairprod.c with parameters\n");
      fprintf (fp, "\n");
      fprintf (fp, "Pars.ppfind=%i\n", Pars.ppfind);
      fprintf (fp, "Pars.ppemin=%f\n", Pars.ppemin);
      fprintf (fp, "Pars.ppde=%f\n", Pars.ppde);;
      fprintf (fp, "Pars.ppeflash=%f\n", Pars.ppeflash);;
      fprintf (fp, "Pars.ppmaxdist=%f\n", Pars.ppmaxdist);;
      fprintf (fp, "Pars.ppdflash=%f\n", Pars.ppdflash);;
      fprintf (fp, "\n");
    };

  ctkStat->ppstat[0]++;
  npp = 0;

//  fprintf (fp, "** event # %i >>\n", ctkStat->nTrackCalled);

  /* find total energy, any chance of pair prod? */

  etot = 0;
  for (jj = 0; jj < shellhit->nhit; jj++)
    etot += shellhit->edet[jj];
//  fprintf (fp, "etot= %f\n", etot);

  /* if too small, just quit out */
  /* this si based in non Doppler corrected energy */
  /* but probably somewhat OK */

  if (etot < Pars.ppemin)
    return (1);
  ctkStat->ppstat[1]++;

  /* find all flash interactions */

//  fprintf (fp, "\n");

  if (nprint > 0)
    {
      fprintf (fp, "\n");
      fprintf (fp, "--------------------\n");
    };

  nflash = 0;

  for (jj = 0; jj < shellhit->nhit; jj++)
    used[jj] = 0;
  ncol = 0;

  for (jj = 0; jj < shellhit->nhit; jj++)
    if (shellhit->edet[jj] > Pars.ppeflash)
      if (used[jj] == 0)
        {
//          fprintf (fp, "%i is flashpoint\n", jj);
          TStrace = shellhit->timestamp[jj];

          /* store */

          flashx[nflash] = shellhit->XX[jj];
          flashy[nflash] = shellhit->YY[jj];
          flashz[nflash] = shellhit->ZZ[jj];
          flashe[nflash] = shellhit->edet[jj];
          flashTS[nflash] = shellhit->timestamp[jj];

          shellhit->used[jj] = 1;
          used[jj] = 1;

          /* check if anyone else is nearby */
          /* and add them */

          for (ll = jj + 1; ll < shellhit->nhit; ll++)
            if (shellhit->edet[ll] > Pars.ppeflash)
              {
//                fprintf (fp,"check %i, other flash\n", ll);

                /* find distance to current */

                d1 = shellhit->XX[ll] - flashx[nflash];
                d2 = d1 * d1;
                d1 = shellhit->YY[ll] - flashy[nflash];
                d2 += d1 * d1;
                d1 = shellhit->ZZ[ll] - flashz[nflash];
                d2 += d1 * d1;
                d2 = sqrtf (d2);
//                fprintf (fp, "dist %i to %i is %f -- lim %f\n", jj, ll, d2, Pars.ppdflash);

                if (d2 < Pars.ppdflash)
                  {
//                    fprintf (fp, "dist %i to %i is %f -- lim %f MERGE\n", jj, ll, d2, Pars.ppdflash);
                    shellhit->used[ll] = 1;
                    used[ll] = 1;
                    ncol++;

                    /* merge this point, bary center */

                    d1 = flashe[nflash] + shellhit->edet[ll];
                    flashx[nflash] = (shellhit->edet[ll] * shellhit->XX[ll] + flashe[nflash] * flashx[nflash]) / d1;
                    flashy[nflash] = (shellhit->edet[ll] * shellhit->YY[ll] + flashe[nflash] * flashy[nflash]) / d1;
                    flashz[nflash] = (shellhit->edet[ll] * shellhit->ZZ[ll] + flashe[nflash] * flashz[nflash]) / d1;
                    flashe[nflash] = d1;

                  };

              };

          /* next flashpoint */

          nflash++;

        };

  /* if no flash points, we are done */

  if (nflash == 0)
    return (2);
  ctkStat->ppstat[2]++;

  /* if we get here, we have an event to process */

  ok = 1;
  if (nprint > 0)
    {
//      fprintf (fp, "\n");
//      fprintf (fp, "--------------------\n");
      fprintf (fp, "we have %i flash point(s) at event # %i TS~%lli >>\n", nflash, ctkStat->nTrackCalled, TStrace);
      fprintf (fp, "we consolidated %i flashpoints\n", ncol);
    };

  /* for each flash, find */
  /* nearby interaction points */

  if (nprint > 0)
    fprintf (fp, "nflash=%i\n", nflash);
  for (ll = 0; ll < nflash; ll++)
    {
      if (nprint > 0)
        fprintf (fp, "Flash point %i:: %6.2f %6.2f %6.2f [%6.3f] TS=%lli\n", ll, flashx[ll], flashy[ll], flashz[ll], flashe[ll], flashTS[ll]);

      /* look for interactions in the vicinity */
      /* of this flash that are not flash points */

      nsto = 0;
      for (kk = 0; kk < shellhit->nhit; kk++)
        if (shellhit->edet[kk] < Pars.ppeflash)
          if (kk != ll)
            {
              d1 = flashx[ll] - shellhit->XX[kk];
              d2 = d1 * d1;
              d1 = flashy[ll] - shellhit->YY[kk];
              d2 += d1 * d1;
              d1 = flashz[ll] - shellhit->ZZ[kk];
              d2 += d1 * d1;
              d2 = sqrtf (d2);

              if (d2 < Pars.ppmaxdist)
                {
                  /* store this candidate hit */
                  /* simply it's shell hit number */

                  pair_can[nsto] = kk;
                  nsto++;

                };
            };

      if (nprint > 0)
        {
          fprintf (fp, "we have %i low energy neighbors to this flashpoint\n", nsto);
          for (kk = 0; kk < nsto; kk++)
            {
              fprintf (fp, "[%2i] %6.2f %6.2f %6.2f ", pair_can[kk], shellhit->XX[pair_can[kk]], shellhit->YY[pair_can[kk]],
                       shellhit->ZZ[pair_can[kk]]);
              fprintf (fp, "[%6.3f]\n", shellhit->edet[pair_can[kk]]);
            };
        };

      /* check neighbors */

      ninterest = 0;
      for (jj = 0; jj < nsto; jj++)
        used[jj] = 0;

      /* search for 0.511s in candicates */

      eexpect = 0.511;
#include "findpairprod2.h"
      if (ninterest == 1)
        ok = 2;
      else if (ninterest == 2)
        ok = 3;

//      for (i = 0; i < nsto; i++)
//        fprintf (fp, "sto = %i, used=%i\n", i, used[i]);


      if (ok == 1)
        {
          /* search for 1.022 in candicates */
          /* if we found no 511s */

          eexpect = 1.022;
#include "findpairprod2.h"
          if (ninterest == 1)
            ok = 4;
        };

      if (nprint > 0)
        fprintf (fp, "ok=%i, ninterest=%i\n", ok, ninterest);

      if (ok > 0)
        ctkStat->ppstat[3 + ok]++;


#if(1)
      /* now require one 511, two 511 or one 1022 */
      /* to actually process as a pp event */

      if (ok == 1)
        {

          /* clear the shell used field */
          /* so event can be used as Compton */

          npp = 0;
          for (ll = 0; ll < shellhit->nhit; ll++)
            shellhit->used[ll] = 0;

          /* quit out */

          if (nprint > 0)
            {
              fprintf (fp, "quit, as OK==1, npp=%i\n", npp);
              printShellHit (fp, shellhit);
            };

          return (0);

        }
#endif

      else if (ok >= 2)
        {

          /* we have a solution */
          /* now not just a flash point */

          /* find the gamma energy */

          efinal = flashe[ll];
          for (ii = 0; ii < nsto; ii++)
            if (used[ii] == 1)
              {
                efinal += shellhit->edet[pair_can[ii]];
                 shellhit->used[pair_can[ii]] = 1;
              };
#if(1)

          /* add in the rest (unused) of the sto points */
          /* this may or may not be the right thing to do */

          for (ii = 0; ii < nsto; ii++)
            if (used[ii] == 0)
              {
                efinal += shellhit->edet[pair_can[ii]];
                shellhit->used[pair_can[ii]] = 1;
              };

#endif

          /* store for writeout */

          ppe[npp] = efinal;
          ppx[npp] = flashx[ll];
          ppy[npp] = flashy[ll];
          ppz[npp] = flashz[ll];
          ppTS[npp] = flashTS[ll];

          /* for now, just encode FOM value */

          if (ok == 1)
            ppfom[npp] = 0.70;
          else if (ok == 2)
            ppfom[npp] = 0.65;
          else if (ok == 3)
            ppfom[npp] = 0.60;
          else if (ok == 4)
            ppfom[npp] = 0.55;
          else
            ppfom[npp] = 0.50;
          npp++;

          /* debug print */

          if (nprint > 0)
            {
              fprintf (fp, "Pair product gamma energy: %6.3f\n", efinal);
              fprintf (fp, "shell hit now >>\n");
              for (ii = 0; ii < shellhit->nhit; ii++)
                {
                  fprintf (fp, "sh - [%2i] %6.2f %6.2f %6.2f ", ii, shellhit->XX[ii], shellhit->YY[ii], shellhit->ZZ[ii]);
                  fprintf (fp, " [%6.3f] {%i} TS=%lli\n", shellhit->edet[ii], shellhit->used[ii], shellhit->timestamp[ii]);
                };
              fprintf (fp, "\n");

              fprintf (fp, "passed on pair tracked gamma ray:\n");
              for (ii = 0; ii < npp; ii++)
                {
                  fprintf (fp, "npp=%i ", ii);
                  fprintf (fp, " %6.2f %6.2f %6.2f ", ppx[ii], ppy[ii], ppz[ii]);
                  fprintf (fp, " e= %6.3f, TS=%lli\n", ppe[ii], ppTS[ii]);
                };

              fprintf (fp, "done\n\n");
            };


        }
      else
        {
          for (ll = 0; ll < shellhit->nhit; ll++)
            shellhit->used[ll] = 0;
          if (nprint > 0)
            fprintf (fp, "reached 'else' point\n");
          return (0);

        }

    };

  /* done */

  if (nprint > 0)
    nprint--;

  return (0);

};
