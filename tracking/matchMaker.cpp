
#include <stdlib.h>
#include <fcntl.h>
#include <stddef.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>

#ifdef SOLARIS
#include <strings.h>
#include <string.h>
#endif

#ifdef LINUX
#include <string.h>
#endif

#include "ctk.h"

/*---------*/
/* globals */
/*---------*/

#define DEBUG 0

volatile extern TRACKINGPAR Pars;       /* tracking parameters */

/*--------------------------------------------------------*/

int
matchMaker (int ii, int jj, int evno, STAT * ctkStat, CLUSTER_INTPTS Clstr[MAXCLUSTERHITS],
            int *nClusters, float target_pos[3])
{
  /* declarations */

  float d1, dist;
  int tryClstr, j;
  static int nn=0;

  nn++;

  if (Pars.nprint > 0 || DEBUG)
    {
      fprintf (stderr,"** matchMaker called at event no %i, ", evno);
      fprintf (stderr,"combine cluster %i and %i?\n", ii, jj);
    };

  /* find distance beween them */
  /* return if they are too far from one another */
  /* Pars.matchmakerMaxDist is already squared */
  /* dist = sqrtf (dist); */
  /* kickout as soon as we can */

  dist = 0;

  d1 = Clstr[ii].intpts[0].xx - Clstr[jj].intpts[0].xx;
  dist += d1 * d1;
  if (Pars.nprint > 0 || DEBUG)
    printf("1,dist=%f, Pars.matchmakerMaxDist=%f\n", dist,Pars.matchmakerMaxDist);
  if (dist > Pars.matchmakerMaxDist)
    return (1);

  d1 = Clstr[ii].intpts[0].yy - Clstr[jj].intpts[0].yy;
  dist += d1 * d1;
  if (Pars.nprint > 0 || DEBUG)
    printf("2,dist=%f, Pars.matchmakerMaxDist=%f\n", dist,Pars.matchmakerMaxDist);
  if (dist > Pars.matchmakerMaxDist)
    return (1);

  d1 = Clstr[ii].intpts[0].zz - Clstr[jj].intpts[0].zz;
  dist += d1 * d1;
  if (Pars.nprint > 0 || DEBUG)
    printf("3,dist=%f, Pars.matchmakerMaxDist=%f\n", dist,Pars.matchmakerMaxDist);
  if (dist > Pars.matchmakerMaxDist)
    return (1);

  /* if we get here, the two single hits */
  /* are close enought to warrent making */
  /* a trial cluster and track it */

  tryClstr = (*nClusters);

  Clstr[tryClstr].ndet = 2;
  Clstr[tryClstr].esum = Clstr[ii].esum + Clstr[jj].esum;
//printf("esum: %f %f --> %f\n", Clstr[ii].esum, Clstr[jj].esum, Clstr[tryClstr].esum);

  Clstr[tryClstr].intpts[0].shellHitPos = Clstr[ii].intpts[0].shellHitPos;
  Clstr[tryClstr].intpts[0].xx = Clstr[ii].intpts[0].xx;
  Clstr[tryClstr].intpts[0].yy = Clstr[ii].intpts[0].yy;
  Clstr[tryClstr].intpts[0].zz = Clstr[ii].intpts[0].zz;
  Clstr[tryClstr].intpts[0].edet = Clstr[ii].intpts[0].edet;
  Clstr[tryClstr].intpts[0].timestamp = Clstr[ii].intpts[0].timestamp;
  Clstr[tryClstr].intpts[0].detno = Clstr[ii].intpts[0].detno;

  Clstr[tryClstr].intpts[1].shellHitPos = Clstr[ii].intpts[0].shellHitPos;
  Clstr[tryClstr].intpts[1].xx = Clstr[jj].intpts[0].xx;
  Clstr[tryClstr].intpts[1].yy = Clstr[jj].intpts[0].yy;
  Clstr[tryClstr].intpts[1].zz = Clstr[jj].intpts[0].zz;
  Clstr[tryClstr].intpts[1].edet = Clstr[jj].intpts[0].edet;
  Clstr[tryClstr].intpts[1].timestamp = Clstr[jj].intpts[0].timestamp;
  Clstr[tryClstr].intpts[1].detno = Clstr[jj].intpts[0].detno;

  /* initialize trial cluster for tracking */

  for (j = 0; j < Clstr[tryClstr].ndet; j++)
    Clstr[tryClstr].intpts[j].order = -1;
  Clstr[tryClstr].valid = 1;
  Clstr[tryClstr].fom = MAXFOM;
  Clstr[tryClstr].tracked = 0;

#if (DEBUG)
  fprintf (stderr,"combined trial cluster: Clstr[tryClstr].ndet=%i\n", Clstr[tryClstr].ndet);
  printCluster (tryClstr, stdout, Clstr);
#endif

  /* track the cluster */

#if (DEBUG)
  fprintf (stderr,"tracking...Pars.trackOps[Clstr[tryClstr].ndet=%i\n",Pars.trackOps[Clstr[tryClstr].ndet]);
  fflush (stdout);
#endif
  ctksort (tryClstr, Clstr, nClusters);
  switch (Pars.trackOps[Clstr[tryClstr].ndet])
    {
    case 0:
      ctktk0 (tryClstr, target_pos, ctkStat, Clstr, nClusters);
      break;
    case 1:
      ctktk1 (tryClstr, target_pos, ctkStat, Clstr, nClusters);
      break;
    case 3:
      ctktk3 (tryClstr, target_pos, ctkStat, Clstr, nClusters);
      break;
    case 4:
      ctktk4 (tryClstr, target_pos, ctkStat, Clstr, nClusters);
      break;
    case 5:
      ctktk5 (tryClstr, target_pos, ctkStat, Clstr, nClusters);
      break;
    default:
      fprintf (stderr,"ctk: tracking option not known!?\n, Quit");
      exit (1);
    };

  /* if OK, promote to valid cluster */

  if (Clstr[tryClstr].fom < Pars.matchmaker_kickoutfom)
    {
      /* validate the new combined cluster */

      Clstr[tryClstr].valid = 1;
      (*nClusters)++;

      /* invalidate the two clusters we combined from */

      Clstr[ii].valid = 0;
      Clstr[jj].valid = 0;

      if (Pars.nprint > 0 || DEBUG)
        {
          fprintf (stderr,"...matchMaker: success! made the combination> \n");
          fflush (stdout);
          printCluster (tryClstr, stdout, Clstr);
        };

      return (0);
    }
  else
    {
      Clstr[tryClstr].valid = 0;

      if (Pars.nprint > 0 || DEBUG)
        {
          fprintf (stderr,"...matchMaker: combination not OK, FOM=%f \n", Clstr[tryClstr].fom);
          fflush (stdout);
        };

      return (2);

    };



  if (1)
    exit (0);

  /* done */

  return (0);

};
