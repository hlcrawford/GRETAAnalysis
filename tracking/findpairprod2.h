  /* is there a more efficient way to do this??? */

  /* see if any single interaction has eexpect */


if (nsto >= 1)
    for (i = 0; i < nsto; i++)
      if (used[i] == 0)
        {
          etest = shellhit->edet[pair_can[i]];
          if (fabs (etest - eexpect) < Pars.ppde)
            {
              ninterest++;
              used[i] =  1;
              if (nprint > 0)
                fprintf (fp, "etest= %6.3f is interesting %i\n", etest, pair_can[i]);
            };
        };

  /* see if any pairs sum up to eexpect */

if (nsto >= 2)
    for (i = 0; i < nsto; i++)
      for (j = i + 1; j < nsto; j++)
        if (used[i] == 0)
          if (used[j] == 0)
            {
              etest = shellhit->edet[pair_can[i]] + shellhit->edet[pair_can[j]];
              if (fabs (etest - eexpect) < Pars.ppde)
                {
                  ninterest++;
                  used[i] =  1;
                  used[j] =  1;
                  if (nprint > 0)
                    fprintf (fp, "etest= %6.3f is interesting %i %i\n", etest, pair_can[i], pair_can[j]);
                };
            };

  /* see if any triples sum up to eexpect */


if (nsto >= 3)
    for (i = 0; i < nsto; i++)
      for (j = i + 1; j < nsto; j++)
        for (k = j + 1; k < nsto; k++)
          if (used[i] == 0)
            if (used[j] == 0)
              if (used[k] == 0)
                {
                  etest = shellhit->edet[pair_can[i]] + shellhit->edet[pair_can[j]] + shellhit->edet[pair_can[k]];
                  if (fabs (etest - eexpect) < Pars.ppde)
                    {
                      ninterest++;
                      used[i] =  1;
                      used[j] =  1;
                      used[k] =  1;
                      if (nprint > 0)
                        fprintf (fp, "etest= %6.3f is interesting %i %i %i\n", etest, pair_can[i], pair_can[j], pair_can[k]);
                    };
                };



  /* see if any quadruples sum up to eexpect */


if (nsto >= 4)
    for (i = 0; i < nsto; i++)
      for (j = i + 1; j < nsto; j++)
        for (k = j + 1; k < nsto; k++)
          for (l = k + 1; l < nsto; l++)
            if (used[i] == 0)
              if (used[j] == 0)
                if (used[k] == 0)
                  if (used[l] == 0)
                    {
                      etest = shellhit->edet[pair_can[i]] + shellhit->edet[pair_can[j]] + shellhit->edet[pair_can[k]] + shellhit->edet[pair_can[l]];
                      if (fabs (etest - eexpect) < Pars.ppde)
                        {
                          ninterest++;
                          used[i] =  1;
                          used[j] =  1;
                          used[k] =  1;
                          used[l] =  1;
                          if (nprint > 0)
                            fprintf (fp, "etest= %6.3f is interesting %i %i %i %i\n", etest, pair_can[i], pair_can[j], pair_can[k], pair_can[l]);
                        };
                    };


  /* see if any femtuples sum up to eexpect */


if (nsto >= 5)
    for (i = 0; i < nsto; i++)
      for (j = i + 1; j < nsto; j++)
        for (k = j + 1; k < nsto; k++)
          for (l = k + 1; l < nsto; l++)
            for (m = l + 1; m < nsto; m++)
              if (used[i] == 0)
                if (used[j] == 0)
                  if (used[k] == 0)
                    if (used[l] == 0)
                      if (used[m] == 0)
                        {
                          etest =
                            shellhit->edet[pair_can[i]] + shellhit->edet[pair_can[j]] + shellhit->edet[pair_can[k]] + shellhit->edet[pair_can[l]] +
                            shellhit->edet[pair_can[m]];
                          if (fabs (etest - eexpect) < Pars.ppde)
                            {
                              ninterest++;
                              used[i] =  1;
                              used[j] =  1;
                              used[k] =  1;
                              used[l] =  1;
                              used[m] =  1;
                              if (nprint > 0)
                                fprintf (fp, "etest= %6.3f is interesting %i %i %i %i %i\n", etest, pair_can[i], pair_can[j], pair_can[k],
                                         pair_can[l], pair_can[m]);
                            };
                        };


  /* see if any sixtuples sum up to eexpect */


if (nsto >= 6)
    for (i = 0; i < nsto; i++)
      for (j = i + 1; j < nsto; j++)
        for (k = j + 1; k < nsto; k++)
          for (l = k + 1; l < nsto; l++)
            for (m = l + 1; m < nsto; m++)
              for (n = m + 1; n < nsto; n++)
                if (used[i] == 0)
                  if (used[j] == 0)
                    if (used[k] == 0)
                      if (used[l] == 0)
                        if (used[m] == 0)
                          if (used[n] == 0)
                            {
                              etest =
                                shellhit->edet[pair_can[i]] + shellhit->edet[pair_can[j]] + shellhit->edet[pair_can[k]] +
                                shellhit->edet[pair_can[l]] + shellhit->edet[pair_can[m]] + shellhit->edet[pair_can[n]];
                              if (fabs (etest - eexpect) < Pars.ppde)
                                {
                                  ninterest++;
                                  used[i] =  1;
                                  used[j] =  1;
                                  used[k] =  1;
                                  used[l] =  1;
                                  used[m] =  1;
                                  used[n] =  1;
                                  if (nprint > 0)
                                    fprintf (fp, "etest= %6.3f is interesting %i %i %i %i %i %i\n", etest, pair_can[i], pair_can[j], pair_can[k],
                                             pair_can[l], pair_can[m], pair_can[n]);
                                };
                            };

