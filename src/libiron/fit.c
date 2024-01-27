// fit
fit() {
eqn_windows = (struct eqn_window *)h_subcalloc(mempool,
               num_filtered,
               sizeof(struct eqn_window));
               if (eqn_windows == NULL)
               {
                 fprintf(stderr, "ERROR -- calloc failed\n");

                 free(signal_pairs);
                 free(temp_ptrs1);
                 free(temp_ptrs2);
                 free(temp_ptrs3);
                 free(temp_ptrs4);

                 return (1);
               }

               /* initialize values for smoothed piecewise linear fit (geometric) */
               /* calculate for ALL points, since we will eventually use them later */
               for (i = 0; i < num_spots; i++)
               {
                 pair_ptr = &signal_pairs[i];
                 pair_ptr->log_adjust = log(pair_ptr->sig1 / pair_ptr->sig2);
               }

               /* windowed linear fits, log(x/y) vs. log(x*y) */
               /* use filt2, since filt1 is already sorted on X which we'll use later */

               num_eqns = fill_geometric_eqn_windows(eqn_windows, filt2, num_filtered,
                                                     flags->iron_fit_window_frac, 1,
                                                     weight_exponent);
               smooth_geometric_fits(eqn_windows, num_eqns, filt2, num_filtered);

#if 0
               /* calculate post-normalized rmsd from center line */
               rmsd = 0.0;
               for (i = 0; i < num_filtered; i++)
                 rmsd += filt1[i]->norm_err * filt1[i]->norm_err;

               if (num_filtered)
                 rmsd = sqrt(rmsd / num_filtered);

               /* filter training points beyond 1 rmsd, store in filt2 */
               /* filt2 will contain final training set, filt1 contains previous set */
               old_num_filtered = num_filtered;

               num_filtered = 0;
               for (i = 0; i < old_num_filtered; i++)
               {
                 if (fabs(filt1[i]->norm_err) < 5*rmsd + 1E-5)
                   filt2[num_filtered++] = filt1[i];
               }

               /* refit on final reduced training set */
               num_eqns = fill_geometric_eqn_windows(eqn_windows, filt2,
                                                     num_filtered, flags->iron_fit_window_frac,
                                                     2, weight_exponent);
               smooth_geometric_fits(eqn_windows, num_eqns, filt2, num_filtered);
#endif

               /* remember which points were in the irank training set */
               for (i = 0; i < num_filtered; i++)
                 filt2[i]->irank_flag = 1;

               /* calculate adjustments for ALL points */
               /* store pointers to ALL points in filt1 */
               for (i = 0; i < num_spots; i++)
                 filt1[i] = signal_pairs + i;

               interpolate_final_scales(filt1, num_spots, filt2, num_filtered,
                                        fit_both_x_y_flag);

               /* use a single global scaling factor, rather than non-linear scaling */
               if (global_scaling_flag)
               {
#if 0
                 /* refit linear line on final reduced training set, single full window */
                 /* appears to give a worse fit than using non-linear fit (?) */
                 num_eqns = fill_geometric_eqn_windows(eqn_windows, filt2,
                                                       num_filtered, 1.0, 1,
                                                       weight_exponent);
                 smooth_geometric_fits(eqn_windows, num_eqns, filt2, num_filtered);
                 interpolate_final_scales(filt1, num_spots, filt2, num_filtered,
                                          fit_both_x_y_flag, err);
#endif

                 int32_t count = 0;

                 for (i = 0; i < num_spots; i++)
                 {
                   if (signal_pairs[i].irank_flag)
                   {
                     global_scale += signal_pairs[i].fit_log_adjust;
                     count++;
                   }
                 }

                 global_scale = exp(global_scale / count);

                 fprintf(stderr, "GlobalScale:\t%s\t%f\t%f\t%d\t%d\t%d\t%d\t%f\n",
                         filestem,
                         global_scale, log(global_scale) / log(2.0),
                         count, num_both_not_weak, num_not_weak, num_spots,
                         (double) count / (double) num_both_not_weak);
               }
               /* use a single line fit to the entire training set */
               else if (flags->iron_untilt_normalization)
               {
                 int32_t count = 0;

                 /* refit linear line on final reduced training set, single full window */
                 num_eqns = fill_geometric_eqn_windows(eqn_windows, filt2,
                                                       num_filtered, 1.0, 1,
                                                       weight_exponent);
                 smooth_geometric_fits(eqn_windows, num_eqns, filt2, num_filtered);
                 interpolate_final_scales(filt1, num_spots, filt2, num_filtered,
                                          fit_both_x_y_flag, err);

                 /* average adjustements together for printing QC info */
                 global_scale = 0.0;
                 for (i = 0; i < num_spots; i++)
                 {
                   if (signal_pairs[i].irank_flag)
                   {
                     global_scale += signal_pairs[i].fit_log_adjust;
                     count++;
                   }
                 }

                 global_scale = exp(global_scale / count);

                 fprintf(stderr, "GlobalFitLine:\t%s\t%f\t%f\t%f\t%d\t%d\t%d\n",
                         filestem,
                         1.0 / global_scale,
                         -log(global_scale) / log(2.0),
                         -(180.0 * atan(eqn_windows[0].slope) / M_PI),
                         num_both_not_weak, num_not_weak, num_spots);
               }

               /* store scaling multipliers in signals2_scales array */
               for (i = 0; i < num_spots; i++)
               {
#if DO_FLOOR
                 /* zero out extremely weak signals */
                 if (signal_pairs[i].sig2 <= MIN_SIGNAL)
                   signals2_scales[i] = 0;
                 else
#endif

                   if (global_scaling_flag)
                     signals2_scales[i] = global_scale;
                   else
                     signals2_scales[i] = exp(signal_pairs[i].fit_log_adjust);
               }


#if DEBUG_FILE
               debug_file = fopen("irank_set.txt", "wb");
               if (!debug_file)
                 fprintf(stderr, "ERROR -- Can't open irank_set.txt\n");

               fprintf(debug_file, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                       "ProbeID", "log10_X", "log10_Y", "log10_Y_norm",
                       "Weight",
                       "IRankSet", "InitialRankSet", "IRankFrac", "IRankFrac_0",
                       "X0_proj", "Y0_proj", "Y0_proj_norm",
                       "log10(X*Y)", "log10(X/Y)", "log10(X/Y_norm)",
                       "log10(X0_proj*Y0_proj)", "log10(X0_proj/Y0_proj)",
                       "log10(X0_proj_norm/Y0_proj_norm)");

               /* find lowest non-zero rank_diff */
               temp = 1.0;
               for (i = 0; i < num_spots; i++)
               {
                 if (signal_pairs[i].irank_frac > 0 &&
                     signal_pairs[i].irank_frac < temp)
                 {
                   temp = signal_pairs[i].irank_frac;
                 }
                 if (signal_pairs[i].irank_frac_0 > 0 &&
                     signal_pairs[i].irank_frac_0 < temp)
                 {
                   temp = signal_pairs[i].irank_frac_0;
                 }
               }

               for (i = 0; i < num_spots; i++)
               {
                 double temp_signal, irank_frac, irank_frac_0;
                 double x_train, y_train, y_train_norm;

                 if (mask_array[i])
                   continue;

                 temp_signal = signal_pairs[i].sig2 * signals2_scales[i];
                 if (temp_signal < MIN_SIGNAL)
                   temp_signal = MIN_SIGNAL;

                 x_train = log(signal_pairs[i].sig1) / log(10.0);
                 y_train = log(signal_pairs[i].sig2) / log(10.0);
                 y_train_norm = log(temp_signal) / log(10.0);

                 if (signal_pairs[i].irank_flag)
                 {
                   x_train = (log(signal_pairs[i].sig1) + 0.5 * signal_pairs[i].norm_err) /
                     log(10.0);
                   y_train = (log(signal_pairs[i].sig2) - 0.5 * signal_pairs[i].norm_err) /
                     log(10.0);
                   y_train_norm = x_train;
                 }

                 irank_frac = signal_pairs[i].irank_frac;
                 irank_frac_0 = signal_pairs[i].irank_frac_0;

                 if (irank_frac < temp) irank_frac = temp;
                 if (irank_frac_0 < temp) irank_frac_0 = temp;

                 fprintf(debug_file, "%d\t%f\t%f\t%f\t%e\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
                         i,
                         log(signal_pairs[i].sig1) / log(10.0),
                         log(signal_pairs[i].sig2) / log(10.0),
                         log(temp_signal) / log(10.0),
#if 0
                         fabs(signal_pairs[i].norm_err_scaled),
#else
                         signal_pairs[i].weight,
#endif

#if DEBUG_COLOR_IRANK
                         signal_pairs[i].irank_flag,
                         signal_pairs[i].initial_set_flag,
                         log10(irank_frac),
                         log10(irank_frac_0),
#else
                         0.0, 0.0, 0.0, 0.0,
#endif
                         x_train,
                         y_train,
                         y_train_norm,
                         log(signal_pairs[i].sig1 * signal_pairs[i].sig2) / log(10.0),
                         log(signal_pairs[i].sig1 / signal_pairs[i].sig2) / log(10.0),
                         log(signal_pairs[i].sig1 / temp_signal) / log(10.0),
                         x_train + y_train,
                         x_train - y_train,
                         x_train - y_train_norm);
                 /*
                  * some editors' auto-indentation freaks out unless you remove the
                  * ); from the #ifdef
                  */
               }
               fclose(debug_file);
#endif


}



static int fill_geometric_eqn_windows(struct eqn_window *eqn_windows,
                                      struct signal_pair **filt_ptrs,
                                      int num_pairs, double window_frac,
                                      int pass, double weight_exponent)
{
  struct signal_pair *pair_ptr;
  int                 n, i;
  int                 w = (int)(window_frac * num_pairs + 0.5);
  int                 wsmall = (int)(0.01 * num_pairs + 0.5);
  double              x, y;
  double              ss_xx = 0.0, ss_xy = 0.0;
  double              x_sum = 0.0, y_sum = 0.0;
  double              x_avg, y_avg;
  double              weight, weight_sum = 0.0;
  double              min_weight = 9E99, max_weight = -9E99;

  assert(filt_ptrs   != NULL);
  assert(eqn_windows != NULL);

  if (w < 100)
    w = 100;
  if (wsmall < 10)
    wsmall = 10;

  if (w > num_pairs)
    w = num_pairs;
  if (wsmall > num_pairs)
    wsmall = num_pairs;

  qsort(filt_ptrs,
        num_pairs,
        sizeof(struct signal_pair *),
        compare_sort_log_xy);


  /* weights */

  /* initialize weights */
  for (i = 0; i < num_pairs; i++)
  {
    filt_ptrs[i]->weight = 0;
    filt_ptrs[i]->n_windows = 0;
  }

  /* slide bin windows */
  // There is a running window of size wsmall, starting from (sorted
  // log xy) the beginning throughout.
  for (n = 0; n <= num_pairs - wsmall; n++)
  {
    x_avg = 0;
    for (i = n; i < n + wsmall; i++)
      x_avg += filt_ptrs[i]->log_xy;
    x_avg /= wsmall;

    weight = 0;
    for (i = n; i < n + wsmall; i++)
    {
      x = filt_ptrs[i]->log_xy - x_avg;
      weight += x*x;
    }
    weight = sqrt(weight / wsmall);

    // x_avg is actually logxy average
    // weight is stddev
    // Both are (unless there is trickery elsewhere) logxy
    for (i = n; i < n + wsmall; i++)
    {
      filt_ptrs[i]->weight += weight;
      filt_ptrs[i]->n_windows++;
    }
    // Since each point is associated with multiple, overlapping
    // windows, apparently the point is cumulative stddev and
    // the total number of windows is tracked.
  }

  // OK, after going through all points, calculate the average stddev
  // for each point. But I'm not sure it's like each point is sort
  // of a model although I'm not positive
  for (i = 0; i < num_pairs; i++)
  {
    filt_ptrs[i]->weight /= filt_ptrs[i]->n_windows;

    if (filt_ptrs[i]->weight >= 1E-5 && filt_ptrs[i]->weight < min_weight)
      min_weight = filt_ptrs[i]->weight;

    if (filt_ptrs[i]->weight > max_weight)
      max_weight = filt_ptrs[i]->weight;
  }

  if (min_weight == 9E99)
    min_weight = max_weight;

#if DEBUG_PRINT
  fprintf(stderr, "Weights:\t%f\t%f\t%f\n",
          min_weight, max_weight, max_weight / min_weight);
#endif

  // This is somehow adjusting the weight (average stddev) by exponent (maybe 4 by default),
  // But it's somehow the proportional stddev (divided by max weight).
  for (i = 0; i < num_pairs; i++)
  {
    if (filt_ptrs[i]->weight < 1E-5)
      filt_ptrs[i]->weight = min_weight;

    /*
     * w^4 = sigma^4 = variance^2; w^8 = (variance^2)^2
     */
    filt_ptrs[i]->weight = pow(filt_ptrs[i]->weight / max_weight, weight_exponent);
  }


  // Near as I can tell, all of the above is computing a value for the weight.
  // Somehow a function of stddev, but not sure how.


  // w is at least 100, if not bigger.
  // I guess under the assumption that w is at least the number of
  // training points.
  // Yes, w will get adjusted down to number of pairs if training
  // points are too small.

  // THe below is a little weird. Cumulative weight, cumulative weighted sum,
  // cumulative weighted sum squared. It must go into the next block, where
  // things are actually calculated. What I don't understand yet is why
  // this is a cumulative across everything (not windowed) except w isn't
  // guaranteed to be all the points?? So like if w is less than num_pairs
  // which I think can happen, only the first w points are considered in
  // this computation?


  /* local windowed fits */
  n = 0;
  weight_sum = 0;

  for (i = 0; i < w; i++)
  {
    pair_ptr = filt_ptrs[i];
    weight = pair_ptr->weight;
    weight_sum += weight;

    x      = pair_ptr->log_xy;
    x_sum += weight * x;
    ss_xx += weight * x*x;

    y      = pair_ptr->log_adjust;
    y_sum += weight * y;
    ss_xy += weight * x*y;
  }


  // Not really sure about this bit, it's an infinte loop with a break
  // which should probably be changed. log(x*y) vs log(x/y)
  while (1)
  {
    double temp;

    x_avg = x_sum / weight_sum;
    y_avg = y_sum / weight_sum;

    temp = ss_xx - weight_sum * x_avg * x_avg;
    if (temp)
      eqn_windows[n].slope = (ss_xy - weight_sum * x_avg * y_avg) / temp;
    else
      eqn_windows[n].slope = 0.0;
    eqn_windows[n].offset = y_avg - eqn_windows[n].slope * x_avg;
    eqn_windows[n].start  = filt_ptrs[n]->log_xy;
    eqn_windows[n].end    = filt_ptrs[n+w-1]->log_xy;

    if (n >= num_pairs - w)
    {
      n++;
      break;
    }

    /* add next point to sums */
    pair_ptr = filt_ptrs[n+w];
    weight = pair_ptr->weight;
    weight_sum += weight;

    x      = pair_ptr->log_xy;
    x_sum += weight * x;
    ss_xx += weight * x*x;

    y      = pair_ptr->log_adjust;
    y_sum += weight * y;
    ss_xy += weight * x*y;


    /* remove first point from sums */
    pair_ptr = filt_ptrs[n];
    weight = pair_ptr->weight;
    weight_sum -= weight;

    x      = pair_ptr->log_xy;
    x_sum -= weight * x;
    ss_xx -= weight * x*x;

    y      = pair_ptr->log_adjust;
    y_sum -= weight * y;
    ss_xy -= weight * x*y;

    n++;
  }

  return (n);
}


static void smooth_geometric_fits(struct eqn_window *eqn_windows,
                                  int num_eqn_windows,
                                  struct signal_pair **filt_ptrs,
                                  int num_pairs)
{
  struct signal_pair *pair_ptr;
  double              x, avg_adjust, sum_slope = 0, sum_offset = 0;
  int                 i, j, min_eqn_idx = 0, end_eqn_idx = 0;
  int                 old_min_eqn_idx = 0, old_end_eqn_idx = 0;

  assert(filt_ptrs   != NULL);
  assert(eqn_windows != NULL);

  for (i = 0; i < num_pairs; i++)
  {
    pair_ptr = filt_ptrs[i];
    x        = pair_ptr->log_xy;

#if 0
    if (isinf(x))
      printf("NAN %d %f %f\n", i, pair_ptr->sig1, pair_ptr->sig2);
#endif

    old_min_eqn_idx = min_eqn_idx;
    old_end_eqn_idx = end_eqn_idx;

    /* find first possibly overlapping eqn window */
    while (min_eqn_idx < num_eqn_windows && eqn_windows[min_eqn_idx].end < x)
      min_eqn_idx++;

    if (end_eqn_idx < min_eqn_idx)
      end_eqn_idx = min_eqn_idx;

    while (end_eqn_idx < num_eqn_windows &&
           x >= eqn_windows[end_eqn_idx].start &&
           x <= eqn_windows[end_eqn_idx].end)
      end_eqn_idx++;

    /* sum slope and offset */
    if (i == 0)
    {
      for (j = min_eqn_idx; j < end_eqn_idx; j++)
      {
        sum_slope += eqn_windows[j].slope;
        sum_offset += eqn_windows[j].offset;
      }
    }
    /* subtract out old window, add in new window */
    else
    {
      for (j = old_min_eqn_idx; j < min_eqn_idx; j++)
      {
        sum_slope -= eqn_windows[j].slope;
        sum_offset -= eqn_windows[j].offset;
      }

      for (j = old_end_eqn_idx; j < end_eqn_idx; j++)
      {
        sum_slope += eqn_windows[j].slope;
        sum_offset += eqn_windows[j].offset;
      }
    }

    /* average the fit adjust value over all overlapping eqns */
    avg_adjust = (sum_slope * x + sum_offset) / (end_eqn_idx - min_eqn_idx);

    pair_ptr->fit_log_adjust = avg_adjust;
    pair_ptr->norm_err       = avg_adjust - pair_ptr->log_adjust;

#if 0
    if (isnan(pair_ptr->fit_log_adjust))
      printf("NAN_LOGADJUST %f\n", x);
    if (isnan(avg_adjust))
      printf("NAN_AVGADJUST %f\n", x);
#endif
  }
}


/* return 1 on error, 0 if no errors*/
static int interpolate_final_scales(struct signal_pair **pair_ptrs,
                                    int num_pairs,
                                    struct signal_pair **filt_ptrs_train,
                                    int num_pairs_train,
                                    int fit_both_x_y_flag)
{
  struct xy_pair    *xy_pairs = NULL;
  struct xy_pair    *xy_pairs2 = NULL;
  struct signal_pair *pair_ptr;
  double              x, y, old_x, old_y, sum;
  int                 num_pairs_train2 = 0;
  int                 num_pairs_train3 = 0;
  int                 i, num;
  int                 min_idx, last_idx;

  assert(pair_ptrs       != NULL);
  assert(filt_ptrs_train != NULL);


  xy_pairs = (struct xy_pair *) calloc(num_pairs_train,
              sizeof(struct xy_pair));
  if (xy_pairs == NULL)
  {
    fprintf(stderr, "ERROR -- calloc failed\n");
    return (1);
  }

  qsort(pair_ptrs,
        num_pairs,
        sizeof(struct signal_pair *),
        compare_sort_sig2);
  qsort(filt_ptrs_train,
        num_pairs_train,
        sizeof(struct signal_pair *),
        compare_sort_sig2);

  /* condense overlapping points */
  /* project onto average best fit line */
  old_x = -9E99;
  old_y = -9E99;
  for (i = 0, num_pairs_train2 = 0; i < num_pairs_train; i++)
  {
    pair_ptr = filt_ptrs_train[i];

    /* project x/y onto best fit avg line using norm_err */
    x = log(pair_ptr->sig1) + 0.5 * pair_ptr->norm_err;
    y = log(pair_ptr->sig2) - 0.5 * pair_ptr->norm_err;

    /* is new point different from last point? */
    /* take round off error into account */
    if (fabs(x - old_x) > 1E-14 || fabs(y - old_y) > 1E-14)
    {
      xy_pairs[num_pairs_train2].x   = x;
      xy_pairs[num_pairs_train2++].y = y;
    }

    old_x = x;
    old_y = y;
  }

  qsort(xy_pairs,
        num_pairs_train2,
        sizeof(struct xy_pair),
        compare_sort_xy_pair_by_y);

  xy_pairs2 = (struct xy_pair *) calloc(num_pairs_train,
               sizeof(struct xy_pair));
  if (xy_pairs2 == NULL)
  {
    fprintf(stderr, "ERROR -- calloc failed\n");
    return (1);
  }


  /* store as [log(y), fit_log_adjust] pairs */
  old_x = -9E99;
  sum = 0;
  num = 0;
  for (i = 0, num_pairs_train3 = -1; i < num_pairs_train2; i++)
  {
    x = xy_pairs[i].y;

#if 0
    if (old_x == x)
    {
      fprintf(stderr, "SAMEXY     %f %f %f     %f %f %f\n",
              xy_pairs[i-1].x,
              xy_pairs[i-1].y,
              xy_pairs[i-1].x - xy_pairs[i-1].y,
              xy_pairs[i].x,
              xy_pairs[i].y,
              xy_pairs[i].x - xy_pairs[i].y);
      printf("SAMEXY     %f %f %f     %f %f %f\n",
             xy_pairs[i-1].x,
             xy_pairs[i-1].y,
             xy_pairs[i-1].x - xy_pairs[i-1].y,
             xy_pairs[i].x,
             xy_pairs[i].y,
             xy_pairs[i].x - xy_pairs[i].y);
    }
#endif

    if (x != old_x)
    {
      num_pairs_train3++;
      sum = 0;
      num = 0;
    }
    old_x = x;

    sum += xy_pairs[i].x - xy_pairs[i].y;
    num++;

    xy_pairs2[num_pairs_train3].x = x;

    /* average scales for identical points */
    if (num > 1)
      xy_pairs2[num_pairs_train3].y = sum / num;
    else
      xy_pairs2[num_pairs_train3].y = sum;
  }

  num_pairs_train3++;
  last_idx = num_pairs_train3 - 1;
  min_idx = 0;
  old_x = -9E99;


#if DEBUG_PRINT
  fprintf(stderr, "TrainingY\t%d\t%d\t%d\n",
          num_pairs_train, num_pairs_train2, num_pairs_train3);
#endif

  for (i = 0; i < num_pairs; i++)
  {
    pair_ptr = pair_ptrs[i];
    x        = log(pair_ptr->sig2);

    /* same X as last time, use the previous fit */
    if (x == old_x && i)
    {
      pair_ptr->fit_log_adjust = sum;
      continue;
    }

    old_x = x;

    /* before first training point */
    /* use average of end 10 adjustments */
    if (x < xy_pairs2[0].x)
    {
      int j;

      sum = 0.0;
      for (j = 0; j < 10 && j < last_idx; j++)
        sum += xy_pairs2[j].y;
      if (j)
        sum /= j;
      pair_ptr->fit_log_adjust = sum;

      continue;
    }

    /* after last training point */
    /* use average of end 10 adjustments */
    if (x > xy_pairs2[last_idx].x)
    {
      int j;

      sum = 0.0;
      for (j = 0; j < 10 && j < last_idx; j++)
        sum += xy_pairs2[last_idx - j].y;
      if (j)
        sum /= j;
      pair_ptr->fit_log_adjust = sum;

      continue;
    }

    /* skip to at or just after current point */
    while(xy_pairs2[min_idx].x < x)
      min_idx++;

    /* Use linear interpolation.
     * WARNING -- curve can have occasional too-rapid changes which
     * greatly confuse Lagrange interpolation and lead to big overshoots.
     * This tends to happen when the adjustment is very close to zero
     * and it oscillates a bit around zero.  Linear interpolation has
     * no issues, since it doesn't require smoothness.
     */
    sum = linear_interp(x, xy_pairs2, last_idx, min_idx);
    pair_ptr->fit_log_adjust = sum;

#if 0
    printf("DEBUG\t%f\t%d\t%f\t%f\t%f\n",
           x, min_idx, xy_pairs2[min_idx].x, xy_pairs2[min_idx].y, pair_ptr->fit_log_adjust);
#endif
  }


  /* Fit scale vs. X, fit scale vs. Y, average the scaling from both fits.
   *
   * Not usually recommended, since it will result in altered intensity rank
   * orders.
   *
   * However, for especially ill-behaved data, it can result in overall better
   * normalizations.
   */
  if (fit_both_x_y_flag)
  {
    qsort(pair_ptrs,
          num_pairs,
          sizeof(struct signal_pair *),
          compare_sort_sig1);
    qsort(filt_ptrs_train,
          num_pairs_train,
          sizeof(struct signal_pair *),
          compare_sort_sig1);

    /* condense overlapping points */
    /* project onto average best fit line */
    old_x = -9E99;
    old_y = -9E99;
    for (i = 0, num_pairs_train2 = 0; i < num_pairs_train; i++)
    {
      pair_ptr = filt_ptrs_train[i];

      /* project x/y onto best fit avg line using norm_err */
      x = log(pair_ptr->sig1) + 0.5 * pair_ptr->norm_err;
      y = log(pair_ptr->sig2) - 0.5 * pair_ptr->norm_err;

      /* is new point different from last point? */
      /* take round off error into account */
      if (fabs(x - old_x) > 1E-14 || fabs(y - old_y) > 1E-14)
      {
        xy_pairs[num_pairs_train2].x   = x;
        xy_pairs[num_pairs_train2++].y = y;
      }

      old_x = x;
      old_y = y;
    }

    qsort(xy_pairs,
          num_pairs_train2,
          sizeof(struct xy_pair),
          compare_sort_xy_pair_by_x);

    /* store as [log(x), fit_log_adjust] pairs */
    old_x = -9E99;
    sum = 0;
    num = 0;
    for (i = 0, num_pairs_train3 = -1; i < num_pairs_train2; i++)
    {
      x = xy_pairs[i].x;

#if 0
      if (old_x == x)
      {
        fprintf(stderr, "SAMEXY     %f %f %f     %f %f %f\n",
                xy_pairs[i-1].x,
                xy_pairs[i-1].y,
                xy_pairs[i-1].x - xy_pairs[i-1].y,
                xy_pairs[i].x,
                xy_pairs[i].y,
                xy_pairs[i].x - xy_pairs[i].y);
        printf("SAMEXY     %f %f %f     %f %f %f\n",
               xy_pairs[i-1].x,
               xy_pairs[i-1].y,
               xy_pairs[i-1].x - xy_pairs[i-1].y,
               xy_pairs[i].x,
               xy_pairs[i].y,
               xy_pairs[i].x - xy_pairs[i].y);
      }
#endif

      if (x != old_x)
      {
        num_pairs_train3++;
        sum = 0;
        num = 0;
      }
      old_x = x;

      sum += xy_pairs[i].x - xy_pairs[i].y;
      num++;

      xy_pairs2[num_pairs_train3].x = x;

      /* average scales for identical points */
      if (num > 1)
        xy_pairs2[num_pairs_train3].y = sum / num;
      else
        xy_pairs2[num_pairs_train3].y = sum;
    }

    num_pairs_train3++;
    last_idx = num_pairs_train3 - 1;
    min_idx = 0;
    old_x = -9E99;


#if DEBUG_PRINT
    fprintf(stderr, "TrainingX\t%d\t%d\t%d\n",
            num_pairs_train, num_pairs_train2, num_pairs_train3);
#endif

    for (i = 0; i < num_pairs; i++)
    {
      pair_ptr = pair_ptrs[i];
      x        = log(pair_ptr->sig1);

      /* same X as last time, use the previous fit */
      if (x == old_x && i)
      {
        pair_ptr->fit_log_adjust = 0.5 * (pair_ptr->fit_log_adjust + sum);
        continue;
      }

      old_x = x;

      /* before first training point */
      /* use average of end 10 adjustments */
      if (x < xy_pairs2[0].x)
      {
        int j;

        sum = 0.0;
        for (j = 0; j < 10 && j < last_idx; j++)
          sum += xy_pairs2[j].y;
        if (j)
          sum /= j;
        pair_ptr->fit_log_adjust = 0.5 * (pair_ptr->fit_log_adjust + sum);

        continue;
      }

      /* after last training point */
      /* use average of end 10 adjustments */
      if (x > xy_pairs2[last_idx].x)
      {
        int j;

        sum = 0.0;
        for (j = 0; j < 10 && j < last_idx; j++)
          sum += xy_pairs2[last_idx - j].y;
        if (j)
          sum /= j;
        pair_ptr->fit_log_adjust = 0.5 * (pair_ptr->fit_log_adjust + sum);

        continue;
      }

      /* skip to at or just after current point */
      while(xy_pairs2[min_idx].x < x)
        min_idx++;

      /* Use linear interpolation.
       * WARNING -- curve can have occasional too-rapid changes which
       * greatly confuse Lagrange interpolation and lead to big overshoots.
       * This tends to happen when the adjustment is very close to zero
       * and it oscillates a bit around zero.  Linear interpolation has
       * no issues, since it doesn't require smoothness.
       */
      sum = linear_interp(x, xy_pairs2, last_idx, min_idx);
      pair_ptr->fit_log_adjust = 0.5 * (pair_ptr->fit_log_adjust + sum);

#if 0
      printf("DEBUG\t%f\t%d\t%f\t%f\t%f\n",
             x, min_idx, xy_pairs2[min_idx].x, xy_pairs2[min_idx].y, pair_ptr->fit_log_adjust);
#endif
    }
  }


  if (xy_pairs)
    free(xy_pairs);

  if (xy_pairs2)
    free(xy_pairs2);

  return (0);
}


static double linear_interp(double x,
                            struct xy_pair *xy_pairs,
                            int last_idx,
                            int idx)
{
  double          a;

  assert(xy_pairs != NULL);

  if (idx >= 1 && xy_pairs[idx].x != xy_pairs[idx-1].x)
  {
    a = (xy_pairs[idx].x - x) / (xy_pairs[idx].x - xy_pairs[idx-1].x);

    return (a * xy_pairs[idx-1].y + (1.0 - a) * xy_pairs[idx].y);
  }

  /* use the value at idx */
  return (xy_pairs[idx].y);
}

