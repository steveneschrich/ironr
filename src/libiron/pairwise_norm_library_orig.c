
/**************************************************************************
 *
 * Filename:  pairwise_norm.c
 *
 * Purpose:   Pairwise normalization method.
 *
 * Creation:  10/07/10
 *
 * Author:    Eric Welsh
 *
 * Copyright: Copyright (C) 2010, Moffitt Cancer Center
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 10/07/10: Creation (EAW)
 * 01/10/13: Prevent iterative pruning from removing all points when N
 *           is small (EAW)
 * 01/10/13: Fixed rare case of reading beyond equation windows bounds,
 *           particularly when N points is small
 * 01/10/13: renamed linear normalization to global scaling, since that's
 *           actually what it is (true intensity-dependent linear
 *           linear normalization is not implemented yet)
 * 01/18/13: Pseudo-density weight exponent is now user definable.
 *           4 works well for microarray, 0 (unweighted) or at most 1 works
 *           best for proteomics.
 * 09/18/13: added iron_fit_both_x_y flag: better norm, may distrub rank order
 * 09/18/13: added iron_fit_window_frac (EAW)
 * 10/04/17: added iron_condense_training flag (EAW)
 * 06/01/18: added support for probeset exclusions during IRON training (EAW)
 * 10/11/19: print GlobalFitLine stats to stderr, similar to GlobalScale (EAW)
 *           stats are for the fit line, rather than the scaling adjustments
 *           (opposite direction/interpretation as GlobalScale)
 * 05/25/21: better handling of empty and near-empty samples (EAW)
 * 12/22/22: define M_PI if not already defined; for Ubuntu/Debian (EAW)
 *
 **************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#define MIN_SIGNAL       1E-5
#define DO_FLOOR         1

#if !defined(M_PI)
#define M_PI 3.14159265358979323846
#endif

/* (Unweighted): cisplatin breast subset works best with NO second pass,
 * rank fraction = 0.01, and window fraction = 0.05
 */

/* Window width fraction:
 *   0.05 too bumpy in UTSouthwestern Illumina 5058818012_E vs. 5067386018_F
 *   0.10 still a little bumpy, but much much smoother
 */


/*
 * flags for development and debug printing
 */
#define DEBUG_PRINT              1
#define DEBUG_FILE               0
#define DEBUG_COLOR_IRANK        0
#define DEBUG_DIE_EARLY          0
#define DEBUG_FIXED_RANK         0


static int compare_sort_xy_pair_by_x(const void *vptr1, const void *vptr2)
{
  struct xy_pair *pptr1 = (struct xy_pair *) vptr1;
  struct xy_pair *pptr2 = (struct xy_pair *) vptr2;

  assert(pptr1 != NULL);
  assert(pptr2 != NULL);

  if (pptr1->x < pptr2->x)
    return (-1);
  if (pptr1->x > pptr2->x)
    return (1);

  if (pptr1->y < pptr2->y)
    return (-1);
  if (pptr1->y > pptr2->y)
    return (1);

  return (0);
}


static int compare_sort_xy_pair_by_y(const void *vptr1, const void *vptr2)
{
  struct xy_pair *pptr1 = (struct xy_pair *) vptr1;
  struct xy_pair *pptr2 = (struct xy_pair *) vptr2;

  assert(pptr1 != NULL);
  assert(pptr2 != NULL);

  if (pptr1->y < pptr2->y)
    return (-1);
  if (pptr1->y > pptr2->y)
    return (1);

  if (pptr1->x < pptr2->x)
    return (-1);
  if (pptr1->x > pptr2->x)
    return (1);

  return (0);
}


static int compare_sort_xy_pair_by_x_plus_y(const void *vptr1,
                                            const void *vptr2)
{
  struct xy_pair *pptr1 = (struct xy_pair *) vptr1;
  struct xy_pair *pptr2 = (struct xy_pair *) vptr2;
  double val1, val2;

  assert(pptr1 != NULL);
  assert(pptr2 != NULL);

  val1 = pptr1->x + pptr1->y;
  val2 = pptr2->x + pptr2->y;

  if (val1 < val2)
    return (-1);
  if (val1 > val2)
    return (1);

  if (pptr1->y < pptr2->y)
    return (-1);
  if (pptr1->y > pptr2->y)
    return (1);

  if (pptr1->x < pptr2->x)
    return (-1);
  if (pptr1->x > pptr2->x)
    return (1);

  return (0);
}










/* return 1 on failure, 0 if no errors */
int fill_normalization_scales(char *filestem,
                               double *signals1,
                               double *signals2,
                               double *signals2_scales,
                               char *mask_array,
                               int num_spots,
                               double rank_frac_cutoff,
                               double rank_frac_cutoff2,
                               int condense_training_flag,
                               AFFY_COMBINED_FLAGS *flags,
                               double *return_training_frac,
                               double *return_rmsd)
{
  double min_sig1 = 9.0E8, min_sig2 = 9.0E8;
  int    bit16_flag1 = 1, bit16_flag2 = 1, i;

  struct signal_pair  *signal_pairs = NULL;
  struct signal_pair  *pair_ptr;
  struct signal_pair **tmp_ptrs1 = NULL, **tmp_ptrs2 = NULL;
  struct signal_pair **tmp_ptrs3 = NULL, **tmp_ptrs4 = NULL;
  struct signal_pair **old_filt1, **old_filt2, **filt1, **filt2;

  int old_num_filtered  = -42;
  int num_filtered      = 0;
  int num_unpruned      = 0;
  int orig_num_unpruned = 0;
  int num_not_weak      = 0;
  int num_both_not_weak = 0;

  int    max_rank_diff;
  int    old_rank_diff_cutoff = num_spots + 1;
  int    rank_diff_cutoff = num_spots + 1;
  int    start, end;
  double rank_diff_cutoff_frac = 999;
  double old_rank_diff_cutoff_frac = 999;

  struct eqn_window *eqn_windows = NULL;
  int    num_eqns = 0;
  double rmsd;
  double global_scale = 0.0;
  int    global_scaling_flag      = flags->iron_global_scaling_normalization;
  int    fit_both_x_y_flag        = flags->iron_fit_both_x_y;
  double weight_exponent          = flags->iron_weight_exponent;

#if DEBUG_FILE
  FILE *debug_file;
  double temp;
#endif

#if DEBUG_FIXED_RANK
  rank_frac_cutoff = 0.005;
#endif




  /* return array similarity metrics */
  rmsd = 0;
  for (i = 0; i < num_spots; i++)
    if (signal_pairs[i].initial_set_flag)
      rmsd += signal_pairs[i].fit_log_adjust * signal_pairs[i].fit_log_adjust;

  if (orig_num_unpruned)
    rmsd = sqrt(rmsd / orig_num_unpruned);

  *return_rmsd          = rmsd / log(10.0);
  *return_training_frac = (double) num_filtered / (double) orig_num_unpruned;

#if DEBUG_PRINT
  fprintf(stderr, "SimilarityMetrics:\tTrain\t%f\tRMSD\t%f\n",
     *return_training_frac, *return_rmsd);
#endif

  free(signal_pairs);
  free(temp_ptrs1);
  free(temp_ptrs2);
  free(temp_ptrs3);
  free(temp_ptrs4);
  free(eqn_windows);

#if DEBUG_DIE_EARLY
  exit(0);
#endif

  return (0);
}
