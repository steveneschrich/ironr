/************************************************************************
 * irof
 * The irof (iterative rank-order filtering) function iteratively identifies
 * values from two parallel vectors that have large differences in rank.
 */
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "iron.h"


/* Calculate the maximum rank difference beyond which filtering occurs */
static int calculate_rankdiff_threshold(int *rank_differences, int nx, int *exclude,
                                 double min_rankdiff_percentile_threshold);

/* filter_iro
 * Iterative rank-order filtering of two vectors
 *
 * @description Iteratively filter out values from two parallel vectors
 * that have large differences in rank.
 *
 * @details The IRON algorithm selects a subset of values to create piecewise
 * linear scaling factors from. This function identifies the set of values to
 * use for changepoints/knees in the piecewise function. It does so by iteratively
 * filtering values from two parallel vectors that have the largest rank differences.
 *
 * IRO filters out values with large rank differences. Large rank differences
 * (threshold) is defined as 5% less than the maximum observed rank difference.
 * When this threshold is 1, filtering stops; or when the percentile of the
 * threshold (relative to all rank differences) is less than 1%.
 *
 * @param x Vector 1
 * @param y Vector 2
 * @param nx Length of vector(s)
 * @param exclude Logical (int) array indicating value should be excluded
 * @param min_rankdiff_percentile_threshold Minimum percentile of rankdiffs that
 * the calculated threshold must exceed.
 *
 * @return A dynamically allocated logical (int) array indicating if the
 * value is to be excluded.
 */
int *filter_iro(double *x, double *y, int nx, int *exclude,
            double min_rankdiff_percentile_threshold) {

  int i;

  /* Filtered is the active exclusion list, initialized from exclude parameter */
  int *filtered = (int *)calloc(nx, sizeof(int));
  if ( filtered == NULL ) return(NULL);
  if ( exclude != NULL ) {
    for ( i = 0; i < nx; i++ )
      filtered[i] = exclude[i];
  }

  /* Rank differences is used to hold the rank differences throughout */
  int *rank_diff = (int *)calloc(nx, sizeof(int));
  if ( rank_diff == NULL ) {
    free(filtered);
    return(NULL);
  }

  int rank_difference_threshold = 0;
  int num_excluded = 0;



  /* Iteratively identify differences in rank between x and ref
     above threshold and exclude them.
  */
  do {
    /* Rank remaining values */
    int *x_ranked = rank(x, nx, filtered);
    int *y_ranked = rank(y, nx, filtered);

    /* Calculate differences in ranks */
    for (i=0; i < nx; i++) {
      if ( filtered[i] ) continue;
      rank_diff[i] = labs(x_ranked[i] - y_ranked[i]);
    }

    /*
     * Calculate rank threshold. Note this threshold can be INT_MAX, in which
     * case nothing will be removed.
     */
    rank_difference_threshold = calculate_rankdiff_threshold(rank_diff, nx,
                                                             filtered, min_rankdiff_percentile_threshold);

    int num_filtered_before = ivec_length(rank_diff, nx, filtered);
    /* Add additional exclusions exceeding threshold */
    for (i = 0; i < nx; i++) {
      if ( filtered[i] ) continue;
      filtered[i] = rank_diff[i] > rank_difference_threshold;
    }
    int num_filtered_after = ivec_length(rank_diff, nx, filtered);
    num_excluded = num_filtered_before - num_filtered_after;
    /* Free the ranks that were calculated in this iteration */
    free(x_ranked);
    free(y_ranked);
  } while (num_excluded > 0);


  /* Free working objects */
  free(rank_diff);

  return(filtered);

}






/* calculate_rankdiff_threshold
 * Based on computed rank differences (rankdiffs), determine a threshold rankdiff
 * for filtering.
 *
 * @description Given a series of rank differences between two vectors, the IRON
 * algorithm selects a minimum rank difference to filter values on. The threshold
 * is 5% less than the percentile of the maximum rank difference.
 *
 * @details IRON works by iteratively calculating value ranks of two vectors and
 * then calculating the difference in these ranks. Differences larger than a threshold
 * are excluded from downstream piece-wise linear models. This function determines
 * the threshold dynamically. Generally the target is 5% less than the percentile
 * of the maximum rank difference. If the rank difference is 1 or if the threshold
 * percentile is below the minimum, a threshold of INT_MAX is returned.
 *
 * @note As described below, there is a deviation from the original IRON code here.
 * This function will return a threshold of INT_MAX in certain cases; this has the
 * effect of ensuring no rank difference filtering. However, in one edge case there
 * is a minimum percentile threshold (parameter min_rankdiff_precentile_threshold).
 * In the IRON C code, data will be filtered once when the that threshold is reached. In
 * this code, when the threshold is reached no filtering is performed (INT_MAX). This
 * simplifies the code and also seems more reasonable with respect to a threshold.
 *
 */
static int calculate_rankdiff_threshold(int *rankdiffs, int nx, int *exclude,
                                        double min_rankdiff_percentile_threshold) {

  /* Length of filtered vector (for calculations) */
  int nf = ivec_length(rankdiffs, nx, exclude);

  /* Maximum observed rank difference */
  int max_rank_difference = ivec_max(rankdiffs, nx, exclude);

  /* The threshold is defined as 5% less than the percentile of the max rank difference */
  double rankdiff_percentile_threshold = ((double) max_rank_difference / nf - 0.005);

  /*
   * Given the percentile threshold (rankdiff_percentile_threshold), back-calculate
   * the corresponding rank difference.
   * Return either the ceiling rank difference (integer value) or
   * INT_MAX if the threshold is too low.
   *
   * NOTE: This differs slightly from libaffy IRON. If rankdiff threshold is below 1.00005, then
   * libaffy backs up an iteration. Here, we use INT_MAX to avoid any filtering which should
   * not end up with functional differences. However, in the libaffy IRON version,
   * if the rankdiff_percentile_threshold is <= min_rankdiff_threshold + 1e-05 then filtering occurs
   * before stopping. The implication of including this condition in this function is that the final
   * filtering step (when rankdiff_percentile_threshold drops to the min threshold) is not
   * performed.
   */
  double rankdiff_threshold = nf * rankdiff_percentile_threshold;

  if ( rankdiff_threshold < 1.0 + 1E-5 )
    return(INT_MAX);
  if ( rankdiff_percentile_threshold < min_rankdiff_percentile_threshold)
    return(INT_MAX);

  return((int)(rankdiff_threshold + 0.5));

}


