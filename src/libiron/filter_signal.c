/**************************************************************************
 * filter
 *
 * This file is used to filter data based on signal
 * level. The goal is to arrive at a reasonable set of rows that are
 * suitable for applying a rank-based selection of invariants.
 */


#include <stdlib.h>
#include "iron.h"


/*
 * Static functions within this file (not exported).
 */

/* Detect 16-bit values (all less then 2^16) within vector */
static int is_16bit(double *x, int nx);

/***********************************************************************
 * filter_signal
 * Filter vector based on signal heuristics.
 *
 * Using two vectors (x and y), apply heuristics to indicate
 * specific values/positions to exclude from further processing. The result is
 * a logical vector indicating remove (TRUE) or keep (FALSE).
 *
 * @param nx The length of x and y.
 *
 * @param exclude A logical vector of same dimension as x/y which indicates
 * whether to exclude the value at the position from any consideration. It can
 * be NULL, in which case all values are considered. If it is not NULL, these
 * values are used as the starting point for filtering.
 *
 * @param exclude_16bit_saturation A very specific use case. If all values are less
 * than 2^16, then exclude values that are near 2^16 in case the data originated
 * from a 16 bit scanner and became saturation near the maximum value. This may
 * be useful in microarrays where scanners are used (and were at one time 16-bit
 * intensity scanners) but should generally not be needed.
 *
 * @param exclude_empirical_minimum If set, exclude the empirical minimum value
 * in the data (independently for both x and y).
 *
 * @param exclude_duplicates If set, the vectors will be filtered for duplicates
 * in both. This avoids an over-emphasis on repeated values during subsequent
 * training.
 */
int *filter_signal(double *x, double *y, int nx, int *exclude,
            int exclude_16bit_saturation,
            int exclude_empirical_minimum,
            int exclude_duplicates) {

  int i;

  /* Create a filter, initialized with exclude if passed as parameter */
  int *filtered = (int *)calloc(nx, sizeof(int));
  if ( filtered == NULL ) return(NULL);
  if ( exclude != NULL ) {
    for (i = 0; i < nx; i++) filtered[i]=exclude[i];
  }


  /* If requested, infer 16-bit saturation and exclude values near
  * the saturation point.
  */
  if ( exclude_16bit_saturation && (is_16bit(x, nx) || is_16bit(y,nx))) {
    for (i = 0; i < nx; i++) {
      if ( filtered[i] ) continue;
      if ( x[i] > 64000 || y[i] > 64000 ) filtered[i]=1;
    }
  }

  /* Exclude positions < MIN_SIGNAL */
  for (i = 0; i < nx; i++) {
    if ( filtered[i] ) continue;
    if (x[i] < MIN_SIGNAL || y[i] < MIN_SIGNAL) filtered[i]=1;
  }

  /* Separately exclude values == MIN_SIGNAL or at the empirical
   * minimum of remaining data (if flag set). Note the the empirical
   * minimum may be MIN_SIGNAL or higher, whereas the default exclusion
   * is MIN_SIGNAL.
   */
  double xmin = MIN_SIGNAL;
  double refmin = MIN_SIGNAL;
  if ( exclude_empirical_minimum ) {
    xmin = vec_min(x, nx, filtered);
    refmin = vec_min(y, nx, filtered);
  }
  for (i = 0; i < nx; i++) {
    if ( filtered[i] ) continue;
    if ( !(x[i] > xmin  && y[i] > refmin) ) filtered[i] = 1;
  }

  if ( exclude_duplicates ) {
    int *duplicates = duplicated(x, y, nx, NULL /*FIXME: filtered*/);
    if ( duplicates != NULL ) {
      for (i = 0; i < nx; i++)
        filtered[i] |= duplicates[i];
      free(duplicates);
    }
  }

  return(filtered);
}



/**
 * Infer data is from 16 bit representation (values within 2^16)
 */
static int is_16bit(double *x, int nx) {
  int i;
  for (i = 0; i < nx; i++) {
    if (x[i] > 1 << 16) return(0);
  }
  return(1);
}

