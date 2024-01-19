#include <stddef.h>
#include <float.h>
#include <limits.h>

#include "iron.h"

/**
 * Identify position of first min value in vector, optionally
 * excluding values indicated in parallel logical array.
 */
int vec_imin(double *x, int nx, int *exclude) {
  int i;
  int min=-1;
  for (i = 0 ; i < nx; i++) {
    if ( exclude != NULL && exclude[i] ) continue;
    if ( min < 0 || x[i] < x[min] ) min = i;
  }
  return(min);
}

/**
 * Calculate minimum value of vector, optionally excluding values
 * indicated in parallel logical array.
 */
double vec_min(double *x, int nx, int *exclude) {
  int imin = vec_imin(x, nx, exclude);
  /* If there was no minimum (all excluded), use DBL_MIN as min. */
  if ( imin == -1 ) return(DBL_MIN);

  return(x[imin]);
}

/**
 * Identify position of first max value in vector, optionally
 * excluding values indicated in parallel logical array.
 */
int vec_imax(double *x, int nx, int *exclude) {
  int i;
  int max = -1;
  for (i = 0; i < nx; i++) {
    if ( exclude != NULL && exclude[i] ) continue;
    if ( max < 0 || x[i] > x[max] ) max = i;
  }
  return(max);
}



/**
 * Calculate maximum value of vector, optionally excluding values
 * indicated in parallel logical array.
 */
double vec_max(double *x, int nx, int *exclude) {
  int imax = vec_imax(x, nx, exclude);
  /* If there was no maximum (all excluded), use DBL_MAX as max. */
  if ( imax == -1 ) return(DBL_MAX);
  return(x[imax]);
}


/**
 * Identify position of first max value in vector, optionally
 * excluding values indicated in parallel logical array.
 */
int ivec_imax(int *x, int nx, int *exclude) {
  int i;
  int max = -1;
  for (i = 0; i < nx; i++) {
    if ( exclude != NULL && exclude[i] ) continue;
    if ( max < 0 || x[i] > x[max] ) max = i;
  }
  return(max);
}

/**
 * Calculate maximum value of vector, optionally excluding values
 * indicated in parallel logical array.
 */
int ivec_max(int *x, int nx, int *exclude) {
  int imax = ivec_imax(x, nx, exclude);
  /* If there was no maximum (all excluded), use INT_MAX as max. */
  if ( imax == -1 ) return(INT_MAX);
  return(x[imax]);
}

/**
 * Calculate the length of the vector excluding specific values.
 *
 * @param x Vector
 * @param nx Length of vector
 * @param exclude Logical vector of exclusions
 *
 * @return The length of x, excluding designated exclusions.
 */
int vec_length(double *x, int nx, int *exclude) {
  int i, length;

  for ( i = 0, length = 0; i < nx; i++ ) {
    if ( exclude != NULL && exclude[i] ) continue;
    length++;
  }
  return(length);
}

/**
 * Calculate the length of the vector excluding specific values.
 *
 * @param x Vector
 * @param nx Length of vector
 * @param exclude Logical vector of exclusions
 *
 * @return The length of x, excluding designated exclusions.
 */
int ivec_length(int *x, int nx, int *exclude) {
  int i, length;

  for ( i = 0, length = 0; i < nx; i++ ) {
    if ( exclude != NULL && exclude[i] ) continue;
    length++;
  }
  return(length);
}
