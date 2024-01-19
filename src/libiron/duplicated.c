#include <stdlib.h>


/* Structure to hold x,y and original index of row */
struct xypair {
  double x;
  double y;
  int index;
};

/* Compare xy by x, then y, then index */
static int xypair_compare(const void *x1, const void *x2);

/* duplicated
 * Determine if x,y values from parallel vectors x and y are
 * repeated.
 *
 * @description Given two parallel vectors x and y, identify
 * if the x,y combination is duplicated at multiple positions
 * within the vectors.
 *
 * @details This function mirrors the R equivalent, in which a
 * logical vector is constructed indicating whether values at
 * position i are duplicates of an earlier position.
 *
 * A logical vector (int) is constructed for the length of x/y (nx).
 * Values are true (non-zero) if the combination of x and y were seen
 * together at a position prior (duplicated).
 *
 * @note The first occurrence
 * of a pair of values would be FALSE (not duplicated) whereas the
 * subsequent occurrences would be TRUE (duplicated).
 *
 * @param x The first vector
 * @param y The second vector
 * @param nx The length of x and y
 * @param filtered An optional logical vector indicating positions to
 * exclude from comparisons.
 *
 * @return A dynamically allocated vector of logical (int) values,
 * indicating if the values that are duplicated.
 */
int *duplicated(double *x,double *y, int nx, int *filtered) {
  int i;

  /* NB: This isn't implemented yet. Not sure it's really needed */
  if ( filtered != NULL ) {
    return(NULL);
  }
  /* Allocate duplicate array */
  int *duplicates = (int *)calloc(nx, sizeof(int));
  if ( duplicates == NULL ) return(NULL);

  /* Load the vectors into a composite structure for sorting */
  struct xypair *p = (struct xypair *)calloc(nx, sizeof(struct xypair));
  if (p == NULL) {
    free(duplicates);
    return(NULL);
  }
  /* Load the xypair structure */
  for (i = 0; i < nx; i++) {
    p[i].x = x[i];
    p[i].y = y[i];
    p[i].index = i;
  }

  /* Sort the pairs by x, then y, then index. */
  qsort(p, nx, sizeof(struct xypair), xypair_compare);

  /* From sorted structure, check the current element against
   * the next element and flag the next element as duplicate
   * if it matches. Note that we flag the **next** element, so
   * that (due to sorting) the first occurence of the value is
   * not flagged as duplicate.
   */
  for (i = 0; i < nx-1; i++) {
    if ( p[i].x == p[i+1].x && p[i].y == p[i+1].y ) {
      int index = p[i+1].index;
      duplicates[index]=1;
    }
  }

  return(duplicates);

}

/* Sorting function
 * Given the xypair structure, sort on x, y, and then index.
 */
static int xypair_compare(const void *x1, const void *x2) {
  struct xypair *x1p = (struct xypair *)x1;
  struct xypair *x2p = (struct xypair *)x2;

  if ( x1p->x > x2p->x ) return 1;
  if ( x1p->x < x2p->x ) return -1;
  if ( x1p->y > x2p->y ) return 1;
  if ( x1p->y < x2p->y ) return -1;
  if ( x1p->index > x2p->index ) return 1;
  if ( x1p->index < x2p->index ) return -1;
  return 0;
}
