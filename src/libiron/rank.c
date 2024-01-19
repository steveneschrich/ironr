#include <stdlib.h>
struct ranked {
  double value;
  int original_order;
};
static int compare_ranked(const void *p1, const void *p2);

/* rank
 * Determine the rank order of the values in the vector.
 *
 * @description Calculate the rank order (1 is lowest) of the vector and
 * return a corresponding integer array of ranks.
 *
 * @param x Vector
 * @param nx Vector length
 * @param exclude Vector of indicators to exclude value in computation
 *
 * @return A dynamically allocated integer array of length |x| with integer
 * ranks corresponding to the values of x.
 */
int *rank(double *x, int nx, int *exclude) {
  int i;
  int nrank;

  struct ranked *sortedx = (struct ranked *)calloc(nx, sizeof(struct ranked));
  if ( sortedx == NULL ) return(NULL);

  int *ranks = (int *)calloc(nx, sizeof(int));
  if ( ranks == NULL ) {
    free(sortedx);
    return(NULL);
  }

  /* Load sortedx into structure and sort it. Note that the length of
     the object can be shorter than nx if the exclude vector has
     exclusions.
  */
  for (i = 0, nrank = 0; i < nx; i++) {
    if ( exclude != NULL && exclude[i] ) continue;
    sortedx[nrank].value = x[i];
    sortedx[nrank].original_order = i;
    nrank++;
  }
  qsort(sortedx, nrank, sizeof(struct ranked), compare_ranked);

  /* Sorted index is the rank, store it in the original (unsorted) index position */
  /* NB: This is using 1-based ranking, not sure if this is consistent with original
   * code
   */
  for (i = 0; i < nrank; i++) {
    int original_order = sortedx[i].original_order;
    ranks[original_order] = i+1;
  }

  /* Free sorted data after using it */
  free(sortedx);

  return(ranks);
}


/* compare_ranked
 * Compare two ranked items based on their value.
 */
static int compare_ranked(const void *p1, const void *p2) {
  struct ranked *d1 = (struct ranked *)p1;
  struct ranked *d2 = (struct ranked *)p2;

  if ( d1->value > d2->value ) return(1);
  if ( d1->value < d2->value ) return(-1);
  return(0);
}
