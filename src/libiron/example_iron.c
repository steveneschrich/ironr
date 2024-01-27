#include <stdlib.h>
#include <stdio.h>
#include "iron.h"
int example_iron(double *x, double *y, int nx) {

  /* Filter on various criteria */
  int *filtered_signal_values;
  int *filtered_iro_values;

  /* Filter based on signal */
  filtered_signal_values = filter_signal(x, y, nx, /* exclude */ NULL,
                                /* exclude_16bit_saturation */ 0,
                                /* exclude_empirical_minimum */ 1,
                                /* exclude_duplicates */ 1);

  if ( filtered_signal_values == NULL ) {
    fprintf(stderr, "Signal filtering was unsuccessful!");
    exit(-1);
  }
  // if ( length(x, nx, filtered_signal_values) == 0 )
  // if ( x == y)

  /* Filter based on rank differences */
  filtered_iro_values = filter_iro(x, y, nx,
                              /* exclude */ filtered_signal_values,
                              /* min_rankdiff_percentile_threshold */ 0.01);

  if ( filtered_iro_values == NULL ) {
    fprintf(stderr, "Iterative filtering was unsuccessful!");
    free(filtered_signal_values);
    exit(-1);
  }

  /* We should end up with a starting set that we can copy to a matrix or
   * something. Then free filtered_signal_values and filtered_irof_values */

  // piecewise_model = fit_plm(m);
  // apply(piecewise_model, y);

}
