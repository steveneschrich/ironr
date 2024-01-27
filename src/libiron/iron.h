#ifndef IRON_H_
#define IRON_H_


#define MIN_SIGNAL       1E-5

/* Data structures
 *
 * There are three main data structures in the code:
 *   - xy_pair?
 *   - signal_pair
 *   - eqn_window
 */
struct xy_pair
{
  double x, y;
};


struct signal_pair
{
  int    index;
  double sig1;
  double sig2;
  int    rank1;
  int    rank2;
  int    rank_diff;
  char   initial_set_flag;
  char   irank_flag;

#if DEBUG_COLOR_IRANK
  double irank_frac_0;
  double irank_frac;
  double norm_err_scaled;
#endif

  double log_xy;
  double log_adjust;
  double fit_log_adjust;
  double norm_err;

  double weight;
  int n_windows;
};

struct eqn_window
{
  double slope;
  double offset;
  double start;
  double end;
};


/* Main entry points for library */

/* Filter parallel vectors by signal characteristics */
int *filter_signal(double *x, double *y, int nx, int *exclude,
            int exclude_16bit_saturation,
            int exclude_empirical_minimum,
            int exclude_duplicates);
/* Iterative rank order filtering */
int *filter_iro(double *x, double *y, int nx, int *exclude,
          double min_rankdiff_percentile_threshold);



/* Double vectors. The functions are named 'vec_'
 * to distinguish from integer vectors (see below).
 */
/* Return index of minimum value in vector */
int vec_imin(double *x, int nx, int *exclude);
/* Return minimum value in vector */
double vec_min(double *x, int nx, int *exclude);
/* Return index of maximum value in vector */
int vec_imax(double *x, int nx, int *exclude);
/* Return maximum value in vector */
double vec_max(double * x, int nx, int *exclude);
/* Return length of x, excluding masked values */
int vec_length(double *x, int nx, int *exclude);

/* Identify duplicated x,y combinations in the vectors */
int *duplicated(double *x,double *y, int nx, int *filtered);

/* Generate ranks for vector x */
int *rank(double *x, int nx, int *exclude);



/* Integer vectors. The functions are named 'ivec_' to
   distinguish for double ('vec_') functions.
 */
/* Return index of maximum value in vector */
int ivec_imax(int *x, int nx, int *exclude);
/* Return maximum value in vector */
int ivec_max(int *x, int nx, int * exclude);
/* Return length of x, excluding masked values */
int ivec_length(int *x, int nx, int *exclude);



#endif
