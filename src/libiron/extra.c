static int compare_sort_sig1(const void *vptr1, const void *vptr2)
{
  struct signal_pair *pptr1 = *((struct signal_pair **) vptr1);
  struct signal_pair *pptr2 = *((struct signal_pair **) vptr2);

  if ( pptr1 == NULL || pptr2 == NULL) return(0);

  if (pptr1->sig1 < pptr2->sig1)
    return (-1);
  if (pptr1->sig1 > pptr2->sig1)
    return (1);

  if (pptr1->sig2 < pptr2->sig2)
    return (-1);
  if (pptr1->sig2 > pptr2->sig2)
    return (1);

  if (pptr1->index < pptr2->index)
    return (-1);
  if (pptr1->index > pptr2->index)
    return (1);

  return (0);
}


static int compare_sort_sig2(const void *vptr1, const void *vptr2)
{
  struct signal_pair *pptr1 = *((struct signal_pair **) vptr1);
  struct signal_pair *pptr2 = *((struct signal_pair **) vptr2);

  if ( pptr1 == NULL || pptr2 == NULL) return 0;

  if (pptr1->sig2 < pptr2->sig2)
    return (-1);
  if (pptr1->sig2 > pptr2->sig2)
    return (1);

  if (pptr1->sig1 < pptr2->sig1)
    return (-1);
  if (pptr1->sig1 > pptr2->sig1)
    return (1);

  if (pptr1->index < pptr2->index)
    return (-1);
  if (pptr1->index > pptr2->index)
    return (1);

  return (0);
}


// This business is to figure out if there is enough to make things
// work. Identical doesn't seem to be a problem, just a short circuit.
// But, the goal is to be identical to what I already have.
// So overall I would like this function to just return the indices or
// something of the chosen ones. And leave the rest to others.
for (i = 0; i < num_spots; i++)
{
  if (signals2[i] > MIN_SIGNAL)
  {
    num_not_weak++;

    if (signals1[i] > MIN_SIGNAL)
      num_both_not_weak++;
  }
}




void extra() {
  // NB: I think this should not be filter, for sure. It can be in the
  // driver function. But as a library element, at best is_identical can
  // be a function.
  /* check to see if we are normalizing vs. self
   * or if there are no good points
   */
  for (i = 0; i < num_spots; i++)
  {
    // Consider using identical here
    if (fabs(signals1[i] - signals2[i]) > 1E-5)
      break;
  }

  signal_pairs = calloc(num_spots, sizeof(struct signal_pair));
  if (signal_pairs == NULL)
  {
    fprintf(stderr, "ERROR -- calloc failed\n");

    return (1);
  }

  tmp_ptrs1 = calloc(num_spots, sizeof(struct signal_pair *));
  if (tmp_ptrs1 == NULL)
  {
    fprintf(stderr, "ERROR -- calloc failed\n");

    free(signal_pairs);

    return (1);
  }

  tmp_ptrs2 = calloc(num_spots, sizeof(struct signal_pair *));
  if (tmp_ptrs2 == NULL)
  {
    fprintf(stderr, "ERROR -- calloc failed\n");

    free(signal_pairs);
    free(temp_ptrs1);

    return (1);
  }

  tmp_ptrs3 = calloc(num_spots, sizeof(struct signal_pair *));
  if (tmp_ptrs3 == NULL)
  {
    fprintf(stderr, "ERROR -- calloc failed\n");

    free(signal_pairs);
    free(temp_ptrs1);
    free(temp_ptrs2);

    return (1);
  }

  tmp_ptrs4 = calloc(num_spots, sizeof(struct signal_pair *));
  if (tmp_ptrs4 == NULL)
  {
    fprintf(stderr, "ERROR -- calloc failed\n");


    free(signal_pairs);
    free(temp_ptrs1);
    free(temp_ptrs2);
    free(temp_ptrs3);

    return (1);
  }

  old_filt1 = tmp_ptrs1;
  old_filt2 = tmp_ptrs2;
  filt1     = tmp_ptrs3;
  filt2     = tmp_ptrs4;




  // This bit is the fail and return consistent but trivial state. That is,
  // scales =1 which does nothing to the data.
  // What I want to do here is assume (need to check) that both spots not
  // weak is a keeper. Or if one is weak, exclude it. Then we can fall through
  // the rest of the function and return an array with all excluded.
  if (i == num_spots || num_both_not_weak == 0)
  {
    for (i = 0; i < num_spots; i++)
      signals2_scales[i] = 1.0;

    *return_training_frac = 1.0;
    *return_rmsd = 0.0;

    if (global_scaling_flag)
      fprintf(stderr, "GlobalScale:\t%s\t%f\t%f\t%d\t%d\t%d\t%d\t%f\n",
              filestem,
              1.0, 0.0,
              num_not_weak,
              num_both_not_weak,
              num_not_weak, num_spots,
              1.0);
    else if (flags->iron_untilt_normalization)
      fprintf(stderr, "GlobalFitLine:\t%s\t%f\t%f\t%f\t%d\t%d\t%d\n",
              filestem,
              1.0, 0.0, 0.0,
              num_both_not_weak, num_not_weak, num_spots);

    return (0);
  }


  // In the driver function,
  // filter(detect_saturation=TRUE, exclude_empirical_minimum=TRUE)
  // if num_available == 0,
  // filter(detect_saturation=FALSE, exclude_empirical_minimum=FALSE).
  // This is only if you've turned them on I guess.

}
