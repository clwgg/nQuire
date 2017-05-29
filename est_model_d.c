#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "d_model.h"

////* Usage and main *////

static int usage(char **argv)
{
  printf("\nUsage: %s %s file.bin\n\n", argv[-1], argv[0]);

  return 1;
}

int est_dmodel_main(int argc, char **argv)
{

  if (argc < 2) {
    return usage(argv);
  }
  char *in = argv[1];

  FILE *file = NULL;
  file = fopen(in, "rb+");

  if (!file) {
    printf("Can't open input file: %s\n\n", argv[1]);
    return usage(argv);
  }

  all_f *arr = NULL;
  arr = read_arr(file);

  const double  a[4] = {0.2, 0.3, 0.5, 0.0};
  const double  m[3] = {0.5, 0.7, 0.2};
  const double  w[3] = {0.01, 0.01, 0.01};
  int f = 0; f |= 7;
  em_data em;
  em.arr = arr;

  setup_EMstr(&em, a, m, w, f);
  run_EM(&em);
  printf("\nfinal loglik: %f\n\n", em.ll);
  int j;
  for( j = 0; j < 3; ++j ) {
    printf("      a = %0.10g\n", em.a[j]);
    printf("      m = %0.10g\n", em.m[j]);
    printf("     sd = %0.10g\n\n", em.w[j]);
  }

  fclose(file);
  free_arr(arr);

  return 0;
}

