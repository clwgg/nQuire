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

int dmodel_main(int argc, char **argv)
{

  if (argc < 2) {
    return usage(argv);
  }
  char *in = argv[1];

  const double a1[4] = {1, 0, 0, 0};
  const double m1[3] = {0.5, 0, 0};

  const double a2[4] = {0.5, 0.5, 0, 0};
  const double m2[3] = {0.33, 0.67, 0};

  const double a3[4] = {0.33, 0.33, 0.33, 0};
  const double m3[3] = {0.25, 0.5, 0.75};

  const double  w[3] = {0.01, 0.01, 0.01};

  int f = 0; f |= 4;

  FILE *file = NULL;
  file = fopen(in, "rb+");

  if (!file) {
    printf("Can't open input file: %s\n\n", argv[1]);
    return usage(argv);
  }

  all_f *arr = NULL;
  arr = read_arr(file);

  em_data em;
  em.arr = arr;

  setup_EMstr(&em, a1, m1, w, f);
  run_EM(&em);
  printf("\nfinal diploid loglik: %f\n", em.ll);
  printf("     sd1 = %f\n", em.w[0]);

  setup_EMstr(&em, a2, m2, w, f);
  run_EM(&em);
  printf("\nfinal triploid loglik: %f\n", em.ll);
  printf("     sd1 = %f\n", em.w[0]);
  printf("     sd2 = %f\n", em.w[1]);

  setup_EMstr(&em, a3, m3, w, f);
  run_EM(&em);
  printf("\nfinal tetraploid loglik: %f\n", em.ll);
  printf("     sd1 = %f\n", em.w[0]);
  printf("     sd2 = %f\n", em.w[1]);
  printf("     sd3 = %f\n", em.w[2]);

  printf("\n");
  fclose(file);
  free_arr(arr);

  return 0;
}

