#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "d_model.h"

static int usage(char **argv)
{
  printf("\nUsage: %s %s file.bin\n\n", argv[-1], argv[0]);

  return 1;
}

int lrdmodel_main(int argc, char **argv)
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

  int f;
  // free model
  const double  a[4] = {0.2, 0.3, 0.5, 0.0};
  const double  m[3] = {0.5, 0.7, 0.2};
  const double  w[3] = {0.01, 0.01, 0.01};
  f = 0; f |= 7;

  em_data em;
  em.arr = arr;
  setup_EMstr(&em, a, m, w, f);
  run_EM(&em);
  double ll_free = em.ll;

  // fixed models
  f = 0; f |= 4;

  const double a1[4] = {1, 0, 0, 0};
  const double m1[3] = {0.5, 0, 0};
  setup_EMstr(&em, a1, m1, w, f);
  run_EM(&em);
  double ll_f1 = em.ll;

  const double a2[4] = {0.5, 0.5, 0, 0};
  const double m2[3] = {0.33, 0.67, 0};
  setup_EMstr(&em, a2, m2, w, f);
  run_EM(&em);
  double ll_f2 = em.ll;

  const double a3[4] = {0.33, 0.33, 0.33, 0};
  const double m3[3] = {0.25, 0.5, 0.75};
  setup_EMstr(&em, a3, m3, w, f);
  run_EM(&em);
  double ll_f3 = em.ll;

  fclose(file);
  free_arr(arr);

  printf("free\t%f\t   diff\n", ll_free);
  printf("dipl\t%f\t%f\n", ll_f1, ll_free - ll_f1);
  printf("trip\t%f\t%f\n", ll_f2, ll_free - ll_f2);
  printf("tetr\t%f\t%f\n", ll_f3, ll_free - ll_f3);

  return 0;

}
