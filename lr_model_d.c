#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <pthread.h>

#include "htslib/htslib/thread_pool.h"

#include "d_model.h"

typedef struct {
  double free;
  double f1;
  double f2;
  double f3;
  all_f *arr;
  char *file;
} ll_out;

typedef struct {
  hts_tpool *p;
  hts_tpool_process *q;
  int n;
  char **files;
} t_optd;

static int usage(char **argv)
{
  printf("\nUsage: %s %s [-t threads] file1.bin [file2.bin ...]\n\n", argv[-1], argv[0]);

  return 1;
}

void *lrdmodel_driver(void *arg)
{
  ll_out *out = (ll_out *)arg;

  if (!out)
    return out;

  all_f *arr = out->arr;

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
  out->free = em.ll;

  // fixed models
  f = 0; f |= 4;

  const double a1[4] = {1, 0, 0, 0};
  const double m1[3] = {0.5, 0, 0};
  setup_EMstr(&em, a1, m1, w, f);
  run_EM(&em);
  out->f1 = em.ll;

  const double a2[4] = {0.5, 0.5, 0, 0};
  const double m2[3] = {0.33, 0.67, 0};
  setup_EMstr(&em, a2, m2, w, f);
  run_EM(&em);
  out->f2 = em.ll;

  const double a3[4] = {0.33, 0.33, 0.33, 0};
  const double m3[3] = {0.25, 0.5, 0.75};
  setup_EMstr(&em, a3, m3, w, f);
  run_EM(&em);
  out->f3 = em.ll;

  return out;
}

static void *lrdmodel_dispatch(void *arg)
{
  t_optd *optd = (t_optd *)arg;

  int i;

  for (i = 0; i < optd->n; i++) {

    FILE *file = NULL;
    file = fopen(optd->files[i], "rb");

    if (!file) {
      printf("Can't open input file: %s\n\n", optd->files[i]);
      // Send NULL pointer as termination signal
      hts_tpool_dispatch(optd->p, optd->q, lrdmodel_driver, NULL);
      pthread_exit((void *)1);
    }

    all_f *arr = NULL;
    arr = read_arr(file);

    fclose(file);

    ll_out *out = NULL;
    out = malloc(sizeof *out);
    memset(out, 0, sizeof *out);
    out->arr = arr;
    out->file = optd->files[i];

    hts_tpool_dispatch(optd->p, optd->q, lrdmodel_driver, out);
  }

  // Send NULL pointer as termination signal
  hts_tpool_dispatch(optd->p, optd->q, lrdmodel_driver, NULL);
  pthread_exit(NULL);

}


int lrdmodel_main(int argc, char **argv)
{

  int n = 1;
  int elem;

  while (( elem = getopt(argc, argv, "t:") ) >= 0) {
    switch(elem) {
    case 't': n = atoi(optarg); break;
    }
  }

  if (argc - optind == 0) {
    return usage(argv);
  }

  printf("file\tfree\tdip\ttri\ttet\td_dip\td_tri\td_tet\n");

  hts_tpool *p = hts_tpool_init(n);
  hts_tpool_process *q = hts_tpool_process_init(p, n*2, 0);
  t_optd optd = {p, q, argc-optind, argv+optind};
  pthread_t tid;

  pthread_create(&tid, NULL, lrdmodel_dispatch, &optd);

  for(;;) {
    hts_tpool_result *r = hts_tpool_next_result_wait(q);

    ll_out *out = (ll_out *)hts_tpool_result_data(r);

    if (!out) {
      hts_tpool_delete_result(r, 1);
      break;
    }

    printf("%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
           out->file, out->free,
           out->f1, out->f2, out->f3,
           out->free - out->f1, out->free - out->f2, out->free - out->f3);

    free_arr(out->arr);
    hts_tpool_delete_result(r, 1);
  }

  hts_tpool_process_flush(q);

  hts_tpool_process_destroy(q);
  hts_tpool_destroy(p);
  pthread_join(tid, NULL);

  return 0;
}
