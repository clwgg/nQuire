#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#include "histo.h"
#include "dump_utils.h"

////* Simulate normal Distribution *////

void sim_norm( double mean, double sd, hist_arr *hist )
{

  int i;
  double *x = hist->x;
  double *y = hist->y;
  int size = hist->size;

  for (i = 0; i < size; ++i) {
    double norm;
    norm = (1/(sd * sqrt(2*M_PI))) * exp( -( pow(x[i]-mean, 2) / (2 * pow(sd, 2)) ));

    y[i] = y[i] + norm;
    hist->total = hist->total + norm;
  }

}

////* Sum sqared Residuals *////

double norm_ssqr( hist_arr *hist, hist_arr *test )
/*
Calculate the sum of squared residuals on normalised y-values.
*/
{
  double ssqr = 0.0;

  double *y1 = hist->y;
  double t1 = hist->total;
  double *y2 = test->y;
  double t2 = test->total;
  int size = hist->size;

  if (hist->size == test->size) {

    int i;
    for(i = 0; i < size; ++i) {
      double norm_y1 = y1[i] / t1;
      double norm_y2 = y2[i] / t2;

      ssqr = ssqr + pow(norm_y1 - norm_y2, 2);
    }
  }

  return ssqr;
}

////* Simple linear regression *////

typedef struct reg_t {

  double alpha;
  double beta;
  double r;
  double beta_se;

} reg_t;

reg_t norm_yls( hist_arr *y, hist_arr *x )
/*
Perform simple linear least squares regression (as described for example in https://en.wikipedia.org/wiki/Simple_linear_regression)
on the normalised y[] values of two hist_arr structs.
*/
{

  reg_t result = {0};

  if (y->size == x->size) {

    int size = y->size;

    double mean = 1.0 / size;

    double sum_dist = 0.0;
    double sum_sqxd = 0.0;

    int i;
    for(i = 0; i < size; ++i) {

      double y_norm = y->y[i] / y->total;
      double x_norm = x->y[i] / x->total;

      sum_dist = sum_dist + ((x_norm - mean) * (y_norm - mean));
      sum_sqxd = sum_sqxd + pow(x_norm - mean, 2);

    }

    result.beta = sum_dist / sum_sqxd;
    result.alpha = mean - result.beta * mean;

    //

    double prod_tot = 0.0;
    double y_sq_tot = 0.0;
    double x_sq_tot = 0.0;

    double sum_err = 0.0;

    for(i = 0; i < size; ++i) {

      double y_norm = y->y[i] / y->total;
      double x_norm = x->y[i] / x->total;

      prod_tot = prod_tot + y_norm * x_norm;
      y_sq_tot = y_sq_tot + pow(y_norm, 2);
      x_sq_tot = x_sq_tot + pow(x_norm, 2);

      sum_err = sum_err + pow(y_norm - (result.alpha + result.beta * x_norm), 2);

    }

    double prod_m = prod_tot / size;
    double y_sq_m = y_sq_tot / size;
    double x_sq_m = x_sq_tot / size;

    result.r = ( prod_m - mean * mean ) / sqrt( (x_sq_m - pow(mean, 2)) * (y_sq_m - pow(mean, 2)) );

    result.beta_se = sqrt((1.0/(size - 2)) * sum_err / sum_sqxd);

  }

  return result;
}

////* Usage and main *////

static int usage(char **argv)
{
  printf("\nUsage: %s %s file.bin\n\n", argv[-1], argv[0]);

  return 1;

}

int htest_main(int argc, char **argv)
{

  if (argc < 2) {
    return usage(argv);
  }

  FILE *file = NULL;
  file = fopen(argv[1], "rb");

  if (!file) {
    printf("Can't open input file: %s\n\n", argv[1]);
    return usage(argv);
  }

  size_t ret = 0;

  dump_h dhead;
  ret = read_header( &dhead, file );
  float minfrac = dhead.minfrac;

  hist_arr hist;

  if (minfrac > 0) {
    int offset = round(minfrac * 100);
    int size = MAX_HIST - (2 * offset);

    hist = alloc_hist( size, offset );
  }
  else {
    hist = alloc_hist( MAX_HIST, 0);
  }

  while (1) {

    dump_p dpos;
    ret = read_pos( &dpos, file, &dhead);

    if (ret == 0) break;

    fill_histo(dpos.c, dpos.pos, &hist);

  }

  fclose(file);

  hist_arr dip, tri, tet;
  dip = alloc_hist( hist.size, hist.offset );
  tri = alloc_hist( hist.size, hist.offset );
  tet = alloc_hist( hist.size, hist.offset );

  sim_norm( 50.0, 5.0, &dip );

  sim_norm( 33.0, 4.0, &tri );
  sim_norm( 67.0, 4.0, &tri );

  sim_norm( 25.0, 4.0, &tet );
  sim_norm( 50.0, 5.0, &tet );
  sim_norm( 75.0, 4.0, &tet );

  reg_t reg;

  printf("Diploid:\n   Norm SSR: %g\n", norm_ssqr(&hist, &dip) );
  reg = norm_yls( &hist, &dip );
  printf("  y-y slope: %g, with std.Err: %g\n        r^2: %g\n\n", reg.beta, reg.beta_se, pow(reg.r, 2) );

  printf("Triploid:\n   Norm SSR: %g\n", norm_ssqr(&hist, &tri) );
  reg = norm_yls( &hist, &tri );
  printf("  y-y slope: %g, with std.Err: %g\n        r^2: %g\n\n", reg.beta, reg.beta_se, pow(reg.r, 2) );

  printf("Tetraploid:\n   Norm SSR: %g\n", norm_ssqr(&hist, &tet) );
  reg = norm_yls( &hist, &tet );
  printf("  y-y slope: %g, with std.Err: %g\n        r^2: %g\n\n", reg.beta, reg.beta_se, pow(reg.r, 2) );

  free_hist(dip);
  free_hist(tri);
  free_hist(tet);

  free_hist(hist);

  return 0;
}

