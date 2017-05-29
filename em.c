#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <float.h>
#include <string.h>
#include <math.h>

#include "d_model.h"

static double comp_ll (em_data *data)
{
  double logsum = 0.0;
  int n,i,j;

  all_f *arr = data->arr;

  for (n = 0; n < arr->poscount; ++n) {
    pos_f p = arr->a[n];

    for (i = 0; i < p.num; ++i) {

      double freq = p.pos[i];

      double possum = 0.0;
      for( j = 0; j < 3; ++j ) {
        possum = possum + ( data->a[j] * prob_norm( data->m[j], data->w[j], freq ));
      }
      possum = possum + data->a[3] * prob_unif( arr->minfrac, freq );
      if (possum == 0.0) possum = DBL_MIN;
      logsum = logsum + log( possum );
    }
  }

  return logsum;
}

static void iter_EM (em_data *data, double *x)
{

  double tmp;

  all_f *arr = data->arr;
  int i, j;
    /* i: number of positions, j: number of normals in mixture */

  double g[3];
    /* stores each row g[j] in M[i][j] matrix while iterating i */
  memset(g, 0, sizeof(double) * 3 );

  double cs[3];
    /* calculates column sums cs[j] of M[i][j] */
  memset(cs, 0, sizeof(double) * 3 );
  double csx[3];
    /* calculates column sum of g[j] * x[i] */
  memset(csx, 0, sizeof(double) * 3 );
  double csd[3];
    /* calculates column sum of g[j] * (x[i] - m[j])^2 */
  memset(csd, 0, sizeof(double) * 3 );

  double sum = 0;
    /* calculates sum of all values in M[i][j] */

  for (i = 0; i < arr->basecount; ++i) {
    /*
     This loop is the E-step of EM and computes the sums of the
     posterior gamma(j) that will be neccessary for the
     calculations done in M.
     */
    double rs = 0;
      /* row sums of g[j] for normalization */

    for( j = 0; j < 3; ++j ) {
      g[j] = data->a[j] * ((tmp = prob_norm( data->m[j], data->w[j], x[i] )) == 0 ? DBL_MIN : tmp);
      rs = rs + g[j];
    }

    double ug;
      /* ug stores the positionwise posterior of the uniform mixture component*/
    ug = data->a[3] * prob_unif( arr->minfrac, x[i] );
      /* ug needs to be included in the row sum, it can be thought of as a fourth column */
    rs = rs + ug;

    for( j = 0; j < 3; ++j ) {
      g[j] = g[j] / rs;
      cs[j] = cs[j] + g[j];
      csx[j] = csx[j] + (g[j] * x[i]);
      csd[j] = csd[j] + (g[j] * pow(x[i] - data->m[j],2));
      sum = sum + g[j];
    }

    /* ug of course also is included in the total sum */
    sum = sum + (ug / rs);
  }

  for( j = 0; j < 3; ++j ) {
    /*
     This loop is the M-step of EM which gets new estimates 
     for a, m and w based on the different sums calculated in E.
     */

    double a = cs[j] / sum;
    double m = (1.0 / cs[j]) * csx[j];
    double w = sqrt( (1.0 / cs[j]) * csd[j] );

    if ((data->f >> 0) & 1 && !isnan(a))
      data->a[j] = a;

    if ((data->f >> 1) & 1 && !isnan(m))
      data->m[j] = m;

    if ((data->f >> 2) & 1 && !isnan(w))
      data->w[j] = w;
  }

  /* ua is the proportion of the uniform mixture that is included in the
   * em_data struct at the fourth position of the "a" array a[3]. */
  double ua = 1 - (data->a[0] + data->a[1] + data->a[2]);
  if ((data->f >> 3) & 1 && !isnan(ua))
    data->a[3] = ua;

}

void run_EM (em_data *data)
{

  double ll = comp_ll(data);
  double delta_ll = 1;
  double last = ll;

  all_f *arr = data->arr;
  double *x = malloc( sizeof(double) * arr->basecount );
  if (flatten_arr( x, arr ) == 0) {
    fprintf(stderr, "problem allocating Data in EM\n");
    return;
  }

  int c = 0;
  while ( delta_ll > 0.01 ) {

    iter_EM(data, x);

    if (c % 10 == 0) {
      last = ll;
      ll = comp_ll(data);
      delta_ll = ll - last;
    }

//    if (c % 100 == 0) fprintf(stderr, "EM iteration %d, delta loglik: %f\n", c, delta_ll);
    ++c;
  }

//  fprintf(stderr, "EM iteration %d, delta loglik: %f\n", c, delta_ll);

  free(x);
  data->ll = ll;
}

void setup_EMstr (em_data *data, 
                  const double a[4], 
                  const double m[3], 
                  const double w[3], 
                  const int f)
{
  memcpy(data->a, a, sizeof a[0] * 4);
  memcpy(data->m, m, sizeof m[0] * 3);
  memcpy(data->w, w, sizeof w[0] * 3);
  data->f = f;
}

