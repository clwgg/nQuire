#ifndef D_MODEL
#define D_MODEL

#include <math.h>
#include <stdint.h>

typedef struct {
  double *pos;
  uint8_t num;
} pos_f;

typedef struct {
  pos_f *a;
  uint32_t poscount;
  int basecount;
  float minfrac;
} all_f;

all_f * read_arr( FILE *file );
int flatten_arr( double *farr, all_f *arr );
void free_arr( all_f *arr );

//

static inline double prob_norm( double mean, double sd, double x )
{

  double norm;
  norm = (1/(sd * sqrt(2*M_PI))) * exp( -( pow(x - mean, 2) / (2 * pow(sd, 2)) ));

  return norm;
}

static inline double prob_unif( float minfrac, double x )
{

  double unif = 0.0;

  float upper = 1.0 - minfrac;

  if ( x >= minfrac && x <= upper ) {
    unif = 1 / (upper - minfrac);
  }

  return unif;
}

//

typedef struct {
  all_f *arr;
  double  a[4];
  double  m[3];
  double  w[3];
  int f; // bitflag: 1-est mix; 2-est mean; 4-est sd; 8-add unif;
  double ll;
} em_data;

// prototypes from em.c
void run_EM (em_data *data); 

void setup_EMstr (em_data *data, 
                  const double a[4], 
                  const double m[3], 
                  const double w[3], 
                  const int f);


#endif
