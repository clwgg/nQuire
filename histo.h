#ifndef HISTO
#define HISTO

////* Histgram data *////

#define MAX_HIST 101

typedef struct hist_arr {

  double *x;
  double *y;
  int size;
  double total;
  int offset;

} hist_arr;

////* Histgram function prototypes *////

void fill_histo( uint16_t n, uint16_t pos[6], hist_arr *hist );
hist_arr alloc_hist( int size, int offset );
void free_hist( hist_arr hist );

#endif
