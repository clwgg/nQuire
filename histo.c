#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#include "histo.h"
#include "dump_utils.h"

void fill_histo( uint16_t n, uint16_t pos[6], hist_arr *hist )
{

  double *y = hist->y;
  int size = hist->size;
  int offset = hist->offset;

  int i;
  for(i = 0; i < 6; ++i) {
    if ( pos[i] != 0 ) {
      double frac = (double)pos[i] / n;
      int val = round( frac * 100 );
      if (0 + offset <= val && val < size + offset) {
        y[val - offset]++;
        hist->total++;
      }
    }
  }
}

hist_arr alloc_hist( int size, int offset )
{

  int i;
  hist_arr hist;
  hist.size = size;
  hist.total = 0.0;
  hist.offset = offset;

  hist.x = malloc( size * sizeof(double) );
  hist.y = malloc( size * sizeof(double) );

  memset( hist.y, 0, size * sizeof(double) );

  for (i = 0; i < size; ++i) {
    hist.x[i] = i + offset;
  }

  return hist;

}

void free_hist( hist_arr hist )
{
  free(hist.x); free(hist.y);
}

void ascii_hist( hist_arr *hist ) 
{
  int i;
  double *x = hist->x;
  double *y = hist->y;
  int size = hist->size;
  double total = hist->total;
  
  for(i = 0; i < size; ++i) {

    fprintf(stdout, "%d\t%.2g\t", (int)x[i], y[i]);

    double hn = (y[i] / total) * 1000;
    if (hn > y[i]) hn = y[i];

    int j;
    for(j = 0; j < hn; ++j) {
      fprintf(stdout, "#");
    }

    fprintf(stdout, "\n");
  }

}

////* Usage and main *////

static int usage(char **argv)
{
  printf("\nUsage: %s %s file.bin\n\n", argv[-1], argv[0]);

  return 1;
}

int histo_main(int argc, char **argv)
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

  ascii_hist( &hist );

  free_hist(hist);

  return 0;
}

