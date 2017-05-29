#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#include "dump_utils.h"
#include "d_model.h"


all_f * read_arr( FILE *file )
{
  dump_h dhead;
  read_header( &dhead, file );

  all_f *arr;
  arr = malloc( sizeof(all_f) );

  arr->a = malloc( sizeof(pos_f) * dhead.poscount );
  arr->poscount = dhead.poscount;
  arr->minfrac = dhead.minfrac;
  arr->basecount = 0;

  int i;
  for (i = 0; i < dhead.poscount; ++i) {

    dump_p dpos;
    read_pos( &dpos, file, &dhead );

    pos_f p;
    p.num = dpos.num;
    p.pos = malloc( sizeof(double) * dpos.num );

    int j;
    for (j = 0; j < dpos.num; ++j) {

      p.pos[j] = (double)dpos.pos[j] / dpos.c;
      ++arr->basecount;
    }

    arr->a[i] = p;

  }

  return arr;
}

int flatten_arr( double *farr, all_f *arr )
{
  int n = 0;
  int i,j;

  for (i = 0; i < arr->poscount; ++i) {
    pos_f p = arr->a[i];

    for (j = 0; j < p.num; ++j) {

      farr[n] = p.pos[j];
      ++n;
    }
  }

  if ( n == arr->basecount )
    return n;
  else
    return 0;
}

void free_arr( all_f *arr )
{
  int i;

  for (i = 0; i < arr->poscount; ++i) {
    free(arr->a[i].pos);
  }

  free(arr->a);
  free(arr);

}

