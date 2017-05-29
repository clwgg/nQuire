#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "dump_utils.h"
#include "d_model.h"
#include "histo.h"

hist_arr build_frac_hist(hist_arr hist, int min)
{

  hist_arr frac_h = alloc_hist( hist.size, hist.offset );

  int i;
  for (i = 0; i < hist.size; ++i) {
    if (hist.y[i] != 0) {
      int s = hist.y[i] - min;
      frac_h.total += s;

      double f = s / hist.y[i];
      if (f < 0.01) {                   // avoid zero and negative probability
        frac_h.y[i] = 0.01;
      }
      else {
        frac_h.y[i] = f;
      }
    }
  }

  return frac_h;
}

static int usage(char **argv)
{
  printf("\nUsage: %s %s -o prefix file.bin\n\n", argv[-1], argv[0]);

  return 1;
}

int denoise_main(int argc, char **argv)
{
  
  int elem;
  char *out = NULL;

  while (( elem = getopt(argc, argv, "o:") ) >= 0) {
    switch(elem) {
      case 'o': out = optarg; break;
    }
  }

  if (argc - optind != 1) {
    return usage(argv);
  }
  else if (!out) {
    return usage(argv);
  }

  FILE *file = NULL;
  file = fopen(argv[optind], "rb+");

  if (!file) {
    printf("Can't open input file: %s\n\n", argv[optind]);
    return usage(argv);
  }

  size_t ret = 0;

  dump_h dhead;
  ret = read_header( &dhead, file );

  // allocate hist_arr
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

  // fill hist_arr
  while (1) {

    dump_p dpos;
    ret = read_pos( &dpos, file, &dhead);
    if (ret == 0) break;

    fill_histo(dpos.c, dpos.pos, &hist);
  }


  // rewind input and estimate uniform noise component to use for filtering
  ret = fseek(file, 0, SEEK_SET);

  all_f *arr = NULL;
  arr = read_arr(file);

  const double  a[4] = {0.2, 0.3, 0.4, 0.1};
  const double  m[3] = {0.5, 0.7, 0.2};
  const double  w[3] = {0.01, 0.01, 0.01};
  int f = 0; f |= 15;

  em_data em;
  em.arr = arr;

  setup_EMstr(&em, a, m, w, f);
  run_EM(&em);

  int min = arr->basecount * em.a[3] / hist.size;
  free_arr(arr);

  // establish filter fractions
  hist_arr frac = build_frac_hist(hist, min);
  free_hist(hist);

  // rewind input and write filtered output
  ret = fseek(file, 0, SEEK_SET);
  if (ret == 0) {

    ret = read_header( &dhead, file );
    int before = dhead.poscount;

    srand48(time(NULL));

    FILE *outf = NULL;
    char *outn = malloc( sizeof(char) * (strlen(out) + 10) );
    sprintf(outn, "%s.bin", out);
    outf = fopen(outn, "wb+");
    free(outn);

    if (!outf) {
      printf("Can't open output file.\n\n");
      free_hist(frac);
      fclose(file);
      return usage(argv);
    }

    dhead.poscount = 0;
    write_header(&dhead, outf);

    while (1) {
      dump_p dpos;
      ret = read_pos( &dpos, file, &dhead);
      if (ret == 0) break;

      int check = 0;

      int i;
      for(i = 0; i < dpos.num ; ++i) {
        double f = (double)dpos.pos[i] / dpos.c;
        int val = round( f * 100 );

        if (val >= frac.offset &&
            val <= 100 - frac.offset &&
            drand48() < frac.y[val - frac.offset])
          check = 1;
      }

      if (check) {
        write_pos(&dpos, outf, &dhead);
        ++dhead.poscount;
      }
    }

    ret = fseek(outf, 0, SEEK_SET);
    if (ret == 0) 
      write_header(&dhead, outf);
    else
      fprintf(stderr, "\n!Header corrupted in writing bindump!\n");

    fclose(outf);

    int after = dhead.poscount;

    printf("Positions to analyse:\n");
    printf("\t\t\tBefore: %d\n", before);
    printf("\t\t\tAfter:  %d (%.1f%%)\n", after, 100 * (double)after/(double)before);
  }
  else {
    printf("Error in reading input file - cannot fseek.\n");
  }

  free_hist(frac);
  fclose(file);
  return ret;

}

