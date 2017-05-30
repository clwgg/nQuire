#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>

#include "htslib/htslib/sam.h"

#include "dump_utils.h"

// prototype hash function from create.c
uint32_t hash(char *str);

static int usage(char **argv)
{
  printf("\nUsage: %s %s [options] file.bin\n\n", argv[-1], argv[0]);
  printf("Options:\n");
  printf("\t-o\toutput prefix if new bindump should be created\n\n");
  printf("\t-c\tmin coverage of position to be reported\n");
  printf("\t-m\tmax coverage of position to be reported\n\n");
  printf("\t-f\tjust show flag of bindump\n");
  printf("\t-a\tbam file for annotation of chromosomes in extended bindump\n\n");

  return 1;
}

int view_main(int argc, char **argv)
{
  
  int elem;
  char *out = NULL;
  int minc = 0;
  int maxc = 0;
  int f = 0;
  char *bfile = NULL;

  int ex_code = 0;

  while (( elem = getopt(argc, argv, "c:m:o:fa:") ) >= 0) {
    switch(elem) {
      case 'c': minc = atoi(optarg); break;
      case 'm': maxc = atoi(optarg); break;
      case 'o': out = optarg; break;
      case 'f': f = 1; break;
      case 'a': bfile = optarg; break;
    }
  }

  if (argc - optind != 1) {
    return usage(argv);
  }

  samFile *b_in = NULL;
  bam_hdr_t *b_h = NULL;
  if (bfile) {
    b_in = sam_open(bfile, "r");
    if (!b_in) return usage(argv);
    b_h = sam_hdr_read(b_in);
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
  int flag = dhead.flag;

  if (f) {
    fprintf(stdout, "%d\n", flag);
    goto exit;
  }
  if (bfile && !flag) {
    fprintf(stdout, "Annotation with bam file only works for extended bindump!\n\n");
    fclose(file);
    bam_hdr_destroy(b_h);
    sam_close(b_in);
    return usage(argv);
  }
  if (flag && bfile) {
    uint32_t hhash2 = hash(b_h->text);
    if (dhead.hhash != hhash2) {
      fprintf(stdout, "\n\tHeader of bam file %s not identical to original source!\n\n", bfile);
      fclose(file);
      bam_hdr_destroy(b_h);
      sam_close(b_in);
      return 1;
    }
  }


  FILE *outf = NULL;
  if (out) {
    char *outn = malloc( sizeof(char) * (strlen(out) + 10) );
    sprintf(outn, "%s.bin", out);
    outf = fopen(outn, "wb+");
    free(outn);

    if (!outf) {
      printf("Can't open output file.\n\n");
      usage(argv);
      ex_code = 1;
      goto exit;
    }

    dhead.poscount = 0;
    write_header(&dhead, outf);

  }

  while (1) {

    dump_p dpos;
    ret = read_pos( &dpos, file, &dhead);

    if (ret == 0) break;

    char *seq = NULL;
    if (flag && bfile) seq = b_h->target_name[dpos.tid];

    if (minc && dpos.c < minc) continue;
    if (maxc && dpos.c > maxc) continue;

    if (out) {
      write_pos(&dpos, outf, &dhead);
      ++dhead.poscount;
    }
    else {

      if (seq) {
        fprintf(stdout, "%s\t%d\t", seq, dpos.p);
      }
      else if (flag) {
        fprintf(stdout, "%d\t%d\t", dpos.tid, dpos.p);
      }

      fprintf(stdout, "%d", dpos.c);
      int i;
      for (i = 0; i < dpos.num; ++i) {
        fprintf(stdout, "\t%d", dpos.pos[i]);
      }
      fprintf(stdout, "\n");

    }
  }

  if (out) {
    ret = fseek(outf, 0, SEEK_SET);
    if (ret == 0) 
      write_header(&dhead, outf);
    else
      fprintf(stderr, "\n!Header corrupted in writing bindump!\n");

    fclose(outf);
  }
exit:
  fclose(file);
  if (bfile) {
    bam_hdr_destroy(b_h);
    sam_close(b_in);
  }

  return ex_code;

}

