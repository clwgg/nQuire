#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>

#include "htslib/htslib/sam.h"

#include "dump_utils.h"

uint32_t hash(char *str)
/*
  djb2 hash after Dan Bernstein, see http://www.cse.yorku.ca/~oz/hash.html
*/
{
    uint32_t hash = 5381;
    int c;

    while ((c = *str++))
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */

    return hash;
}


typedef struct aux_t {

  samFile *in;
  bam_hdr_t *h;
  int min_mapQ;

} aux_t;

int read_bam(void *data, bam1_t *b)
{

  aux_t *aux = (aux_t*)data;
  int ret;

  while (1) {
    ret = sam_read1(aux->in, aux->h, b);
    if (ret < 0) break;
    if ( (int)b->core.qual < aux->min_mapQ ) continue;
    if ( b->core.flag & (BAM_FUNMAP | BAM_FDUP) ) continue;
    break;
  }
  return ret;

}

void count_pos( int n, const bam_pileup1_t *plp, uint16_t pos[6] )
{

  int i;
  for (i = 0; i < n; ++i) {

    const bam_pileup1_t *p = plp + i;
    uint8_t *seq = bam_get_seq(p->b);
    uint8_t nuc = bam_seqi(seq, p->qpos);

    if (p->is_del) nuc = 32;

    switch (nuc) {
      case 1:  pos[0]++ ; break;
      case 2:  pos[1]++ ; break;
      case 4:  pos[2]++ ; break;
      case 8:  pos[3]++ ; break;
      case 15: pos[4]++ ; break;
      case 32: pos[5]++ ; break;
    }
  }

}

static inline void netw_sort6( uint16_t pos[6] )
{
#define min(x, y) (x<y?x:y)
#define max(x, y) (x<y?y:x) 
#define SWAP(x,y) { const int a = min(pos[x], pos[y]); \
                    const int b = max(pos[x], pos[y]); \
                    pos[x] = a; pos[y] = b; }
    SWAP(1, 2);
    SWAP(4, 5);
    SWAP(0, 2);
    SWAP(3, 5);
    SWAP(0, 1);
    SWAP(3, 4);
    SWAP(1, 4);
    SWAP(0, 3);
    SWAP(2, 5);
    SWAP(1, 3);
    SWAP(2, 4);
    SWAP(2, 3);
#undef SWAP
#undef min
#undef max
}

static int usage(char **argv)
{
  printf("\nUsage: %s %s [options] -b file.bam -o prefix\n\n", argv[-1], argv[0]);
  printf("Options:\n");
  printf("\t-b\tInput BAM file\n");
  printf("\t-o\toutfile prefix\n\n");
  printf("\t-a\treport all positions, instead of just biallelic ones (not recommended for ploidy analysis)\n");
  printf("\t-f\tmin fraction of read support (default: 0.2)\n");
  printf("\t-q\tmin map quality (default: 1)\n");
  printf("\t-c\tmin coverage of position to be reported\n");
  printf("\t-m\tmax coverage of position to be reported\n\n");
  printf("\t-x\tuse extended format for bindump\n\n");

  return 1;

}

int create_main(int argc, char **argv)
{

  char *bfile = NULL;
  float minfrac = 0.2;
  int het = 1;
  int minc = 0;
  int maxc = 0;
  int elem;
  char *out = NULL;

  uint8_t flag = 0;

  bam_plp_t iter;
  const bam_pileup1_t *plp;
  int tid, p, n;

  aux_t data;
  data.min_mapQ = 1;

  while (( elem = getopt(argc, argv, "b:f:q:ac:m:o:x") ) >= 0) {
    switch(elem) {
      case 'b': bfile = optarg; break;
      case 'f': minfrac = atof(optarg); break;
      case 'q': data.min_mapQ = atoi(optarg); break;
      case 'a': het = 0; break;
      case 'c': minc = atoi(optarg); break;
      case 'm': maxc = atoi(optarg); break;
      case 'o': out = optarg; break;
      case 'x': flag = 1; break;
    }
  }

  if (!bfile || !out) return usage(argv);

  data.in = sam_open(bfile, "r");

  if (!data.in) return usage(argv);

  data.h = sam_hdr_read(data.in);

  iter = bam_plp_init(read_bam, (void*)&data);

  FILE *file = NULL;

  char *outn = malloc( sizeof(char) * (strlen(out) + 10) );
  sprintf(outn, "%s.bin", out);
  file = fopen(outn, "wb+");
  free(outn);

  if (!file) {
    printf("Can't open output file.\n\n");
    return usage(argv);
  }

  dump_h dhead;
  dhead.poscount = 0;
  dhead.flag = flag;
  dhead.minfrac = minfrac;
  dhead.hhash = hash(data.h->text);

  write_header(&dhead, file);

  while ((plp = bam_plp_auto(iter, &tid, &p, &n)) != 0) {

    if (minc && n < minc) continue;
    if (maxc && n > maxc) continue;

    if ( n > UINT16_MAX ) {
      fprintf(stderr, 
          "Position %d has coverage of %d, which is higher than max of %d... skipping\n",
          p, n, UINT16_MAX);
      continue;
    }

    uint16_t pos[6];
    memset(pos, 0, sizeof(pos));

    count_pos( n, plp, pos );
    netw_sort6(pos);

    uint8_t numbase = 0;
    int i;
    for (i = 0; i < 6; ++i) {
      if (pos[i] == 0) continue;
      else if (pos[i] < (int)(minfrac * n)) {
        pos[i] = 0;
      }
      else numbase++;
    }

    if (het && numbase != 2) continue; 
    else {

      dump_p dpos;
      dpos.tid = (uint32_t)tid;
      dpos.p = (uint32_t)p;
      dpos.c = (uint16_t)n;
      dpos.num = numbase;

      memset(dpos.pos, 0, sizeof(dpos.pos));
      for (i = 0; i < numbase; ++i) {
        dpos.pos[i] = pos[6 - (numbase - i)];
      }

      write_pos(&dpos, file, &dhead);

      ++dhead.poscount;
    }
  }

  int ret = fseek(file, 0, SEEK_SET);
  if (ret == 0) 
    write_header(&dhead, file);
  else
    fprintf(stderr, "\n!Header corrupted in creating bindump!\n");

  bam_plp_destroy(iter);
  bam_hdr_destroy(data.h);
  sam_close(data.in);

  fclose(file);

  return 0;
}

