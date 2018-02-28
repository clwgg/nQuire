#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>

#include "htslib/htslib/sam.h"

#include "dump_utils.h"

typedef struct aux_t {
  samFile *in;
  bam_hdr_t *h;
  hts_idx_t *idx;
  hts_itr_t *itr;
  int min_mapQ;
} aux_t;

typedef struct {
  int pid;
  int start;
  int end;
  char name[256];
} bedline;

int read_bed(FILE *bed, bedline *line, aux_t *data)
{

  char *buf = NULL;
  size_t len = 0;
  int ret;

  memset(line, 0, sizeof *line);

  if ((ret = getline(&buf, &len, bed)) != -1) {
    buf = strtok(buf, "\n"); // remove trailing newline
    char *token = NULL;
    token = strtok(buf, " \t");
    if (!token) goto bederr;
    line->pid = bam_name2id(data->h, token);
    if (line->pid < 0) goto bederr;
    token = strtok(0, " \t");
    if (!token) goto bederr;
    line->start = atoi(token);
    token = strtok(0, " \t");
    if (!token) goto bederr;
    line->end = atoi(token);
    token = strtok(0, " \t");
    if (!token) goto bederr;
    strncpy(line->name, token, 255);
  }

  free(buf);
  return ret;

 bederr:
  free(buf);
  return -2;
}

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

int read_bam(void *data, bam1_t *b)
{

  aux_t *aux = (aux_t*)data;
  int ret;

  while (1) {
    ret = aux->itr ? sam_itr_next(aux->in, aux->itr, b) : sam_read1(aux->in, aux->h, b);
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
  printf("\t-c\tmin coverage of position to be reported (default: 10)\n");
  printf("\t-m\tmax coverage of position to be reported\n\n");
  printf("\t-x\tuse extended format for bindump\n\n");
  printf("\t-r\tbedfile for region based analysis (Format: Chr Start End Name)\n");
  printf("\t-y\tconcatenate bed regions into one file (only with -r)\n\n");

  return 1;

}

typedef struct {
  float minfrac;
  int het;
  int minc;
  int maxc;
  char *out;
  uint8_t flag;
  int bedcc;
  FILE *bed;
} c_opts;

int write_bam_iter(c_opts *opt, aux_t *data, bedline *line, FILE *file, dump_h *dhead)
{
  bam_plp_t iter;
  const bam_pileup1_t *plp;
  int tid, p, n;

  iter = bam_plp_init(read_bam, (void*)data);

  while ((plp = bam_plp_auto(iter, &tid, &p, &n)) != 0) {

    if (line) {
      if (p < line->start || p >= line->end)
        continue;
    }

    if (opt->minc && n < opt->minc) continue;
    if (opt->maxc && n > opt->maxc) continue;

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
      else if (pos[i] < (int)(opt->minfrac * n)) {
        pos[i] = 0;
      }
      else numbase++;
    }

    if (opt->het && numbase != 2) continue;
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

      write_pos(&dpos, file, dhead);

      ++dhead->poscount;
    }
  }

  bam_plp_destroy(iter);
  return 0;
}

int write_dump(c_opts *opt, aux_t *data, bedline *line)
{

  FILE *file = NULL;

  if (!line) {
    data->itr = NULL;

    char *outn = malloc( sizeof(char) * (strlen(opt->out) + 20) );
    if (opt->bedcc)
      sprintf(outn, "%s-bedcc.bin", opt->out);
    else
      sprintf(outn, "%s.bin", opt->out);
    file = fopen(outn, "wb");
    free(outn);
  } else {
    if (data->itr) {
      hts_itr_destroy(data->itr);
      data->itr = NULL;
    }
    data->itr = sam_itr_queryi(data->idx, line->pid, line->start, line->end);

    char *outn = malloc( sizeof(char) * (strlen(opt->out) + 10) * strlen(line->name) );
    sprintf(outn, "%s-%s.bin", opt->out, line->name);
    file = fopen(outn, "wb");
    free(outn);
  }

  if (!file) {
    printf("Can't open output file.\n\n");
    return -1;
  }

  dump_h dhead;
  dhead.poscount = 0;
  dhead.flag = opt->flag;
  dhead.minfrac = opt->minfrac;
  dhead.hhash = hash(data->h->text);

  write_header(&dhead, file);

  int ret = 0;
  if (opt->bedcc) {
    line = malloc(sizeof *line);
    while((ret = read_bed(opt->bed, line, data)) >= 0) {
      if (data->itr) {
        hts_itr_destroy(data->itr);
        data->itr = NULL;
      }
      data->itr = sam_itr_queryi(data->idx, line->pid, line->start, line->end);
      write_bam_iter(opt, data, line, file, &dhead);
    }
    free(line);
    if (ret == -2) { // Error handling for read_bed
      printf("Error while parsing bed file.\n");
      return -2;
    }
  }
  else {
    write_bam_iter(opt, data, line, file, &dhead);
  }

  ret = fseek(file, 0, SEEK_SET);
  if (ret == 0)
    write_header(&dhead, file);
  else
    fprintf(stderr, "\n!Header corrupted in creating bindump!\n");

  fclose(file);

  return 0;
}

int create_main(int argc, char **argv)
{

  char *bfile = NULL;
  char *bedf = NULL;
  int elem;
  int ret = 0;

  c_opts opt;
  opt.minfrac = 0.2;
  opt.het = 1;
  opt.minc = 10;
  opt.maxc = 0;
  opt.out = NULL;
  opt.flag = 0;
  opt.bedcc = 0;
  opt.bed = NULL;

  aux_t data;
  memset(&data, 0, sizeof data);
  data.min_mapQ = 1;

  bedline *line = NULL;

  while (( elem = getopt(argc, argv, "b:f:q:ac:m:o:xr:y") ) >= 0) {
    switch(elem) {
      case 'b': bfile = optarg; break;
      case 'f': opt.minfrac = atof(optarg); break;
      case 'q': data.min_mapQ = atoi(optarg); break;
      case 'a': opt.het = 0; break;
      case 'c': opt.minc = atoi(optarg); break;
      case 'm': opt.maxc = atoi(optarg); break;
      case 'o': opt.out = optarg; break;
      case 'x': opt.flag = 1; break;
      case 'r': bedf = optarg; break;
      case 'y': opt.bedcc = 1; break;
    }
  }

  if (!bfile || !opt.out) return usage(argv);

  if (!bedf && opt.bedcc) return usage(argv);

  data.in = sam_open(bfile, "r");

  if (!data.in) {
    ret = -1;
    goto err;
  }

  data.h = sam_hdr_read(data.in);

  if (bedf) {
    opt.bed = fopen(bedf, "r");
    if (!opt.bed) {
      printf("Could not open bedfile: %s\n", bedf);
      ret = -1;
      goto err;
    }

    data.idx = sam_index_load(data.in, bfile);
    if (!data.idx) {
      printf("Could not open bam index\n");
      ret = -1;
      goto err;
    }

    if (opt.bedcc) {
      ret = write_dump(&opt, &data, NULL);
      if (ret < 0)
        goto err;
    }
    else {
      line = malloc(sizeof *line);

      while((ret = read_bed(opt.bed, line, &data)) >= 0) {
        ret = write_dump(&opt, &data, line);
        if (ret < 0) { // Error handling for write_dump
          goto err;
        }
      }
      if (ret == -2) { // Error handling for read_bed
        printf("Error while parsing bed file.\n");
        goto err;
      }
      ret = 0;
    }

  }
  else {
    ret = write_dump(&opt, &data, NULL);
  }

 err:
  if(data.itr) hts_itr_destroy(data.itr);
  if(data.idx) hts_idx_destroy(data.idx);
  if(line) free(line);
  if(opt.bed) fclose(opt.bed);
  if(data.h) bam_hdr_destroy(data.h);
  if(data.in) sam_close(data.in);

  if (ret < 0)
    return usage(argv);

  return 0;
}
