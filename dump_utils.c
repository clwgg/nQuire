#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "dump_utils.h"

#include <arpa/inet.h>

/* Reading */

size_t read_header(dump_h *h, FILE *file)
{
  size_t ret = 0;

  uint32_t poscount;
  uint8_t flag;
  uint32_t minfrac;

  ret = fread(&poscount, sizeof poscount, 1, file);
  ret = fread(&flag, sizeof flag, 1, file);
  ret = fread(&minfrac, sizeof minfrac, 1, file);

  h->poscount = ntohl(poscount);
  h->flag     = flag;
  h->minfrac  = (float)ntohl(minfrac) / 1000;

  if (h->flag) {
    uint32_t hhash;
    ret = fread(&hhash, sizeof hhash, 1, file);
    h->hhash = ntohl(hhash);
  }

  return ret;
}

size_t read_pos(dump_p *p, FILE *file, dump_h *h)
{
  size_t ret = 0;

  if (h->flag) {

    uint32_t tid;
    uint32_t pr;

    ret = fread(&tid, sizeof tid, 1, file);
    ret = fread(&pr, sizeof pr, 1, file);

    p->tid = ntohl(tid);
    p->p = ntohl(pr);
  }

  uint16_t c;
  uint8_t num;

  ret = fread(&c, sizeof c, 1, file);
  ret = fread(&num, sizeof num, 1, file);

  p->c = ntohs(c);
  p->num = num;

  if (ret == 0) // break early to avoid accessing undefined p->num
    return ret;

  uint16_t pos[6];
  memset(pos, 0, sizeof(pos));

  ret = fread(&pos, sizeof pos[0], p->num, file);

  memset(p->pos, 0, sizeof(p->pos));

  int i;
  for(i = 0; i < p->num; ++i)
    p->pos[i] = ntohs(pos[i]);

  return ret;
}

/* Writing */

void write_header(dump_h *h, FILE *file)
{

  uint32_t poscount = htonl(h->poscount);
  uint8_t flag      = h->flag;
  uint32_t minfrac  = (uint32_t)htonl(round(h->minfrac * 1000));
  uint32_t hhash    = htonl(h->hhash);

  fwrite(&poscount, sizeof poscount, 1, file);
  fwrite(&flag, sizeof flag, 1, file);
  fwrite(&minfrac, sizeof minfrac, 1, file);

  if (h->flag) 
    fwrite(&hhash, sizeof hhash, 1, file);

}

void write_pos(dump_p *p, FILE *file, dump_h *h)
{

  uint32_t tid = htonl(p->tid);
  uint32_t pr  = htonl(p->p);
  uint16_t c   = htons(p->c);
  uint8_t num  = p->num;

  if (h->flag) {
    fwrite(&tid, sizeof tid, 1, file);
    fwrite(&pr, sizeof pr, 1, file);
  }

  fwrite(&c, sizeof c, 1, file);
  fwrite(&num, sizeof num, 1, file);

  uint16_t pos[6];
  memset(pos, 0, sizeof(pos));

  int i;
  for(i = 0; i < p->num; ++i)
    pos[i] = htons(p->pos[i]);

  fwrite(&pos, sizeof pos[0], p->num, file);

}


