#ifndef DUMP_UTILS
#define DUMP_UTILS

#include <stdint.h>

typedef struct {
  uint32_t poscount;
  uint8_t flag;
  float minfrac;
  uint32_t hhash;
} dump_h;

typedef struct {
  uint32_t tid;
  uint32_t p;
  uint16_t c;
  uint8_t num;
  uint16_t pos[6];
} dump_p;


size_t read_header(dump_h *h, FILE *file);
size_t read_pos(dump_p *p, FILE *file, dump_h *h);
void write_header(dump_h *h, FILE *file);
void write_pos(dump_p *p, FILE *file, dump_h *h);

#endif
