#ifndef MM_BSEQ_H
#define MM_BSEQ_H

#include <stdint.h>

#include "main.h"
#include "overlapping.h"

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t
{
    size_t l, m;
    char *s;
} kstring_t;
#endif

struct bseq_file_s;
typedef struct bseq_file_s bseq_file_t;

bseq_file_t *bseq_open(const char *fn);
void bseq_close(bseq_file_t *fp);
void bseq_seek_cur(bseq_file_t *fp, long int offset);
len_t *bseq_read_l(bseq_file_t *fp, uint64_t batch, uint32_t *_n);
READ_t *bseq_read(bseq_file_t *fp, uint64_t batch, uint32_t *_n);
READ_t * bseq_read_idx(bseq_file_t *fp, uint64_t batch, uint32_t *_n, uint32_t *_sta, uint32_t *_sta_idx, uint32_t *iter_idx, file_idx_t *file_idx);
READ_t *bseq_read_map(bseq_file_t *fp, uint64_t batch, uint32_t *_n, uint32_t *_sta, uint32_t *_sta_idx, uint32_t *iter_map, file_idx_t *file_idx);
int bseq_rewrite_read(bseq_file_t *fp, FILE *fp_un_ove, uint32_t *symn, file_idx_t *file_idx, char file_type);

#endif
