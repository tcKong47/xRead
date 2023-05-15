#ifndef INDEX_H
#define INDEX_H

#include "main.h"

mini_idx_t *indexing_input_read(bseq_file_t *fp, uint32_t *rep_n, uint32_t *n_sta, uint32_t *in_idx_n_sta, read_stat_t *read_stat, file_idx_t *file_idx, file_add_p_t *file_add_p);
uint32_t sketching_core(READ_t *seq, int ridx, MINIMIZER_t *mi, int w, int k);

#endif