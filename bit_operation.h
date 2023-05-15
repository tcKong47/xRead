#ifndef BIT_OPERATION_H_
#define BIT_OPERATION_H_

#include "main.h"

int compare_rlen(const void *a, const void *b);
int compare_mi(const void *a, const void *b);
int compare_rep(const void *a, const void *b);
int compare_tid(const void *a, const void *b);
int compare_qs(const void *a, const void *b);
int compare_ove_len(const void *a, const void *b);
int compare_ove_score(const void *a, const void *b);
// int compare_ave_cov(const void *a, const void *b);
int compare_uint32_t(const void *a, const void *b);
int compare_double(const void *a, const void *b);

uint64_t hash64(uint64_t key, uint64_t mask);

extern uint8_t seq_nt4_table[256];

#endif
