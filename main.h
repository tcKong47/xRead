/*************************************************************************
	> File Name: 
	> Author: 
	> Mail: 
 ************************************************************************/

#ifndef _ASSEMBLY_H
#define _ASSEMBLY_H

#include <stdio.h>
#include <stdint.h>
#include <getopt.h> 
#include <zlib.h>

typedef struct
{
	float x_read; //x% longest reads
	int k; //k-mer length
	int w; //window size
	int l; //l-mer length
	float r_n; //repeat n threshood
	int s_s;
	int m_b; //matching bases
	int min_ove;
    int max_hang;
	int top_n; // keep top n overlaps
	float cov_ratio;
	int part_size;
	int X_read; //x% longest reads
	int split_len;
	int iter;
	int trans_iter;
    int trans_file;
    char *output_file;
	uint64_t batch_size;
    int read_type;
	int thread_n;
	int memory;
}param_map;

typedef struct READ
{
	uint32_t rid; //true read id
	char* read_name;
	char* read_seq;
	uint32_t read_length;
}READ_t;

typedef struct MINIMIZER
{
	uint64_t x; // minimizer hash value
	uint64_t y; // rid<<32 | pos<<1 | strand
}MINIMIZER_t;

typedef struct
{
	uint64_t mi_m, mi_n;
	MINIMIZER_t *mi_v;
} mini_t;

typedef struct
{
	int k, w, l;
	uint32_t bucket_num;
	uint32_t *mi_count;
	double mi_step;
	mini_t mm;
	uint64_t mi_mask;
} mini_idx_t;

typedef struct
{
	uint64_t n, m;
	int *offs;
} file_idx_t;


extern float x_read;
extern int k;
extern int w;
extern int l;
extern float r_n;
extern int s_s;

extern int m_b;
extern int min_ove;
extern int max_hang;
extern int top_n;

extern float cov_ratio;
extern int part_size;
extern int X_read;
extern int split_len;
extern int iter;
extern int trans_iter;
extern int trans_file;

extern int waitingLen;
extern uint64_t batch_size_base;
extern char *output_file;
extern int thread_n;
extern int memory;
extern double realtime0;

#endif