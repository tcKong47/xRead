#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <assert.h>
#include <math.h>

#include "paf.h"
#include "bseq.h"
#include "ktime.h"
#include "index.h"
#include "graph.h"
#include "thread.h"
#include "overlapping.h"
#include "bit_operation.h"

#define INTER_RESULT
// #define LOAD_INTER_RESULT

// shared data of all threads
typedef struct
{
	uint32_t map_batch;
	uint32_t idx_batch;
	file_idx_t file_idx;
	file_add_p_t file_add_p;

	int iter;
	read_stat_t *read_stat;
	uint32_t max_len;
	uint32_t x_len;
	uint32_t rep_n;
	uint32_t batch_i;
	uint32_t n_threads;
	uint32_t n_processed;
	uint32_t n_total_read;
	uint64_t index_batch_size;
	uint64_t query_batch_size;
	bseq_file_t *fp;
	mini_idx_t *mi;
	uint32_t top_n;
	uint32_t index_n_sta;     // for iter_idx, start rid of current batch
	uint32_t index_idx_n_sta; // for iter_idx, start index of current batch
	uint32_t index_idx_n_end; // for iter_idx, end index of current batch, start index of next batch
	uint32_t map_n_sta;     // for iter_map, start rid of current batch
	uint32_t map_idx_n_sta; // for iter_map, start index of current batch
	uint32_t map_idx_n_end; // for iter_map, end index of current batch, start index of next batch
	OVE_C_t *ove_cl;
	uint32_t *symm;
	int m_b;
	int min_ove;
	int max_hang;

	uint32_t search_step;
	uint32_t split_len;
	read_cov_t *read_cov;
	double *ave_cov;
	uint16_t *blkn;
	double ave_median;
	uint32_t **visited_rid;
} pipeline_t;

typedef struct
{
	uint32_t n_seq;
	len_t *len;
} step_sort_t;

typedef struct
{
	const pipeline_t *p;
	uint32_t n_seq;
	uint32_t n_sta;
	uint32_t n_sta_idx;
	READ_t *read_info;
	thread_buf_t **buf;
	uint32_t ove_alloc;
} step_query_t;

static void *sort_read_length_pipeline(void *shared, int step, void *in)
{
	int i;
	pipeline_t *p = (pipeline_t *)shared;
	if (step == 0)
	{
		step_sort_t *s = NULL;
		s = (step_sort_t *)calloc(1, sizeof(step_sort_t));
		if (s == NULL)
		{
			fprintf(stderr, "[%s] calloc s memory error!\n", __func__);
			exit(1);
		}
		s->len = bseq_read_l(p->fp, p->index_batch_size, &s->n_seq);
		if (s->n_seq)
		{
			return s;
		}
		else
		{
			if (s != NULL) {free(s); s = NULL;}
		}
	}
	else if (step == 1)
	{
		step_sort_t *s = (step_sort_t *)in;
		if (s->n_seq)
		{
			p->file_idx.offs = (int *)realloc(p->file_idx.offs, (p->n_processed + s->n_seq) * sizeof(int));
			p->read_stat->rlen = (uint32_t *)realloc(p->read_stat->rlen, (p->n_processed + s->n_seq) * sizeof(uint32_t));
			if (p->read_stat->rlen == NULL)
			{
				fprintf(stderr, "[%s] calloc %ldGB p->read_stat->rlen memory error!\n", __func__, (p->n_processed + s->n_seq) * sizeof(uint32_t) / 1024 / 1024 / 1024);
				exit(1);
			}
			for (i = 0; i < s->n_seq; i++)
			{
				p->read_stat->rlen[p->n_processed++] = s->len[i].rlen;
				if (s->len[i].rlen > p->max_len) p->max_len = s->len[i].rlen;
			}
			for (i = 0; i < s->n_seq; i++)
			{
				long int offset = 0;
				offset += s->len[i].rlen * sizeof(char);
				offset += (s->len[i].nlen) * sizeof(char);
				p->file_idx.offs[p->file_idx.n++] = offset;
			}
		}
		qsort(s->len, s->n_seq, sizeof(len_t), compare_rlen);
		uint32_t read_n = s->n_seq * x_read * 0.01;
		p->x_len += s->len[read_n - 1].rlen;
		p->batch_i += 1;
		if (s->len != NULL) {free(s->len); s->len = NULL;}
		if (s != NULL) {free(s); s = NULL;}
		// fprintf(stderr, "[2%s: %.3fs, %.3fGBG]\tReal time:%.3f sec\tCPU:%.3f sec\tMemory peak:%.3f GB\n", __func__, realtime() - realtime0, peak_memory() / 1024.0 / 1024.0, realtime() - realtime0, peak_memory() / 1024.0 / 1024.0);
	}
	return 0;
}

static int binary_search(MINIMIZER_t minmer, int k, int l, uint64_t *pos, mini_idx_t *mi)
{
	uint64_t key;
	uint64_t pos_s = 0, pos_e = 0, pos_m = 0;
	uint64_t mask_l = (1ULL << 2 * l) - 1; // (111...111) 20

	pos[0] = pos[1] = -1;
	key = (minmer.x >> (k-l)*2) & mask_l;
	pos_s = mi->mi_count[key];
	pos_e = mi->mi_count[key+1] - 1;

	if (minmer.x < mi->mm.mi_v[pos_s].x || minmer.x > mi->mm.mi_v[pos_e].x)	return -1;
	while(pos_s <= pos_e)
	{
		pos_m = (pos_s + pos_e)/2;
		if (mi->mm.mi_v[pos_m].x == minmer.x)
		{
			pos[0] = pos[1] = pos_m;
			//left extend
			int sub_s = pos_s, sub_e = pos_m - 1, sub_m;
			while(sub_s <= sub_e)
			{
				sub_m = (sub_s + sub_e) / 2;
				if (mi->mm.mi_v[sub_m].x == minmer.x)
				{
					pos[0] = sub_m;
					sub_e = sub_m - 1;
				}
				else if (mi->mm.mi_v[sub_m].x < minmer.x)	sub_s = sub_m + 1;
				else if (mi->mm.mi_v[sub_m].x > minmer.x)	exit(1);
			}
			//right extend
			sub_s = pos_m + 1;
			sub_e = pos_e;
			while(sub_s <= sub_e)
			{
				sub_m = (sub_s + sub_e) / 2;
				if (mi->mm.mi_v[sub_m].x == minmer.x)
				{
					pos[1] = sub_m;
					sub_s = sub_m + 1;
				}
				else if (mi->mm.mi_v[sub_m].x > minmer.x)	sub_e = sub_m -1 ;
				else if (mi->mm.mi_v[sub_m].x < minmer.x)	exit(1);
			}
			return 1;
		}
		else if (mi->mm.mi_v[pos_m].x > minmer.x)	pos_e = pos_m - 1;
		else if (mi->mm.mi_v[pos_m].x < minmer.x)	pos_s = pos_m + 1;
	}
	return -1;
}

static uint32_t merge_colinear_mr(MR_t *vertex_mr, uint32_t *v_n)
{
	uint32_t su_i = 0, j, tid_tmp;
	uint32_t s1, e1, Eindel = 20, tstr_tmp, qstr_tmp, cov_tmp;
	tid_tmp = vertex_mr[0].t_id;
	tstr_tmp = vertex_mr[0].tstr;
	qstr_tmp = vertex_mr[0].qstr;
	j = 0;
	while (j < (*v_n))
	{
		s1 = j++;
		cov_tmp = vertex_mr[s1].cov;
		while ((tid_tmp == vertex_mr[j].t_id) && ((tstr_tmp ^ vertex_mr[j].tstr) == (qstr_tmp ^ vertex_mr[j].qstr)) && (vertex_mr[j].qs >= vertex_mr[j - 1].qs) && (j < (*v_n)))
		{
			int diff = (int)(vertex_mr[j].qs - vertex_mr[j - 1].qe - 1);
			int diff_q = (int)(vertex_mr[j].qs - vertex_mr[j - 1].qs);
			int diff_tt = vertex_mr[j].tstr == vertex_mr[j].qstr ? (int)(vertex_mr[j].ts - vertex_mr[j - 1].ts) : (int)(vertex_mr[j - 1].te - vertex_mr[j].te);
			// printf("s1 %d,%d %d, diff %d %d, %d, %d, %d\n",s1,j-1,j,diff,waitingLen,diff_q,diff_tt,abs(diff_tt - diff_q));
			if (diff > waitingLen)
				break;
			if (abs(diff_tt - diff_q) < Eindel)
			{
				cov_tmp += (diff > 0) ? vertex_mr[j++].cov : (diff + vertex_mr[j++].cov);
			}
			else
				break;
		}
		e1 = j - 1;

		vertex_mr[su_i].t_id = vertex_mr[s1].t_id;
		vertex_mr[su_i].qs = vertex_mr[s1].qs;
		vertex_mr[su_i].qe = vertex_mr[e1].qe;
		vertex_mr[su_i].ts = vertex_mr[s1].tstr == vertex_mr[s1].qstr ? vertex_mr[s1].ts : vertex_mr[e1].ts;
		vertex_mr[su_i].te = vertex_mr[e1].tstr == vertex_mr[e1].qstr ? vertex_mr[e1].te : vertex_mr[s1].te;
		vertex_mr[su_i].tstr = vertex_mr[e1].tstr;
		vertex_mr[su_i].qstr = vertex_mr[e1].qstr;
		vertex_mr[su_i++].cov = cov_tmp;

		tid_tmp = vertex_mr[j].t_id;
		tstr_tmp = vertex_mr[j].tstr;
		qstr_tmp = vertex_mr[j].qstr;
		// if (su_i >= v_m) {printf("[%s] use unlegal memory, %d-%d, something wrong with it...", __func__, su_i, v_m); exit(0);}
	}
	(*v_n) = su_i;
	return su_i;
}

static uint32_t finding_hits(void *shared_data, int64_t qidx, uint32_t pid, mini_t *seq_mi)
{
	int result;
	step_query_t *data = (step_query_t *)shared_data;
	READ_t *query = &data->read_info[qidx];
	const pipeline_t *p = data->p;
	mini_idx_t *mi = data->p->mi;
	uint32_t v_n = 0, *v_m = &data->buf[pid]->v_m;
	uint32_t i, ii, hits_num = 0;
	uint64_t range[2] = {UINT64_MAX, UINT64_MAX};
	uint32_t pos_q, pos_tt, strand_mask = 1ULL, pos_mask = (1ULL << 32) - 1;
	uint32_t su_i = 0, j, tid_tmp, cov_tmp;
	// uint32_t s1, e1, Eindel = 20, tstr_tmp, qstr_tmp;
	uint32_t rid;
	int k_i;

	assert(query->read_length <= seq_mi->mi_m);
	seq_mi->mi_n = sketching_core(&data->read_info[qidx], qidx, seq_mi->mi_v, mi->w, mi->k);
	if (seq_mi->mi_n >= seq_mi->mi_m) printf("[%s] use unlegal memory, %ld-%ld, something wrong with it...", __func__, seq_mi->mi_n, seq_mi->mi_m);

	for (i = 0; i < seq_mi->mi_n; ++i)
	{
		// rep_mini_m++;
		result = binary_search(seq_mi->mi_v[i], p->mi->k, p->mi->l, range, mi);
		if (result != -1)
		{
			hits_num = range[1] - range[0] + 1;
			if (hits_num > 0 && hits_num <= p->rep_n) //meeting repeat region, drop these hits
			{
				for (ii = range[0]; ii <= range[1]; ++ii)
				{
					rid = p->read_stat->iter_idx[mi->mm.mi_v[ii].y >> 32];
					if (rid != query->rid)
					{
						pos_q = ((seq_mi->mi_v[i].y & pos_mask) >> 1);
						pos_tt = ((mi->mm.mi_v[ii].y & pos_mask) >> 1);

						if (v_n >= (*v_m))
						{
							(*v_m) = (*v_m) == 0 ? 256 : (*v_m) << 1;
							data->buf[pid]->vertex_mr = (MR_t *)realloc(data->buf[pid]->vertex_mr, (*v_m) * sizeof(MR_t));
							if (data->buf[pid]->vertex_mr == NULL)
							{
								fprintf(stderr, "[%s] calloc %ldGB data->buf[pid]->vertex_mr memory error!\n", __func__, (*v_m) * sizeof(MR_t) / 1024 / 1024 / 1024);
								exit(1);
							}
							fprintf(stderr, "[%s] will not use these code, something wrong with it...",__func__);
						}
						data->buf[pid]->vertex_mr[v_n].qs = pos_q + 1 - k;
						data->buf[pid]->vertex_mr[v_n].qe = pos_q;
						data->buf[pid]->vertex_mr[v_n].ts = pos_tt + 1 - k;
						data->buf[pid]->vertex_mr[v_n].te = pos_tt;
						data->buf[pid]->vertex_mr[v_n].cov = k;
						data->buf[pid]->vertex_mr[v_n].t_id = mi->mm.mi_v[ii].y >> 32;
						data->buf[pid]->vertex_mr[v_n].tstr = (mi->mm.mi_v[ii].y & strand_mask);
						data->buf[pid]->vertex_mr[v_n++].qstr = (seq_mi->mi_v[i].y & strand_mask);
					}
				}
			}
		}
	}
	
	// if ((double)rep_mini_n / seq_mi->mi_n > 0.8) data->buf[pid]->is_ove[qidx] = 1;
	data->buf[pid]->v_n = v_n;
	if (v_n == 0)	return v_n;

	// merge co-linear match region in the same target read
	qsort(data->buf[pid]->vertex_mr, v_n, sizeof(MR_t), compare_tid);

	su_i = merge_colinear_mr(data->buf[pid]->vertex_mr, &v_n);
	data->buf[pid]->v_n = su_i;

	uint32_t count = 1, su_ii = 0;
	cov_tmp = data->buf[pid]->vertex_mr[0].cov;
	tid_tmp = data->buf[pid]->vertex_mr[0].t_id;
	for (j = 1; j < su_i; ++j)
	{
		if (data->buf[pid]->vertex_mr[j].t_id != tid_tmp)
		{
			if (count == 1 && cov_tmp <= k)
			{
				count = 1;
				cov_tmp = data->buf[pid]->vertex_mr[j].cov;
				tid_tmp = data->buf[pid]->vertex_mr[j].t_id;
			}
			else
			{
				// printf("%d %d %d %d\n",count,cov_tmp,tid_tmp,j);
				for (k_i = count; k_i > 0; k_i--)
				{
					// printf("%d %d %d %d\n", su_ii, j, k_i, j-k_i);
					data->buf[pid]->vertex_mr[su_ii++] = data->buf[pid]->vertex_mr[j - k_i];
				}
				count = 1;
				cov_tmp = data->buf[pid]->vertex_mr[j].cov;
				tid_tmp = data->buf[pid]->vertex_mr[j].t_id;
			}
		}
		else
		{
			cov_tmp += data->buf[pid]->vertex_mr[j].cov;
			count++;
		}
		if (su_ii >= data->buf[pid]->v_m) printf("[%s] use unlegal memory, %d-%d, something wrong with it...", __func__, su_ii, data->buf[pid]->v_m);
	}

	// the last one
	if (!(count == 1 && cov_tmp <= k))
	{
		for (k_i = count; k_i > 0; k_i--)
		{
			data->buf[pid]->vertex_mr[su_ii++] = data->buf[pid]->vertex_mr[j - k_i];
		}
	}

	if (su_ii >= data->buf[pid]->v_m) printf("[%s] use unlegal memory, %d-%d, something wrong with it...", __func__, su_ii, data->buf[pid]->v_m);
	data->buf[pid]->v_n = su_ii;

	// detect repetitive regions
	tid_tmp = data->buf[pid]->vertex_mr[0].t_id;
	uint32_t qs_tmp = data->buf[pid]->vertex_mr[0].qs;
	uint32_t qe_tmp = data->buf[pid]->vertex_mr[0].qe;
	uint32_t vi, rep_mini_n = 0, rep_mini_m = 0, split_len = 1000, cur_window_pos = split_len;
	for (vi = 1; vi < su_ii; ++vi)
	{
		if (data->buf[pid]->vertex_mr[j].t_id != tid_tmp)
		{
			if ((double)rep_mini_n / rep_mini_m > 0.8) {data->buf[pid]->is_ove[qidx] = 1; break;}
			tid_tmp = data->buf[pid]->vertex_mr[j].t_id;
			rep_mini_n = rep_mini_m = 0;
		}
		else
		{
			rep_mini_m++;
			if (data->buf[pid]->vertex_mr[j].qs == qs_tmp && data->buf[pid]->vertex_mr[j].qe == qe_tmp) hits_num++;
			else
			{
				if (hits_num >= 2) rep_mini_n++;
				hits_num = 0;
				qs_tmp = data->buf[pid]->vertex_mr[j].qs;
				qe_tmp = data->buf[pid]->vertex_mr[j].qe;
			}
			if (qe_tmp >= cur_window_pos)
			{
				if ((double)rep_mini_n / rep_mini_m > 0.8) {data->buf[pid]->is_ove[qidx] = 1; break;}
				cur_window_pos += split_len;
				rep_mini_n = rep_mini_m = 0;
			}
		}
	}
	// the last one
	if ((double)rep_mini_n / rep_mini_m > 0.8) {data->buf[pid]->is_ove[qidx] = 1;}

	// if (query->rid == 9184)
	// {
	// 	printf("read %d\t\n", query->rid);
	// 	int vi;
	// 	printf("qid\thid\tcov\tqstrand\tqs\tqe\ttid\ttstrand\tts\tte\tlen\n");
	// 	for (vi = 0; vi < su_ii; ++vi)
	// 	{
	// 		// if (data->buf[pid]->vertex_mr[vi].t_id == 92)
	// 		printf("%ld\t%d\t%d\t%d\t%ld\t%ld\t%d\t%d\t%ld\t%ld\n", qidx, vi, data->buf[pid]->vertex_mr[vi].cov, data->buf[pid]->vertex_mr[vi].qstr, data->buf[pid]->vertex_mr[vi].qs, data->buf[pid]->vertex_mr[vi].qe, data->buf[pid]->vertex_mr[vi].t_id, data->buf[pid]->vertex_mr[vi].tstr, data->buf[pid]->vertex_mr[vi].ts, data->buf[pid]->vertex_mr[vi].te);
	// 	}
	// 	printf("\n");
	// }

	return v_n;
}

static void print_single_read_overlap(OVE_t ove)
{
	printf("qid %d\tql %d, tid %d tl %d\n", ove.qid, ove.ql, ove.tid, ove.tl);
	printf("qs qe %d %d; ts te %d %d\n", ove.qs , ove.qe, ove.ts, ove.te);
	printf("rev %d\t mbp %d\t mln %d\t score %.6f\n", ove.rev, ove.mbp, ove.mln, ove.score);
}

static int get_overlap_info(step_query_t *data, thread_buf_t *buf, uint32_t q_idx, uint32_t pid, uint32_t max_index_n, uint32_t *max_index, PATH_t *dist_path, mini_t *seq_mi)
{
	READ_t *query = &data->read_info[q_idx];
	int32_t j, k;
	uint32_t i = 0, matching_bases = 0;
	OVE_t ove_tmp;
	OVE_t *ove_cl_tmp = NULL;
	uint32_t n_tmp = 0, n_rep = 0;
	OVE_t *ove_rep = NULL;
	int32_t *tmp_path_idx = NULL, idx_n = 0;
	path_idx_t *path_idx = NULL;
	int32_t sta_idx, end_idx, tid_tmp, rev_tmp, diff_q, diff_tt;
	uint32_t qpre = 0, qsuf = 0, tpre = 0, tsuf = 0;
	uint32_t left_overhang = 0, right_overhang = 0, overhang_tr = 0, overhang = 0;
	MR_t *vertex_mr = buf->vertex_mr;
	int mln_ratio, mln_th;
	// uint32_t tmp_tpos;

	if (max_index_n == 0) return 0;
	ove_cl_tmp = (OVE_t*)calloc(max_index_n,sizeof(OVE_t));
	ove_rep = (OVE_t*)calloc(max_index_n,sizeof(OVE_t));
	tmp_path_idx = (int32_t*)calloc(buf->v_n,sizeof(int32_t));
	path_idx = (path_idx_t*)calloc(max_index_n,sizeof(path_idx_t));
	if (ove_cl_tmp == NULL)
	{
		fprintf(stderr, "[%s] calloc %ldGB ove_cl_tmp memory error!\n", __func__, max_index_n * sizeof(OVE_t) / 1024 / 1024 / 1024);
		exit(1);
	}

	for (i = 0; i < max_index_n; ++i)
	{
		j = (int32_t)max_index[i];
		matching_bases = 0;
		idx_n = 0;
		ove_tmp.score = dist_path[j].dist;
		// first iter vertex_mr.t_id is the rid of index read; second iter vertex_mr.t_id is the index of index read in array read_stat->iter_idx
		ove_tmp.tid  = vertex_mr[j].t_id;
		ove_tmp.qid = query->rid;
		ove_tmp.ql = query->read_length;
		ove_tmp.tl = data->p->read_stat->rlen[data->p->read_stat->iter_idx[ove_tmp.tid]];
		ove_tmp.rev = vertex_mr[j].tstr == vertex_mr[j].qstr ? 0 : 1;
		ove_tmp.qe = vertex_mr[j].qe;
		if (ove_tmp.rev == 0)
			ove_tmp.te = vertex_mr[j].te;
		else
			ove_tmp.ts = vertex_mr[j].ts;
		// ove_tmp.te = ove_tmp.rev == 0 ? vertex_mr[j].te : vertex_mr[j].ts;

		while(j != -1)
		{
			tmp_path_idx[idx_n++] = j;
			ove_tmp.qs = vertex_mr[j].qs;
			if (ove_tmp.rev == 0)
				ove_tmp.ts = vertex_mr[j].ts;
			else
				ove_tmp.te = vertex_mr[j].te;
			// ove_tmp.ts = ove_tmp.rev == 0 ? vertex_mr[j].ts : vertex_mr[j].te;
			matching_bases += vertex_mr[j].cov;
			// printf("path %d\t%d\t%d\t%d\t%ld\t%ld\t%d\t%d\t%ld\t%ld\n", query->rid, j, vertex_mr[j].cov, vertex_mr[j].qstr, vertex_mr[j].qs, vertex_mr[j].qe, vertex_mr[j].t_id, vertex_mr[j].tstr, vertex_mr[j].ts, vertex_mr[j].te);
			j = (int32_t)dist_path[j].pre_node;
			if (j >= (int32_t)buf->v_n) printf("[%s] 1 use unlegal memory, %d-%d, something wrong with it...", __func__, j, buf->v_n);
		}

		qpre = ove_tmp.qs; qsuf = ove_tmp.ql - ove_tmp.qe;
		if (ove_tmp.rev == 0)
			{tpre = ove_tmp.ts; tsuf = ove_tmp.tl - ove_tmp.te;}
		else if (ove_tmp.rev == 1)
			{tpre = ove_tmp.tl - ove_tmp.te; tsuf = ove_tmp.ts;}

		left_overhang  = qpre < tpre ? qpre : tpre;
		right_overhang = qsuf < tsuf ? qsuf : tsuf;
		overhang = left_overhang + right_overhang;
		ove_tmp.mbp = matching_bases;
		ove_tmp.mln = ove_tmp.te - ove_tmp.ts + 1;

		overhang_tr = data->p->max_hang < ove_tmp.mln * 0.8 ? data->p->max_hang : ove_tmp.mln * 0.8;
		mln_th = ove_tmp.ql < ove_tmp.tl ? ove_tmp.ql : ove_tmp.tl;
		mln_ratio = data->p->min_ove > mln_th * 0.3 ? data->p->min_ove : mln_th * 0.3;
		if (matching_bases >= data->p->m_b && overhang < overhang_tr && ove_tmp.mln > mln_ratio)
		{
			ove_cl_tmp[n_tmp++] = ove_tmp;
			if (n_tmp > max_index_n) printf("[%s]2 use unlegal memory, %d-%d, something wrong with it...", __func__, n_tmp, max_index_n);
		}
		else if (ove_tmp.score > data->p->m_b)
		{
			path_idx[n_rep].idx_n = idx_n;
			path_idx[n_rep].idx = (int32_t *)calloc(idx_n, sizeof(int32_t));
			for (k = 0; k < path_idx[n_rep].idx_n; ++k)
			{
				path_idx[n_rep].idx[k] = tmp_path_idx[k];
			}
			ove_rep[n_rep++] = ove_tmp;
		}
		// if (ove_tmp.qid == 9184)
		// 	print_single_read_overlap(ove_tmp);
	}
	int Eindel = 100;
	// int result;
	// const pipeline_t *p = data->p;
	// mini_idx_t *mi = data->p->mi;
	// uint64_t range[2] = {UINT64_MAX, UINT64_MAX};
	// uint32_t m_i, ii, t_id;
	// uint32_t pos_q, pos_tt, pos_ts, pos_qs, strand_mask = 1ULL, pos_mask = (1ULL << 32) - 1;
	// uint8_t tstr, qstr;
	// int add_flag = 0;
	if (n_tmp < data->p->top_n && n_rep > n_tmp)
	{
		// if (ove_tmp.qid == 82902 || ove_tmp.qid == 453976 || ove_tmp.qid == 604087 || ove_tmp.qid == 203541 || ove_tmp.qid == 1647 || ove_tmp.qid == 2076 || ove_tmp.qid == 3921 || ove_tmp.qid == 15250)
		// printf("[%s] Processing read %d with repeat regions, n_tmp %d ...\n", __func__, query->rid, n_tmp);
		// qsort(ove_rep, n_rep, sizeof(OVE_t), compare_ove_score);
		for (i = 0; i < n_rep; ++i)
		{
			sta_idx = path_idx[i].idx[path_idx[i].idx_n-1];
			end_idx = path_idx[i].idx[0];
			tid_tmp = vertex_mr[sta_idx].t_id;
			rev_tmp = vertex_mr[sta_idx].qstr ^ vertex_mr[sta_idx].tstr;
			// printf("sta idx %d, end idx %d, tid %d\n", sta_idx, end_idx, tid_tmp);
			uint32_t v_n = 0, v_m = 64;
			MR_t *path_mr = (MR_t *)calloc(v_m, sizeof(MR_t));
			// for (m_i = 0; m_i < seq_mi->mi_n; ++m_i)
			// {
			// 	result = binary_search(seq_mi->mi_v[m_i], p->mi->k, p->mi->l, range, mi);
			// 	if (result != -1)
			// 	{
			// 		for (ii = range[0]; ii <= range[1]; ++ii)
			// 		{
			// 			t_id = mi->mm.mi_v[ii].y >> 32;
			// 			tstr = (mi->mm.mi_v[ii].y & strand_mask);
			// 			qstr = (seq_mi->mi_v[m_i].y & strand_mask);
			// 			if (p->read_stat->iter_idx[t_id] != query->rid && t_id == tid_tmp && (tstr ^ qstr) == rev_tmp)
			// 			{
			// 				// printf("%d %d %d %d %d %d\n", i, ii, t_id, tid_tmp, tstr^qstr, rev_tmp);
			// 				// printf("%ld %ld %d\n", vertex_mr[sta_idx].qs, vertex_mr[end_idx].qs, p->mi->k);
			// 				pos_q = ((seq_mi->mi_v[m_i].y & pos_mask) >> 1);
			// 				pos_tt = ((mi->mm.mi_v[ii].y & pos_mask) >> 1);
			// 				pos_qs = pos_q + 1 - p->mi->k;
			// 				pos_ts = pos_tt + 1 - p->mi->k;
			// 				if (pos_qs < vertex_mr[sta_idx].qs) // looking forward
			// 				{
			// 					diff_q = (int32_t)(vertex_mr[sta_idx].qs - pos_qs);
			// 					if (rev_tmp == 0) // the same strand
			// 					{
			// 						if (pos_ts >= vertex_mr[sta_idx].ts)	continue;
			// 						diff_tt = (int32_t)(vertex_mr[sta_idx].ts - pos_ts);
			// 						if (abs(diff_q - diff_tt) < Eindel)	add_flag = 1;
			// 					}
			// 					else
			// 					{
			// 						if (pos_ts <= vertex_mr[sta_idx].ts)	continue;
			// 						diff_tt = (int32_t)(pos_tt - vertex_mr[sta_idx].te);
			// 						if (abs(diff_q - diff_tt) < Eindel)	add_flag = 1;

			// 					}
			// 				}
			// 				else if (pos_qs > vertex_mr[end_idx].qs) // looking backward
			// 				{
			// 					diff_q = (int32_t)(pos_qs - vertex_mr[end_idx].qs);
			// 					if (rev_tmp == 0) // the same strand
			// 					{
			// 						if (pos_ts <= vertex_mr[end_idx].ts)	continue;
			// 						diff_tt = (int32_t)(pos_ts - vertex_mr[end_idx].ts);
			// 						if (abs(diff_q - diff_tt) < Eindel)	add_flag = 1;
			// 					}
			// 					else // the oppsite strand
			// 					{
			// 						if (pos_ts >= vertex_mr[end_idx].ts)	continue;
			// 						diff_tt = (int32_t)(vertex_mr[end_idx].te - pos_tt);
			// 						if (abs(diff_q - diff_tt) < Eindel)	add_flag = 1;
			// 					}
			// 				}
			// 				if (add_flag == 0)	continue;
			// 				// printf("%d %d %d %d\n", pos_q, pos_tt, pos_qs, pos_ts);
			// 				// printf("%d %d\n", diff_q, diff_tt);
			// 				if (v_n >= v_m)
			// 				{
			// 					v_m = v_m << 1;
			// 					path_mr = (MR_t *)realloc(path_mr, v_m * sizeof(MR_t));
			// 					if (path_mr == NULL)
			// 					{
			// 						fprintf(stderr, "[%s] calloc %ldGB path_mr memory error!\n", __func__, v_m * sizeof(MR_t) / 1024 / 1024 / 1024);
			// 						exit(1);
			// 					}
			// 				}
			// 				path_mr[v_n].qs = pos_qs;
			// 				path_mr[v_n].qe = pos_q;
			// 				path_mr[v_n].ts = pos_ts;
			// 				path_mr[v_n].te = pos_tt;
			// 				path_mr[v_n].cov = p->mi->k;
			// 				path_mr[v_n].t_id = t_id;
			// 				path_mr[v_n].tstr = tstr;
			// 				path_mr[v_n++].qstr = qstr;
			// 				add_flag = 0;
			// 			}
			// 		}
			// 	}
			// }

			// looking forward
			// if (sta_idx > 0)
			for (j = sta_idx; j >= 0; j--)
			{
				// printf("%d %d\n", j, vertex_mr[j].t_id);
				if (vertex_mr[j].t_id != tid_tmp)	break;
				if ((vertex_mr[j].qstr ^ vertex_mr[j].tstr) != rev_tmp)	continue;
				diff_q = (int32_t)(vertex_mr[sta_idx].qs - vertex_mr[j].qs);
				if (rev_tmp == 0) // the same strand
				{
					if (vertex_mr[j].ts >= vertex_mr[sta_idx].ts)	continue;
					diff_tt = (int32_t)(vertex_mr[sta_idx].ts - vertex_mr[j].ts);
					if (abs(diff_q - diff_tt) < Eindel)
					{
						if (v_n >= v_m)
						{
							v_m = v_m << 1;
							path_mr = (MR_t *)realloc(path_mr, v_m * sizeof(MR_t));
						}
						path_mr[v_n++] = vertex_mr[j];
					}
				}
				else // the oppsite strand
				{
					if (vertex_mr[j].ts <= vertex_mr[sta_idx].ts)	continue;
					diff_tt = (int32_t)(vertex_mr[j].te - vertex_mr[sta_idx].te);
					if (abs(diff_q - diff_tt) < Eindel)
					{
						if (v_n >= v_m)
						{
							v_m = v_m << 1;
							path_mr = (MR_t *)realloc(path_mr, v_m * sizeof(MR_t));
						}
						path_mr[v_n++] = vertex_mr[j];
					}
				}
			}
			// scanning the path
			for (k = 0; k < path_idx[i].idx_n; ++k)
			{
				if (v_n >= v_m)
				{
					v_m = v_m << 1;
					path_mr = (MR_t *)realloc(path_mr, v_m * sizeof(MR_t));
				}
				path_mr[v_n++] = vertex_mr[path_idx[i].idx[k]];
			}
			// looking backward
			// if (end_idx < buf->v_n)
			for (j = end_idx; j < buf->v_n ; j++)
			{
				if (vertex_mr[j].t_id != tid_tmp)	break;
				if ((vertex_mr[j].qstr ^ vertex_mr[j].tstr) != rev_tmp)	continue;
				diff_q = (int32_t)(vertex_mr[j].qs - vertex_mr[end_idx].qs);
				if (rev_tmp == 0) // the same strand
				{
					if (vertex_mr[j].ts <= vertex_mr[end_idx].ts)	continue;
					diff_tt = (int32_t)(vertex_mr[j].ts - vertex_mr[end_idx].ts);
					if (abs(diff_q - diff_tt) < Eindel)
					{
						if (v_n >= v_m)
						{
							v_m = v_m << 1;
							path_mr = (MR_t *)realloc(path_mr, v_m * sizeof(MR_t));
						}
						path_mr[v_n++] = vertex_mr[j];
					}
				}
				else // the oppsite strand
				{
					if (vertex_mr[j].ts >= vertex_mr[end_idx].ts)	continue;
					diff_tt = (int32_t)(vertex_mr[end_idx].te - vertex_mr[j].te);
					if (abs(diff_q - diff_tt) < Eindel)
					{
						if (v_n >= v_m)
						{
							v_m = v_m << 1;
							path_mr = (MR_t *)realloc(path_mr, v_m * sizeof(MR_t));
						}
						path_mr[v_n++] = vertex_mr[j];
					}
				}
			}
			// printf("\n");
			if (v_n > path_idx[i].idx_n)
			{
				qsort(path_mr, v_n, sizeof(MR_t), compare_qs);
				// for (j = 0; j < v_n; ++j)
				// {
				// 	printf("%d\t%d\t%d\t%d\t%ld\t%ld\t%d\t%d\t%ld\t%ld\n", ove_rep[i].qid, j, path_mr[j].cov, path_mr[j].qstr, path_mr[j].qs, path_mr[j].qe, path_mr[j].t_id, path_mr[j].tstr, path_mr[j].ts, path_mr[j].te);
				// }
				// printf("\n");
				merge_colinear_mr(path_mr, &v_n);
				// if (ove_tmp.qid == 82902 || ove_tmp.qid == 453976 || ove_tmp.qid == 604087 || ove_tmp.qid == 203541 || ove_tmp.qid == 1647 || ove_tmp.qid == 2076 || ove_tmp.qid == 3921 || ove_tmp.qid == 15250)
				// {
				// 	for (j = 0; j < v_n; ++j)
				// 	{
				// 		printf("%d\t%d\t%d\t%d\t%ld\t%ld\t%d\t%d\t%ld\t%ld\n", ove_rep[i].qid, j, path_mr[j].cov, path_mr[j].qstr, path_mr[j].qs, path_mr[j].qe, path_mr[j].t_id, path_mr[j].tstr, path_mr[j].ts, path_mr[j].te);
				// 	}
				// 	printf("\n");
				// }
				uint32_t *max_index2, max_index2_n;
				PATH_t *path;
				Graph Dp_graph;

				max_index2 = (uint32_t *)calloc(v_n, sizeof(uint32_t)); // v_n*4 bytes
				path = (PATH_t *)calloc(v_n, sizeof(PATH_t)); // v_n*12 bytes
				init_graph(&Dp_graph, v_n);  // v_n*16 bytes

				max_index2_n = create_graph(data->p->search_step, path_mr, v_n, max_index2, path, &Dp_graph);
				ove_cl_tmp = (OVE_t *)realloc(ove_cl_tmp, (max_index_n+max_index2_n) * sizeof(OVE_t));

				for (k = 0; k < max_index2_n; ++k)
				{
					j = (int32_t)max_index2[k];
					matching_bases = 0;
					ove_tmp.score = path[j].dist;
					// first iter vertex_mr.t_id is the rid of index read; second iter vertex_mr.t_id is the index of index read in array read_stat->iter_idx
					ove_tmp.tid  = path_mr[j].t_id;
					ove_tmp.qid = query->rid;
					ove_tmp.ql = query->read_length;
					ove_tmp.tl = data->p->read_stat->rlen[data->p->read_stat->iter_idx[ove_tmp.tid]];
					ove_tmp.rev = path_mr[j].tstr == path_mr[j].qstr ? 0 : 1;
					ove_tmp.qe = path_mr[j].qe;
					// ove_tmp.te = ove_tmp.rev == 0 ? path_mr[j].te : path_mr[j].ts;
					if (ove_tmp.rev == 0)
						ove_tmp.te = path_mr[j].te;
					else
						ove_tmp.ts = path_mr[j].ts;
					while(j != -1)
					{
						ove_tmp.qs = path_mr[j].qs;
						// ove_tmp.ts = ove_tmp.rev == 0 ? path_mr[j].ts : path_mr[j].te;
						if (ove_tmp.rev == 0)
							ove_tmp.ts = path_mr[j].ts;
						else
							ove_tmp.te = path_mr[j].te;
						matching_bases += path_mr[j].cov;
						j = (int32_t)path[j].pre_node;
						if (j >= (int32_t)buf->v_n) printf("[%s] 1 use unlegal memory, %d-%d, something wrong with it...", __func__, j, buf->v_n);
					}

					qpre = ove_tmp.qs; qsuf = ove_tmp.ql - ove_tmp.qe;
					if (ove_tmp.rev == 0)
						{tpre = ove_tmp.ts; tsuf = ove_tmp.tl - ove_tmp.te;}
					else if (ove_tmp.rev == 1)
						{tpre = ove_tmp.tl - ove_tmp.te; tsuf = ove_tmp.ts;}

					left_overhang  = qpre < tpre ? qpre : tpre;
					right_overhang = qsuf < tsuf ? qsuf : tsuf;
					overhang = left_overhang + right_overhang;
					ove_tmp.mbp = matching_bases;
					ove_tmp.mln = ove_tmp.te - ove_tmp.ts + 1;

					overhang_tr = data->p->max_hang < ove_tmp.mln * 0.8 ? data->p->max_hang : ove_tmp.mln * 0.8;
					mln_th = ove_tmp.ql < ove_tmp.tl ? ove_tmp.ql : ove_tmp.tl;
					mln_ratio = data->p->min_ove > mln_th * 0.3 ? data->p->min_ove : mln_th * 0.3;
					// mln_ratio = mln_ratio > 5000 ? 5000 : mln_ratio;
					if (matching_bases >= data->p->m_b && overhang < overhang_tr && ove_tmp.mln > mln_ratio)
					{
						ove_cl_tmp[n_tmp++] = ove_tmp;
						if (n_tmp > (max_index_n+max_index2_n)) printf("[%s]3 use unlegal memory, %d-%d, something wrong with it...", __func__, n_tmp, (max_index_n+max_index2_n));
						// print_single_read_overlap(ove_tmp);
					}
					// if (ove_tmp.qid == 82902 || ove_tmp.qid == 453976 || ove_tmp.qid == 604087 || ove_tmp.qid == 203541 || ove_tmp.qid == 1647 || ove_tmp.qid == 2076 || ove_tmp.qid == 3921 || ove_tmp.qid == 15250)
						// print_single_read_overlap(ove_tmp);
				}

				if (max_index2 != NULL) {free(max_index2);max_index2 = NULL;}
				if (path != NULL) {free(path);path = NULL;}
				free_graph(&Dp_graph);
			}

			if (path_mr != NULL) {free(path_mr); path_mr = NULL;}
		}
		// printf("\n");
	}

	if (n_tmp == 0)
	{
		for (i = 0; i < n_rep; ++i)	if (path_idx[i].idx != NULL) {free(path_idx[i].idx); path_idx[i].idx = NULL;}
		if (path_idx != NULL) {free(path_idx); path_idx = NULL;}
		if (tmp_path_idx != NULL) {free(tmp_path_idx); tmp_path_idx = NULL;}
		if (ove_rep != NULL) {free(ove_rep); ove_rep = NULL;}
		if (ove_cl_tmp != NULL) {free(ove_cl_tmp); ove_cl_tmp = NULL;}
		return 0;
	}

	qsort(ove_cl_tmp, n_tmp, sizeof(OVE_t), compare_ove_score);
	n_tmp = n_tmp > data->p->top_n ? data->p->top_n : n_tmp;

	for (i = 0; i < n_tmp; i++)
	{
		if (buf->ove_cl[q_idx].n >= buf->ove_cl[q_idx].m)
		{
			buf->ove_cl[q_idx].m = buf->ove_cl[q_idx].m << 1;
			buf->ove_cl[q_idx].ove = (OVE_t *)realloc(buf->ove_cl[q_idx].ove, buf->ove_cl[q_idx].m * sizeof(OVE_t));
			if (buf->ove_cl[q_idx].ove == NULL)
			{
				fprintf(stderr, "[%s] calloc %ldGB buf->ove_cl[q_idx].ove memory error!\n", __func__, buf->ove_cl[q_idx].m * sizeof(OVE_t) / 1024 / 1024 / 1024);
				exit(1);
			}
			printf("[%s] will not use these code, something wrong with it...",__func__);
		}
		// cp_ove(&buf->ove_cl[q_idx].ove[buf->ove_cl[q_idx].n], &ove_cl_tmp[i]);
		buf->ove_cl[q_idx].ove[buf->ove_cl[q_idx].n] = ove_cl_tmp[i];
		buf->ove_cl[q_idx].ove[buf->ove_cl[q_idx].n++].tid = data->p->read_stat->iter_idx[ove_cl_tmp[i].tid];
		if (buf->ove_cl[q_idx].n > buf->ove_cl[q_idx].m) printf("[%s]buf[%d]->ove_cl[%d] memory leak..,%d > %d\n",__func__,pid,q_idx,buf->ove_cl[q_idx].n,buf->ove_cl[q_idx].m);
	}

	for (i = 0; i < n_rep; ++i)	if (path_idx[i].idx != NULL) {free(path_idx[i].idx); path_idx[i].idx = NULL;}
	if (path_idx != NULL) {free(path_idx); path_idx = NULL;}
	if (tmp_path_idx != NULL) {free(tmp_path_idx); tmp_path_idx = NULL;}
	if (ove_rep != NULL) {free(ove_rep); ove_rep = NULL;}
	if (ove_cl_tmp != NULL) {free(ove_cl_tmp); ove_cl_tmp = NULL;}

	return n_tmp;
}

// the i is the read array's index
static void overlapping_core(void *data, int64_t i, int pid)
{
	step_query_t *step = (step_query_t *)data;
	mini_t *seq_mi;
	uint32_t *max_index, max_index_n;
	PATH_t *path;
	Graph Dp_graph;
	// READ_t *read_info = &step->read_info[i];
	seq_mi = (mini_t *)calloc(1, sizeof(mini_t));
	seq_mi->mi_n = 0;
	seq_mi->mi_m = step->read_info[i].read_length;
	seq_mi->mi_v = (MINIMIZER_t *)calloc(seq_mi->mi_m, sizeof(MINIMIZER_t)); // max_len*16 bytes

	finding_hits(step, i, pid, seq_mi);
	if (step->buf[pid]->v_n > 0)
	{
		max_index = (uint32_t *)calloc(step->buf[pid]->v_n, sizeof(uint32_t)); // v_n*4 bytes
		path = (PATH_t *)calloc(step->buf[pid]->v_n, sizeof(PATH_t)); // v_n*12 bytes
		init_graph(&Dp_graph, step->buf[pid]->v_n);  // v_n*16 bytes

		max_index_n = create_graph(step->p->search_step, step->buf[pid]->vertex_mr, step->buf[pid]->v_n, max_index, path, &Dp_graph);
		get_overlap_info(step, step->buf[pid], i, pid, max_index_n, max_index, path, seq_mi);

		if (max_index != NULL) {free(max_index);max_index = NULL;}
		if (path != NULL) {free(path);path = NULL;}
		free_graph(&Dp_graph);
	}

	if (seq_mi->mi_v != NULL) {free(seq_mi->mi_v);seq_mi->mi_v = NULL;}
	if (seq_mi != NULL) {free(seq_mi);seq_mi = NULL;}

	return;
}

static void *mapping_pipeline(void *shared, int step, void *in)
{
	int i, j, k;
	pipeline_t *p = (pipeline_t *)shared;
	if (step == 0)
	{
		step_query_t *s = (step_query_t *)calloc(1, sizeof(step_query_t));
		// fprintf(stderr, "[%s step0: %.3fs, %.3fGB] start loading query read...\n", __func__, realtime() - realtime0, peak_memory() / 1024.0 / 1024.0);
		s->read_info = bseq_read_map(p->fp, p->query_batch_size, &s->n_seq, &p->map_n_sta, &p->map_idx_n_sta, p->read_stat->iter_map, &p->file_idx);
		if (s->n_seq)
		{
			s->p = p;
			s->ove_alloc = s->n_seq;
			p->n_processed += s->n_seq;
			p->map_batch++;

			s->buf = (thread_buf_t **)calloc(p->n_threads, sizeof(thread_buf_t *));
			for (i = 0; i < p->n_threads; i++)
			{
				s->buf[i] = (thread_buf_t *)calloc(1, sizeof(thread_buf_t));
				// s->buf[i]->seq_mi = (mini_t *)calloc(1, sizeof(mini_t));
				// s->buf[i]->seq_mi->mi_n = 0;
				// s->buf[i]->seq_mi->mi_m = p->max_len;
				// s->buf[i]->seq_mi->mi_v = (MINIMIZER_t *)calloc(s->buf[i]->seq_mi->mi_m, sizeof(MINIMIZER_t)); // max_len*16 bytes
				s->buf[i]->max_v_num = p->rep_n * p->max_len / p->mi->mi_step;
				s->buf[i]->v_n = 0;
				s->buf[i]->v_m = s->buf[i]->max_v_num; // max_v_num = rep_n*max_len/mi_step = [5-100]*[50k-100k]/[1-3] = [80k,10M]
				s->buf[i]->vertex_mr = (MR_t *)calloc(s->buf[i]->v_m, sizeof(MR_t)); // max_v_num*48 bytes
				// init_graph(&s->buf[i]->Dp_graph, s->buf[i]->max_v_num);  // max_v_num*16 bytes
				// s->buf[i]->path = (PATH_t *)calloc(s->buf[i]->max_v_num, sizeof(PATH_t)); // max_v_num*12 bytes
				// s->buf[i]->max_index = (uint32_t *)calloc(s->buf[i]->max_v_num, sizeof(uint32_t)); // max_v_num*4 bytes
				s->buf[i]->ove_cl = (OVE_C_t *)calloc(s->ove_alloc, sizeof(OVE_C_t)); // #reads*(16+top_n*48)bytes
				for (j = 0; j < s->ove_alloc; j++)
				{
					s->buf[i]->ove_cl[j].n = 0;
					s->buf[i]->ove_cl[j].m = p->top_n;
					s->buf[i]->ove_cl[j].ove = (OVE_t *)calloc(s->buf[i]->ove_cl[j].m, sizeof(OVE_t));
				}
				// graph.c: vertex_num*(search_step*16 + 4 + 1) bytes
				s->buf[i]->is_ove = (uint32_t *)calloc(s->ove_alloc, sizeof(uint32_t));
			} 
			fprintf(stderr, "[Mapping: %.3fs, %.3fGB] Loaded %d query read...\n", realtime() - realtime0, peak_memory() / 1024.0 / 1024.0, s->n_seq);
			return s;
		}
		else
		{
			if (s != NULL) {free(s); s = NULL;}
		}
	}
	else if (step == 1)
	{
		kt_for(p->n_threads, overlapping_core, in, ((step_query_t *)in)->n_seq);
		fprintf(stderr, "[Mapping: %.3fs, %.3fGB] Mapped %d query read...\n", realtime() - realtime0, peak_memory() / 1024.0 / 1024.0, ((step_query_t *)in)->n_seq);
		return in;
	}
	else if (step == 2)
	{
		uint32_t id;
		step_query_t *s = (step_query_t *)in;

		for (i = 0; i < p->n_threads; i++)
		{
			for (j = 0; j < s->ove_alloc; j++)
			{
				id = s->read_info[j].rid;
				if (s->buf[i]->is_ove[j] == 1) {p->read_stat->is_ove[id] = 1;}
				if (s->buf[i]->ove_cl[j].n == 0) continue;
				id = s->buf[i]->ove_cl[j].ove[0].qid;

				for (k = 0; k < s->buf[i]->ove_cl[j].n; k++)
				{
					if (p->ove_cl[id].n >= p->ove_cl[id].m)
					{
						p->ove_cl[id].m = p->ove_cl[id].m << 1;
						p->ove_cl[id].ove = (OVE_t *)realloc(p->ove_cl[id].ove, p->ove_cl[id].m * sizeof(OVE_t));
					}
					// cp_ove(&p->ove_cl[id].ove[p->ove_cl[id].n++], &s->buf[i]->ove_cl[j].ove[k]);
					p->ove_cl[id].ove[p->ove_cl[id].n++] = s->buf[i]->ove_cl[j].ove[k];
					if (p->ove_cl[id].n > p->ove_cl[id].m) printf("[%s]p->ove_cl[%d] memory leak..,%d > %d\n",__func__,id,p->ove_cl[id].n,p->ove_cl[id].m);
				}
				if (p->ove_cl[id].n - p->symm[id] > p->top_n)
				{
					qsort(&p->ove_cl[id].ove[p->symm[id]], p->ove_cl[id].n - p->symm[id], sizeof(OVE_t), compare_ove_score);
					p->ove_cl[id].n = p->symm[id] + p->top_n;
				}
			}
		}

		for (i = 0; i < s->n_seq; i++)
		{
			READ_t *read = &s->read_info[i];
			if (read->read_seq != NULL) {free(read->read_seq); read->read_seq = NULL;}
			if (read->read_name != NULL) {free(read->read_name); read->read_name = NULL;}
		}
		if (s->read_info != NULL) {free(s->read_info); s->read_info = NULL;}

		for (i = 0; i < p->n_threads; i++)
		{
			// if (s->buf[i]->seq_mi->mi_v != NULL) {free(s->buf[i]->seq_mi->mi_v);s->buf[i]->seq_mi->mi_v = NULL;}
			// if (s->buf[i]->seq_mi != NULL) {free(s->buf[i]->seq_mi);s->buf[i]->seq_mi = NULL;}
			if (s->buf[i]->vertex_mr != NULL) {free(s->buf[i]->vertex_mr);s->buf[i]->vertex_mr = NULL;}
			// if (s->buf[i]->path != NULL) {free(s->buf[i]->path);s->buf[i]->path = NULL;}
			// if (s->buf[i]->max_index != NULL) {free(s->buf[i]->max_index);s->buf[i]->max_index = NULL;}
			// free_graph(&s->buf[i]->Dp_graph);
			for (j = 0; j < s->ove_alloc; j++)
			{
				if (s->buf[i]->ove_cl[j].ove != NULL) {free(s->buf[i]->ove_cl[j].ove);s->buf[i]->ove_cl[j].ove = NULL;}
			}
			if (s->buf[i]->ove_cl != NULL) {free(s->buf[i]->ove_cl);s->buf[i]->ove_cl = NULL;}
			if (s->buf[i]->is_ove != NULL) {free(s->buf[i]->is_ove); s->buf[i]->is_ove = NULL;}
			if (s->buf[i] != NULL) {free(s->buf[i]);s->buf[i] = NULL;}
		}
		if (s->buf != NULL) {free(s->buf); s->buf = NULL;}
		if (s != NULL) {free(s); s = NULL;}
	}
	return 0;
}


static void estimate_coverage_of_seed_blocks(void *data, int64_t i, int pid)
// static void estimate_coverage_of_seed_blocks(void *data, int64_t i, uint8_t *visited_rid)
{
	pipeline_t *p = (pipeline_t *)data;
	uint32_t j, k, oi;
	uint32_t seed_rid, read_rid, trans_rid;
	int32_t qs, qe, ql, ts, te, tl;
	int32_t t_qs, t_qe, t_ql, t_ts, t_te, t_tl;
	uint8_t rev, t_rev;
	int32_t q_sta_pos, q_end_pos;
	int32_t t_q_sta_pos, t_q_end_pos;
	uint32_t cov_q_sta, cov_q_end;
	uint32_t size_tr = p->split_len / 2, sum_cov;

	if (p->read_stat->is_idx[i] == 1)
	{
		seed_rid = i;
		// printf("seed read %d\n", seed_rid);
		for(j = 0; j < p->symm[seed_rid]; ++j)
		{
			read_rid = p->ove_cl[seed_rid].ove[j].tid;
			qs = p->ove_cl[seed_rid].ove[j].qs;
			qe = p->ove_cl[seed_rid].ove[j].qe;
			ql = p->ove_cl[seed_rid].ove[j].ql;
			rev = p->ove_cl[seed_rid].ove[j].rev;
			ts = p->ove_cl[seed_rid].ove[j].ts;
			te = p->ove_cl[seed_rid].ove[j].te;
			tl = p->ove_cl[seed_rid].ove[j].tl;
			if (p->visited_rid[pid][read_rid] != i)
			{
				if (rev == 0)
				{
					q_sta_pos = qs > ts ? qs - ts : 0;
					q_end_pos = qe + (tl - te) < ql ? qe + (tl - te) : ql;
				}
				else
				{
					q_sta_pos = qs > (tl - te) ? qs - (tl - te) : 0;
					q_end_pos = qe + ts < ql ? qe + ts : ql;
				}

				if (q_sta_pos / p->split_len == p->blkn[seed_rid] - 1)
					cov_q_sta = p->blkn[seed_rid] - 1;
				else
					cov_q_sta = q_sta_pos % p->split_len < size_tr ? q_sta_pos / p->split_len : q_sta_pos / p->split_len + 1;
				if (q_end_pos / p->split_len == p->blkn[seed_rid] - 1) // the last region of seed read < p->split_len bp
					cov_q_end = p->blkn[seed_rid] - 1;
				else if (q_end_pos / p->split_len == 0)
					cov_q_end = 0;
				else
					cov_q_end = q_end_pos % p->split_len > size_tr ? q_end_pos / p->split_len : q_end_pos / p->split_len - 1;
				for(k = cov_q_sta; k <= cov_q_end; ++k)
				{
					p->read_cov[seed_rid].cov[k] += 1;
				}
				p->visited_rid[pid][read_rid] = i;
			}

			if (p->read_stat->is_idx[read_rid] == 1) // processing seed reads, transitive one round
			{
				// printf("seed read2 %d, %d\n", read_rid, ql);
				// printf("%d %d\n", q_sta_pos, q_end_pos);
				// printf("%d %d\n", cov_q_sta, cov_q_end);
				for(oi = 0; oi < p->symm[read_rid]; ++oi)
				{
					trans_rid = p->ove_cl[read_rid].ove[oi].tid;
					if (p->visited_rid[pid][trans_rid] == i)	continue;

					t_qs = p->ove_cl[read_rid].ove[oi].qs;
					t_qe = p->ove_cl[read_rid].ove[oi].qe;
					t_ql = p->ove_cl[read_rid].ove[oi].ql;
					t_rev = p->ove_cl[read_rid].ove[oi].rev;
					t_ts = p->ove_cl[read_rid].ove[oi].ts;
					t_te = p->ove_cl[read_rid].ove[oi].te;
					t_tl = p->ove_cl[read_rid].ove[oi].tl;
					if (t_rev == 0)
					{
						t_q_sta_pos = t_qs > t_ts ? t_qs - t_ts : 0;
						t_q_end_pos = t_qe + (t_tl - t_te) < t_ql ? t_qe + (t_tl - t_te) : t_ql;
					}
					else
					{
						t_q_sta_pos = t_qs > (t_tl - t_te) ? t_qs - (t_tl - t_te) : 0;
						t_q_end_pos = t_qe + t_ts < t_ql ? t_qe + t_ts : t_ql;
					}
					if (rev == 0)
					{
						q_sta_pos = qs + t_q_sta_pos - ts;
						q_end_pos = qe + t_q_end_pos - te;
					}
					else
					{
						q_sta_pos = qs + te - t_q_end_pos;
						q_end_pos = qe + ts - t_q_sta_pos;
					}
					if (q_sta_pos > ql || q_end_pos < 0)	continue;
					// printf("1 %d %d\n", t_q_sta_pos, t_q_end_pos);
					// printf("2 %d %d\n", q_sta_pos, q_end_pos);
					q_sta_pos = q_sta_pos < 0 ? 0 : q_sta_pos;
					q_end_pos = q_end_pos > ql ? ql : q_end_pos;
					// printf("3 %d %d\n", q_sta_pos, q_end_pos);

					if (q_sta_pos / p->split_len == p->blkn[seed_rid] - 1)
						cov_q_sta = p->blkn[seed_rid] - 1;
					else
						cov_q_sta = q_sta_pos % p->split_len < size_tr ? q_sta_pos / p->split_len : q_sta_pos / p->split_len + 1;
					if (q_end_pos / p->split_len == p->blkn[seed_rid] - 1) // the last region of seed read < p->split_len bp
						cov_q_end = p->blkn[seed_rid] - 1;
					else if (q_end_pos / p->split_len == 0)
						cov_q_end = 0;
					else
						cov_q_end = q_end_pos % p->split_len > size_tr ? q_end_pos / p->split_len : q_end_pos / p->split_len - 1;

					for(k = cov_q_sta; k <= cov_q_end; ++k)
					{
						p->read_cov[seed_rid].cov[k] += 1;
					}
					p->visited_rid[pid][trans_rid] = i;
				}
			}
		}
		sum_cov = 0;
		// printf("seed read %d\t", seed_rid);
		for(k = 0; k < p->blkn[seed_rid]; ++k)
		{
			sum_cov += p->read_cov[seed_rid].cov[k];
			// printf("%d\t",p->read_cov[seed_rid].cov[k]);
		}
		p->ave_cov[seed_rid] = sum_cov / (double)p->blkn[seed_rid];
		// printf("ave %.3f\n", p->ave_cov[seed_rid]);
	}
}

static void estimate_coverage_of_read_blocks(void *data, int64_t i, int pid)
{
	pipeline_t *p = (pipeline_t *)data;
	uint32_t seed_rid, sum_cov;
	// uint32_t last_size;
	double ave;
	int32_t cov_t_sta, cov_t_end, cov_t_sta_pos, cov_t_end_pos;
	int32_t k, cov_q_sta, cov_q_end, cov_q_sta_pos, cov_q_end_pos;
	uint32_t j, cov_nn = p->max_len / p->split_len + 1;
	uint32_t *cov_n = (uint32_t *)calloc(cov_nn, sizeof(uint32_t));
	uint32_t *tmp_cov = (uint32_t *)calloc(cov_nn, sizeof(uint32_t));
	uint8_t rev;
	uint32_t ts, te, qs, qe;
	int32_t ql, tl;
	uint32_t size_tr = p->split_len / 2;
	if (p->read_stat->is_idx[i] != 1)
	{
		// step 2: estimate coverage of other read
		for (j = 0; j < cov_nn; ++j) cov_n[j] = 0;
		for (j = 0; j < p->symm[i]; ++j)
		{
			seed_rid = p->ove_cl[i].ove[j].tid;
			rev = p->ove_cl[i].ove[j].rev;
			ts = p->ove_cl[i].ove[j].ts;
			te = p->ove_cl[i].ove[j].te;
			qs = p->ove_cl[i].ove[j].qs;
			qe = p->ove_cl[i].ove[j].qe;
			tl = (int32_t)p->read_stat->rlen[seed_rid];
			ql = (int32_t)p->read_stat->rlen[i];
			
			if (rev == 0)
			{
				cov_t_sta_pos = ts > qs ? ts - qs : 0;
				cov_t_end_pos = te + (ql - qe) < tl ? te + (ql - qe) : tl;
				cov_q_sta_pos = ts < qs ? qs - ts : 0;
				cov_q_end_pos = qe + (tl - te) < ql ? qe + (tl - te) : ql;
			}
			else
			{
				cov_t_sta_pos = ts > (ql - qe) ? ts - (ql - qe) : 0;
				cov_t_end_pos = te + qs < tl ? te + qs : tl;
				cov_q_sta_pos = qs > (tl - te) ? qs - (tl - te) : 0;
				cov_q_end_pos = qe + ts < ql ? qe + ts : ql;
			}

			// cov_t_sta_pos = rev == 0 ? ts - qs : te - qs;
			// cov_t_end_pos = rev == 0 ? te + (ql - qe) : ts + (ql - qe);
			// cov_t_sta_pos = cov_t_sta_pos < 0 ? 0 : cov_t_sta_pos;
			// cov_t_end_pos = cov_t_end_pos > tl ? tl : cov_t_end_pos;
			cov_t_sta = cov_t_sta_pos / p->split_len;
			cov_t_end = cov_t_end_pos / p->split_len;
			sum_cov = 0;
			for (k = cov_t_sta; k <= cov_t_end; ++k)
			{
				sum_cov += p->read_cov[seed_rid].cov[k];
			}
			ave = sum_cov / (double)(cov_t_end - cov_t_sta + 1);

			// cov_q_sta_pos = rev == 0 ? qs - ts : qs - (p->read_stat->rlen[seed_rid] - ts);
			// cov_q_end_pos = rev == 0 ? qe + (p->read_stat->rlen[seed_rid] - te) : qe + te;
			// cov_q_sta_pos = cov_q_sta_pos < 0 ? 0 : cov_q_sta_pos;
			// cov_q_end_pos = cov_q_end_pos > ql ? ql : cov_q_end_pos;
			// cov_q_sta = cov_q_sta_pos / p->split_len;
			// cov_q_end = cov_q_end_pos / p->split_len;
			if (cov_q_sta_pos / p->split_len == p->blkn[i] - 1)
				cov_q_sta = p->blkn[i] - 1;
			else
				cov_q_sta = cov_q_sta_pos % p->split_len < size_tr ? cov_q_sta_pos / p->split_len : cov_q_sta_pos / p->split_len + 1;
			if (cov_q_end_pos / p->split_len == p->blkn[i] - 1) // the last region of seed read < p->split_len bp
				cov_q_end = p->blkn[i] - 1;
			else if (cov_q_end_pos / p->split_len == 0)
				cov_q_end = 0;
			else
				cov_q_end = cov_q_end_pos % p->split_len > size_tr ? cov_q_end_pos / p->split_len : cov_q_end_pos / p->split_len - 1;
			for (k = cov_q_sta; k <= cov_q_end; ++k)
			{
				tmp_cov[k] += ave;
				cov_n[k] += 1;
			}
		}
		sum_cov = 0;
		// printf("map read %ld\t", i);
		for (k = 0; k < p->blkn[i]; ++k)
		{
			if (cov_n[k] == 0)
			{
				// printf("%d\t", p->read_cov[i].cov[k]);
				continue;
			}
			p->read_cov[i].cov[k] = tmp_cov[k] / cov_n[k];
			sum_cov += p->read_cov[i].cov[k];
			// printf("%d\t", p->read_cov[i].cov[k]);
		}
		p->ave_cov[i] = sum_cov / (double)p->blkn[i];
		// printf("ave %.3f\n", p->ave_cov[i]);
	}
	if (cov_n != NULL) {free(cov_n); cov_n = NULL;}
	if (tmp_cov != NULL) {free(tmp_cov); tmp_cov = NULL;}

}

static int choose_idx_read_for_next_iteration(pipeline_t *pl, read_cov_t *read_cov, double ave_median)
{
	int ret, r_i, k, tmp;
	uint32_t buf_size = part_size;
	uint32_t *buf = (uint32_t *)calloc(buf_size, sizeof(uint32_t));
	uint32_t buf_pos, buf_cov, ele;
	double buf_ave_cov;
	uint32_t last_map_n = pl->read_stat->iter_map_n;

	pl->read_stat->iter_idx_n = 0;
	pl->read_stat->iter_map_n = 0;

	double cov_tr = ave_median * cov_ratio;
	for (r_i = 0; r_i < pl->n_total_read; r_i++)
	{
		srand(r_i * (pl->iter + 1));
		// printf("new read %d, %d block\n", r_i, pl->blkn[r_i]);
		ret = 0;
		buf_pos = 0;
		buf_cov = 0;
		for (k = 0; k < pl->blkn[r_i]; k++)
		{
			ele = buf[buf_pos];
			buf[buf_pos] = read_cov[r_i].cov[k];
			buf_cov += buf[buf_pos++];
			if (k >= buf_size)
			{
				buf_cov -= ele;
				buf_ave_cov = (double)buf_cov / buf_size;
				// printf("%d %f\n", k, buf_ave_cov);
				if (buf_ave_cov <= cov_tr)	{ret = 1; break;}
			}
			if (buf_pos == buf_size) buf_pos = 0;
		}
		if (pl->blkn[r_i] <= buf_size)
		{
			buf_ave_cov = buf_cov / pl->blkn[r_i];
			// printf("%d %f\n", pl->blkn[r_i], buf_ave_cov);
			if (buf_ave_cov <= cov_tr)	ret = 1;
		}
		if (ret == 1)
		{
			pl->read_stat->iter_map[pl->read_stat->iter_map_n++] = r_i;
			tmp = rand() % 100;
			if (tmp < X_read)
			{
				pl->read_stat->iter_idx[pl->read_stat->iter_idx_n++] = r_i;
				pl->read_stat->is_idx[r_i] = 1;
			}
		}
	}
	fprintf(stderr, "[To next iteration: %.3fs, %.3fGB] %d reads do not have enough coverage( < %f), median coverage %.6f, index %d read...\n", realtime() - realtime0, peak_memory() / 1024.0 / 1024.0, pl->read_stat->iter_map_n, cov_tr, ave_median, pl->read_stat->iter_idx_n);
	if (pl->read_stat->iter_map_n < last_map_n * 0.9) ret = 1;
	else ret = 0;

	if (buf != NULL) {free(buf); buf = NULL;}
	return ret;
}

static int is_redundant_ove(OVE_C_t *seed_ove_cl, uint32_t r_i)
{
	uint32_t i;
	for(i = 0; i < seed_ove_cl->n; ++i)
	{
		if (seed_ove_cl->ove[i].tid == r_i)	return 1;
	}
	return 0;
}

static int construct_symmetrical_graph(uint32_t n_total_read, uint32_t *symm, OVE_C_t *ove_cl, read_stat_t *read_stat)
{
	uint32_t r_i, o_i;
	for (r_i = 0; r_i < n_total_read; ++r_i) symm[r_i] = ove_cl[r_i].n;
	for (r_i = 0; r_i < n_total_read; ++r_i)
	{
		if (ove_cl[r_i].n == 0) continue;
		for (o_i = 0; o_i < ove_cl[r_i].n; ++o_i)
		{
			int seed_rid = ove_cl[r_i].ove[o_i].tid;
			if (symm[seed_rid] >= ove_cl[seed_rid].m)
			{
				ove_cl[seed_rid].m = ove_cl[seed_rid].m << 1;
				ove_cl[seed_rid].ove = (OVE_t *)realloc(ove_cl[seed_rid].ove, ove_cl[seed_rid].m * sizeof(OVE_t));
			}

			if (read_stat->is_idx[r_i] == 1 && is_redundant_ove(&ove_cl[seed_rid], r_i) == 1) continue;
			ove_cl[seed_rid].ove[symm[seed_rid]].qid = ove_cl[r_i].ove[o_i].tid;
			ove_cl[seed_rid].ove[symm[seed_rid]].ql = ove_cl[r_i].ove[o_i].tl;
			ove_cl[seed_rid].ove[symm[seed_rid]].qs = ove_cl[r_i].ove[o_i].ts;
			ove_cl[seed_rid].ove[symm[seed_rid]].qe = ove_cl[r_i].ove[o_i].te;
			ove_cl[seed_rid].ove[symm[seed_rid]].tid = ove_cl[r_i].ove[o_i].qid;
			ove_cl[seed_rid].ove[symm[seed_rid]].tl = ove_cl[r_i].ove[o_i].ql;
			ove_cl[seed_rid].ove[symm[seed_rid]].ts = ove_cl[r_i].ove[o_i].qs;
			ove_cl[seed_rid].ove[symm[seed_rid]].te = ove_cl[r_i].ove[o_i].qe;
			ove_cl[seed_rid].ove[symm[seed_rid]].rev = ove_cl[r_i].ove[o_i].rev;
			ove_cl[seed_rid].ove[symm[seed_rid]].mbp = ove_cl[r_i].ove[o_i].mbp;
			ove_cl[seed_rid].ove[symm[seed_rid]].mln = ove_cl[r_i].ove[o_i].mln;
			ove_cl[seed_rid].ove[symm[seed_rid]++].score = ove_cl[r_i].ove[o_i].score;
		}
	}
	return 0;
}

static int generate_transitive_edges_of_graph(void *data, int64_t i, FILE *fp_trans_ove, uint32_t *visited_rid)
{
	pipeline_t *p = (pipeline_t *)data;
	uint32_t read_rid, j;
	int32_t q_sta_pos, q_end_pos;
	cov_data_t last_read, cur_read;

	int32_t q_n = 0, q_m = trans_iter * 3, q_cur = 0;
	cov_data_t *queue = (cov_data_t *)calloc(q_m, sizeof(cov_data_t));

	int cur_iter = 0, cur_iter_read_num = 0;
	int *trans_iter_read_num = (int *)calloc(trans_iter + 2, sizeof(int));

	read_rid = i;
	cur_read.rid = read_rid;
	cur_read.sta_pos = cur_read.end_pos = 0;
	cur_read.ts = cur_read.te = 0;
	cur_read.rev = 0;
	queue[q_n++] = cur_read;
	trans_iter_read_num[cur_iter] = 1;
	visited_rid[cur_read.rid] = i;
	while ((q_n-q_cur) > 0)
	{
		last_read = queue[q_cur++];
		cur_iter_read_num += 1;
		// printf("new round %d\n", last_read.rid);
		// printf("%d\t%d\n", cur_iter, cur_iter_read_num);
		// for(j = 0; j < trans_iter+2; ++ j) printf("%d\t",trans_iter_read_num[j]);
		// printf("\n");
		for(j = 0; j < p->symm[last_read.rid]; ++ j)
		{
			if (visited_rid[p->ove_cl[last_read.rid].ove[j].tid] == i)	continue;
			visited_rid[p->ove_cl[last_read.rid].ove[j].tid] = i;
			if (last_read.rev == 1)
			{
				q_sta_pos = last_read.sta_pos - (p->ove_cl[last_read.rid].ove[j].qe - last_read.te);
				q_end_pos = last_read.end_pos - (p->ove_cl[last_read.rid].ove[j].qs - last_read.ts);
			}
			else
			{
				q_sta_pos = last_read.sta_pos + p->ove_cl[last_read.rid].ove[j].qs - last_read.ts;
				q_end_pos = last_read.end_pos + p->ove_cl[last_read.rid].ove[j].qe - last_read.te;
			}
			if (q_sta_pos > (int32_t)p->read_stat->rlen[read_rid] || q_end_pos < 0)	continue;

			if (p->read_stat->is_idx[p->ove_cl[last_read.rid].ove[j].tid] == 1)
			{
				cur_read.rid = p->ove_cl[last_read.rid].ove[j].tid;
				cur_read.sta_pos = q_sta_pos;
				cur_read.end_pos = q_end_pos;
				cur_read.ts = p->ove_cl[last_read.rid].ove[j].ts;
				cur_read.te = p->ove_cl[last_read.rid].ove[j].te;
				cur_read.rev = last_read.rev ^ p->ove_cl[last_read.rid].ove[j].rev;
				if (q_n >= q_m)
				{
					q_m = q_m << 1;
					queue = (cov_data_t *)realloc(queue, q_m * sizeof(cov_data_t));
				}
				queue[q_n++] = cur_read;
				// printf("push %d\n", p->ove_cl[last_read.rid].ove[j].tid);
				trans_iter_read_num[cur_iter+1] += 1;
			}

			OVE_t trans_ove;
			trans_ove.qid = read_rid;
			trans_ove.ql = p->read_stat->rlen[read_rid];
			trans_ove.tid = p->ove_cl[last_read.rid].ove[j].tid;
			trans_ove.tl = p->read_stat->rlen[p->ove_cl[last_read.rid].ove[j].tid];
			trans_ove.rev = last_read.rev ^ p->ove_cl[last_read.rid].ove[j].rev;
			trans_ove.qs = q_sta_pos;
			trans_ove.qe = q_end_pos;
			trans_ove.ts = p->ove_cl[last_read.rid].ove[j].ts;
			trans_ove.te = p->ove_cl[last_read.rid].ove[j].te;
			if (trans_ove.rev == 0)
			{
				if (q_sta_pos < 0)
				{
					trans_ove.qs = 0;
					trans_ove.ts = trans_ove.ts - q_sta_pos;
				}
				if (q_end_pos > trans_ove.ql)
				{
					trans_ove.qe = trans_ove.ql;
					trans_ove.te = trans_ove.te - (q_end_pos - trans_ove.ql);
				}
			}
			else
			{
				if (q_sta_pos < 0)
				{
					trans_ove.qs = 0;
					trans_ove.te = trans_ove.te + q_sta_pos;
				}
				if (q_end_pos > trans_ove.ql)
				{
					trans_ove.qe = trans_ove.ql;
					trans_ove.ts = trans_ove.ts + (q_end_pos - trans_ove.ql);
				}
			}
			trans_ove.mbp = 0;
			trans_ove.mln = ((trans_ove.qe - trans_ove.qs) + (trans_ove.te - trans_ove.ts)) / 2;
			trans_ove.score = 0;

			// fprintf(fp_trans_ove, "%d,%d\n", trans_ove.qid, trans_ove.tid);
			fprintf(fp_trans_ove, "%d\t%d\t%d\t%d\t%c\t", trans_ove.qid, trans_ove.ql, trans_ove.qs, trans_ove.qe, "+-"[trans_ove.rev]);
			fprintf(fp_trans_ove, "%d\t%d\t%d\t%d\t", trans_ove.tid, trans_ove.tl, trans_ove.ts, trans_ove.te);
			fprintf(fp_trans_ove, "%d\t%d\t%.2f\n", trans_ove.mbp, trans_ove.mln, trans_ove.score);
		}
		if (cur_iter_read_num == trans_iter_read_num[cur_iter])
		{
			cur_iter++;
			if (cur_iter > trans_iter)	break;
			cur_iter_read_num = 0;
		}
	}

	if (trans_iter_read_num != NULL) {free(trans_iter_read_num); trans_iter_read_num = NULL;}
	if (queue != NULL) {free(queue); queue = NULL;}
	return 0;
}

#ifdef LOAD_INTER_RESULT
uint32_t load_inter_result_info(const char *temp_batch_dir, pipeline_t *pl, uint32_t *batch_iter)
{
	uint32_t i;
	FILE *fp = NULL;
	fp = fopen(temp_batch_dir, "r");

	if(fscanf(fp, "%d %d %d %d\n", batch_iter,&pl->index_n_sta,&pl->index_idx_n_sta,&pl->index_idx_n_end) == 4)
	{
		fprintf(stderr,"%d %d %d %d\n",(*batch_iter),pl->index_n_sta,pl->index_idx_n_sta,pl->index_idx_n_end);
	}
	else exit(1);

	if (fscanf(fp, "%d\n", &pl->read_stat->iter_idx_n) == 1)
	{
		fprintf(stderr,"%d\n",pl->read_stat->iter_idx_n);
		for (i = 0; i < pl->read_stat->iter_idx_n; i++)
		{
			if (fscanf(fp, "%d ", &pl->read_stat->iter_idx[i]) == 1)
			{
				// fprintf(stderr,"%d\n",pl->read_stat->iter_idx[i]);
			}
			else exit(1);
		}
	}
	else exit(1);

	if (fscanf(fp, "%d\n", &pl->n_total_read) == 1)
	{
		fprintf(stderr,"%d\n",pl->n_total_read);
		for (i = 0; i < pl->n_total_read; i++)
		{
			if (fscanf(fp, "%d ", &pl->read_stat->is_idx[i]) == 1)
			{
				// fprintf(stderr,"%d\n",pl->read_stat->iter_idx[i]);
			}
			else exit(1);
		}
	}
	else exit(1);

	if (fscanf(fp, "%d\n", &pl->read_stat->iter_map_n) == 1)
	{
		fprintf(stderr,"%d\n",pl->read_stat->iter_map_n);
		for (i = 0; i < pl->read_stat->iter_map_n; i++)
		{
			if (fscanf(fp, "%d ", &pl->read_stat->iter_map[i]) == 1)
			{
				// fprintf(stderr,"%d\n",pl->read_stat->iter_map[i]);
			}
			else exit(1);
		}
	}
	else exit(1);

	if (fscanf(fp, "%d\n", &pl->n_total_read) == 1)
	{
		fprintf(stderr,"%d\n",pl->n_total_read);
		for (i = 0; i < pl->n_total_read; i++)
		{
			if (fscanf(fp, "%d ", &pl->read_stat->is_ove[i]) == 1)
			{
				// fprintf(stderr,"%d\n",pl->read_stat->iter_idx[i]);
			}
			else exit(1);
		}
	}
	else exit(1);

	if (fscanf(fp, "%ld %ld\n", &pl->file_idx.n, &pl->file_idx.m) == 2)
	{
		fprintf(stderr, "%ld %ld\n", pl->file_idx.n, pl->file_idx.m);
		for (i = 0; i < pl->file_idx.n; ++i)
		{
			if (fscanf(fp, "%d ", &pl->file_idx.offs[i]) == 1)
			{
				// fprintf(stderr,"%d\n",pl->file_idx.offs[i]);
			}
			else exit(1);
		}
	}
	else exit(1);

	if (fscanf(fp, "%d %d %d %d\n", &pl->file_add_p.last_rid, &pl->file_add_p.upper_v, &pl->file_add_p.add_n, &pl->file_add_p.add_m) == 4)
	{
		fprintf(stderr, "%d %d %d %d\n", pl->file_add_p.last_rid, pl->file_add_p.upper_v, pl->file_add_p.add_n, pl->file_add_p.add_m);
		pl->file_add_p.add_p = (uint32_t*)calloc(pl->file_add_p.add_m, sizeof(uint32_t));
		for (i = 0; i < pl->file_add_p.add_n; ++i)
		{
			if (fscanf(fp, "%d ", &pl->file_add_p.add_p[i]) == 1)
			{
				fprintf(stderr,"%d\n",pl->file_add_p.add_p[i]);
			}
			else exit(1);
		}
	}
	else exit(1);

	fclose(fp);
	return 0;
}
#endif
OVE_C_t* load_paf_file(const char *temp_batch_dir, OVE_C_t *ove_cl, uint32_t *read_n)
{
    paf_file_t *fp;
    paf_rec_t r;
	uint32_t qid, tid;
	uint32_t new_read_n, i;

    fp = paf_open(temp_batch_dir);
    if (!fp)
    {
        fprintf(stderr, "[%s] could not open PAF file %s\n", __func__, temp_batch_dir);
        exit(1);
    }

    while (paf_read(fp, &r) >= 0) //r: a paf row
    {
		qid = r.qid;
		tid = r.tid;
		if (qid >= (*read_n) || tid >= (*read_n))
		{
			new_read_n = (*read_n) == 0 ? 65536 : (*read_n) << 1;
			ove_cl = (OVE_C_t *)realloc(ove_cl, new_read_n * sizeof(OVE_C_t));
			for (i = (*read_n); i < new_read_n; i++)
			{
				ove_cl[i].n = 0;
				ove_cl[i].m = 2;
				ove_cl[i].ove = (OVE_t *)calloc(ove_cl[i].m, sizeof(OVE_t));
			}
			printf("read n %d->%d\n",(*read_n),new_read_n);
			(*read_n) = new_read_n;
		}
        if (ove_cl[qid].n >= ove_cl[qid].m)
        {
            ove_cl[qid].m = ove_cl[qid].m == 0 ? 2 : ove_cl[qid].m << 1;
            ove_cl[qid].ove = (OVE_t *)realloc(ove_cl[qid].ove, ove_cl[qid].m * sizeof(OVE_t));
        }
		// printf("%d %d\n",qid,tid);
        ove_cl[qid].ove[ove_cl[qid].n].qid = r.qid; ove_cl[qid].ove[ove_cl[qid].n].tid = r.tid;
        ove_cl[qid].ove[ove_cl[qid].n].ql = r.ql; ove_cl[qid].ove[ove_cl[qid].n].tl = r.tl;
        ove_cl[qid].ove[ove_cl[qid].n].qs = r.qs; ove_cl[qid].ove[ove_cl[qid].n].qe = r.qe;
        ove_cl[qid].ove[ove_cl[qid].n].ts = r.ts; ove_cl[qid].ove[ove_cl[qid].n].te = r.te;
        ove_cl[qid].ove[ove_cl[qid].n].rev = r.rev; ove_cl[qid].ove[ove_cl[qid].n].mbp = r.mb;
		ove_cl[qid].ove[ove_cl[qid].n].mln = r.ml; ove_cl[qid].ove[ove_cl[qid].n++].score = r.ms;
    }
    // printf("\n");

    paf_close(fp);
    return ove_cl;
}

uint32_t generate_transitive_overlaps_file_core(const char *temp_file_perfix, pipeline_t *pl)
{
	uint32_t r_i;
	FILE *fp_trans_ove = NULL;
	char temp_trans_ove_dir[1024];
	memset(temp_trans_ove_dir, 0, 1024);
	if (temp_file_perfix == NULL)
	{
		strcpy(temp_trans_ove_dir, "./trans.paf");
	}
	else
	{
		strcpy(temp_trans_ove_dir, temp_file_perfix);
		strcat(temp_trans_ove_dir, "_trans.paf");
	}

	fp_trans_ove = fopen(temp_trans_ove_dir, "w");
	if (fp_trans_ove == NULL)
	{
		fprintf(stderr, "[%s Wrong] Failed to open file %s!!!\n",  __func__, temp_trans_ove_dir);
		exit(1);
	}
	fprintf(stderr, "[%s] %d transitive iterations...\n",  __func__, trans_iter);
	fprintf(stderr, "[%s] Output transitive overlaps to file %s...\n",  __func__, temp_trans_ove_dir);

	uint32_t *visited_rid = (uint32_t *)calloc(pl->n_total_read, sizeof(uint32_t));
	// for (r_i = 80000000; r_i < pl->n_total_read; ++r_i)
	for (r_i = 0; r_i < pl->n_total_read; ++r_i)
		generate_transitive_edges_of_graph(pl, r_i, fp_trans_ove, visited_rid);
	if (visited_rid != NULL) {free(visited_rid); visited_rid = NULL;}

	fprintf(stderr, "[%s] Transitive overlaps generated...\n",  __func__);
	fclose(fp_trans_ove);
	return 0;
}

uint32_t generate_transitive_overlaps_file(const char *trans_file)
{
	pipeline_t pl;
	uint32_t r_i, o_i;
	pl.n_total_read = 65536;
	pl.ove_cl = (OVE_C_t *)calloc(pl.n_total_read, sizeof(OVE_C_t));
	for (r_i = 0; r_i < pl.n_total_read; r_i++)
	{
		pl.ove_cl[r_i].n = 0;
		pl.ove_cl[r_i].m = 2;
		pl.ove_cl[r_i].ove = (OVE_t *)calloc(pl.ove_cl[r_i].m, sizeof(OVE_t));
	}

	pl.ove_cl = load_paf_file(trans_file, pl.ove_cl, &pl.n_total_read);
	fprintf(stderr, "[%s] Loaded read overlaps...\n", __func__);

	pl.symm = (uint32_t *)calloc(pl.n_total_read, sizeof(uint32_t));
	for (r_i = 0; r_i < pl.n_total_read; ++r_i) pl.symm[r_i] = pl.ove_cl[r_i].n;
	pl.read_stat = (read_stat_t *)calloc(1, sizeof(read_stat_t));
	pl.read_stat->rlen = NULL;
	pl.read_stat->is_idx = NULL;
	pl.read_stat->rlen = (uint32_t *)calloc(pl.n_total_read, sizeof(uint32_t));
	pl.read_stat->is_idx = (uint32_t *)calloc(pl.n_total_read, sizeof(uint32_t));
	for (r_i = 0; r_i < pl.n_total_read; ++r_i)
	{
		if (pl.ove_cl[r_i].n == 0) continue;
		for (o_i = 0; o_i < pl.ove_cl[r_i].n; ++o_i)
		{
			pl.read_stat->is_idx[pl.ove_cl[r_i].ove[o_i].tid] = 1;
			pl.read_stat->rlen[pl.ove_cl[r_i].ove[o_i].tid] = pl.ove_cl[r_i].ove[o_i].tl;
			pl.read_stat->rlen[pl.ove_cl[r_i].ove[o_i].qid] = pl.ove_cl[r_i].ove[o_i].ql;
		}
	}
	construct_symmetrical_graph(pl.n_total_read, pl.symm, pl.ove_cl, pl.read_stat);
	fprintf(stderr, "[%s] Symm read overlaps...\n", __func__);

	int filename_len = strlen(trans_file);
	char *file_prefix = (char *)calloc(filename_len, sizeof(char));
	char gz_type = trans_file[filename_len-1];
	if (gz_type == 'z') strncpy(file_prefix,trans_file,filename_len-7);
	else strncpy(file_prefix,trans_file,filename_len-4);
	char *str_trans_iter = (char *)calloc(4, sizeof(char));
	sprintf(str_trans_iter, "%d", trans_iter);
	strcat(file_prefix, str_trans_iter);

	generate_transitive_overlaps_file_core(file_prefix, &pl);

	for (r_i = 0; r_i < pl.n_total_read; r_i++)
	{
		if (pl.ove_cl[r_i].ove != NULL) {free(pl.ove_cl[r_i].ove); pl.ove_cl[r_i].ove = NULL;}
	}
	if (pl.ove_cl != NULL) {free(pl.ove_cl); pl.ove_cl = NULL;}
	if (pl.symm != NULL) {free(pl.symm); pl.symm = NULL;}
	if (pl.read_stat->rlen != NULL) {free(pl.read_stat->rlen); pl.read_stat->rlen = NULL;}
	if (pl.read_stat->is_idx != NULL) {free(pl.read_stat->is_idx); pl.read_stat->is_idx = NULL;}
	if (pl.read_stat != NULL) {free(pl.read_stat); pl.read_stat = NULL;}
	return 0;
}

uint32_t finding_overlapping(const char *read_fastq, const char *index_fastq, const char *temp_file_perfix)
{
	uint32_t i, r_i, o_i;
	// step 1: sorting reads' length
	bseq_file_t *bf = NULL;
	bf = bseq_open(read_fastq);
	if (bf == NULL)
	{
		fprintf(stderr, "[Warning] Wrong input file route or name:%s \n", read_fastq);
		exit(1);
	}

	pipeline_t pl;
	pl.file_idx.n = 0;
	pl.file_idx.m = 256;
	pl.file_idx.offs = (int*)calloc(pl.file_idx.m, sizeof(int));
	pl.x_len = 0;
	pl.max_len = 0;
	pl.batch_i = 0;
	pl.n_processed = 0;
	pl.n_total_read = 0;
	pl.index_batch_size = 1000000;
	pl.fp = bf;
	pl.read_stat = (read_stat_t *)calloc(1, sizeof(read_stat_t));
	pl.read_stat->rlen = NULL;
	pl.read_stat->is_idx = NULL;
	pl.read_stat->iter_idx = NULL;
	pl.read_stat->iter_idx_n = 0;
	pl.read_stat->iter_map = NULL;
	pl.read_stat->iter_map_n = 0;
	kt_pipeline(thread_n < 2 ? thread_n : 2, sort_read_length_pipeline, &pl, 2);

	pl.x_len = pl.x_len / pl.batch_i;
	pl.n_total_read = pl.n_processed;
	if (bf != NULL) bseq_close(bf);

	// TODO: a better method to estimate the number of vertex to replace (pl.x_len * rep)
	double rep = thread_n / 8;
	batch_size_base = (double)(memory - (double)pl.n_total_read  / 1024 / 1024 * (58 + 96 * top_n + 2 * (double)pl.x_len / split_len) / 1024 - (double)thread_n / 1024 * (20 + 16 * pl.x_len + 240 * pl.x_len * rep) / 1024 / 1024) * (double)pl.x_len / (double)(19 * pl.x_len + thread_n * 48 * top_n) * 1024 * 1024 * 1024;
	// batch_size_base = (memory / (sizeof(MINIMIZER_t) + 3.0 + 2.0 / thread_n) * 1024 * 1024 * 1024) - (pl.n_total_read*(sizeof(OVE_C_t)+top_n*sizeof(OVE_t)+5*sizeof(uint32_t)) / (sizeof(MINIMIZER_t) + 3.0 + 2.0 / thread_n));
	fprintf(stderr, "[Batch Size: %.3fs, %.3fGB] Total %d reads, modify batch size %ld\n", realtime() - realtime0, peak_memory() / 1024.0 / 1024.0, pl.n_total_read, batch_size_base);
	fprintf(stderr, "[Index: %.3fs, %.3fGB] Sorted %d read length, indexing read > %d, max read length %d\n", realtime() - realtime0, peak_memory() / 1024.0 / 1024.0, pl.n_total_read, pl.x_len, pl.max_len);

	// step 2: iteratively indexing and mapping
	pl.iter = 0;
	pl.idx_batch = 0;
	pl.n_threads = thread_n;
	pl.query_batch_size = batch_size_base;
	pl.m_b = m_b;
	pl.min_ove = min_ove;
	pl.max_hang = max_hang;
	pl.top_n = top_n;
	pl.search_step = s_s;
	pl.split_len = split_len;

	pl.read_stat->is_idx = (uint32_t *)calloc(pl.n_total_read, sizeof(uint32_t));
	pl.read_stat->iter_idx = (uint32_t *)calloc(pl.n_total_read, sizeof(uint32_t));
	pl.read_stat->iter_map = (uint32_t *)calloc(pl.n_total_read, sizeof(uint32_t));
	pl.read_stat->is_ove = (uint32_t *)calloc(pl.n_total_read, sizeof(uint32_t));
	for (r_i = 0; r_i < pl.n_total_read; ++r_i)
	{
		if (pl.read_stat->rlen[r_i] >= pl.x_len)
		{
			pl.read_stat->is_idx[r_i] = 1;
			pl.read_stat->iter_idx[pl.read_stat->iter_idx_n++] = r_i; // mark seed read for the first iteration
		}
		pl.read_stat->iter_map[pl.read_stat->iter_map_n++] = r_i; // mark seed read for the first iteration
	}
	// fprintf(stderr, "[%s: %.3fs, %.3fGB] calloc %fGB pl.read_stat memory...\n", __func__, realtime() - realtime0, peak_memory() / 1024.0 / 1024.0, pl.n_total_read * (2 * sizeof(uint32_t) + 2 * sizeof(uint32_t)) / 1024.0 / 1024.0 / 1024.0);

	pl.read_cov = (read_cov_t *)calloc(pl.n_total_read, sizeof(read_cov_t));
	pl.ave_cov = (double *)calloc(pl.n_total_read, sizeof(double));
	pl.blkn = (uint16_t *)calloc(pl.n_total_read, sizeof(uint16_t));
	for (r_i = 0; r_i < pl.n_total_read; ++r_i)
	{
		pl.ave_cov[r_i] = 0.;
		pl.blkn[r_i] = pl.read_stat->rlen[r_i] / split_len + 1;
		pl.read_cov[r_i].cov = (uint16_t *)calloc(pl.blkn[r_i], sizeof(uint16_t));
	}
	// fprintf(stderr, "[%s: %.3fs, %.3fGB] calloc %fGB pl.read_cov memory...\n", __func__, realtime() - realtime0, peak_memory() / 1024.0 / 1024.0, pl.n_total_read * (sizeof(read_cov_t) + (pl.x_len / split_len + 1) * sizeof(uint16_t) + sizeof(uint16_t) + sizeof(double)) / 1024.0 / 1024.0 / 1024.0);
	pl.visited_rid = (uint32_t **)calloc(pl.n_threads, sizeof(uint32_t*));
	for (r_i = 0; r_i < pl.n_threads; ++r_i)
	{
		pl.visited_rid[r_i] = (uint32_t *)calloc(pl.n_total_read, sizeof(uint32_t));
	}

	int is_next = 1;
	uint32_t extra_ove_num = 0;
	uint32_t extra_ove_num2 = 0;
	bseq_file_t *fp = NULL;
	pl.file_add_p.upper_v = INT32_MAX - pl.max_len * 3;

	while(is_next)
	{
		is_next = 0;
		pl.file_add_p.last_rid = 0;
		pl.file_add_p.add_n = 1;
		pl.file_add_p.add_m = 2;
		pl.file_add_p.add_p = (uint32_t*)calloc(pl.file_add_p.add_m, sizeof(uint32_t));
		pl.index_n_sta = 0;
		pl.index_idx_n_sta = 0;
		pl.index_idx_n_end = 0;
		pl.ove_cl = (OVE_C_t *)calloc(pl.n_total_read, sizeof(OVE_C_t));
		if (pl.ove_cl == NULL)
		{
			fprintf(stderr, "[%s] calloc %ldGB pl.ove_cl memory error!\n", __func__, pl.n_total_read * sizeof(OVE_C_t) / 1024 / 1024 / 1024);
			exit(1);
		}
		for (i = 0; i < pl.n_total_read; i++)
		{
			pl.ove_cl[i].n = 0;
			pl.ove_cl[i].m  = pl.top_n;
			pl.ove_cl[i].ove = (OVE_t *)calloc(pl.ove_cl[i].m, sizeof(OVE_t));
		}
		pl.symm = (uint32_t *)calloc(pl.n_total_read, sizeof(uint32_t));

		for(;;)
		{
			// step 2.1: indexing reads in read_stat->iter_idx
			mini_idx_t *mi = NULL;
			fp = bseq_open(index_fastq);
			pl.index_idx_n_sta = pl.index_idx_n_end;
			mi = indexing_input_read(fp, &pl.rep_n, &pl.index_n_sta, &pl.index_idx_n_end, pl.read_stat, &pl.file_idx, &pl.file_add_p);
			if (fp != NULL) bseq_close(fp);
			// fprintf(stderr, "index idx sta %d, end %d\n",pl.index_idx_n_sta,pl.index_idx_n_end);

			if (mi->mm.mi_n == 0)
			{
				if (mi->mm.mi_v != NULL) {free(mi->mm.mi_v); mi->mm.mi_v = NULL;}
				if (mi->mi_count != NULL) {free(mi->mi_count); mi->mi_count = NULL;}
				if (mi != NULL) {free(mi); mi = NULL;}
				break;
			}
			pl.idx_batch++;

			// step 2.2: mapping reads to indexed reads
			bf = bseq_open(read_fastq);
			pl.mi = mi;
			pl.fp = bf;
			pl.map_batch = 0;
			pl.n_processed = 0;
			pl.map_n_sta = 0;
			pl.map_idx_n_sta = 0;
			pl.map_idx_n_end = 0;

			kt_pipeline(thread_n < 3 ? thread_n : 3, mapping_pipeline, &pl, 3);

			if (bf != NULL)	bseq_close(bf);

			if (mi->mm.mi_v != NULL) {free(mi->mm.mi_v); mi->mm.mi_v = NULL;}
			if (mi->mi_count != NULL) {free(mi->mi_count); mi->mi_count = NULL;}
			if (mi != NULL) {free(mi); mi = NULL;}
		}

		if (pl.iter > 0)
		{
			char temp_iter_dir[1024];
			memset(temp_iter_dir, 0, 1024);

			char iter[64] = {0};
			sprintf(iter, "%d", pl.iter-1);
			if (temp_file_perfix == NULL)
			{
				strcpy(temp_iter_dir, "./oves_iter_");
				strcat(temp_iter_dir, iter);
				strcat(temp_iter_dir, ".paf");
			}
			else
			{
				strcpy(temp_iter_dir, temp_file_perfix);
				strcat(temp_iter_dir, "_oves_iter_");
				strcat(temp_iter_dir, iter);
				strcat(temp_iter_dir, ".paf");
			}
			fprintf(stderr, "[Coverage: %.3fs, %.3fGB] the intermediate result loaded from %s\n", realtime() - realtime0, peak_memory() / 1024.0 / 1024.0, temp_iter_dir);
			load_paf_file(temp_iter_dir, pl.ove_cl, &pl.n_total_read);
		}
		construct_symmetrical_graph(pl.n_total_read, pl.symm, pl.ove_cl, pl.read_stat);

		for (r_i = 0; r_i < pl.n_total_read; ++r_i)
		{
			pl.ave_cov[r_i] = 0.;
			memset(pl.read_cov[r_i].cov, 0, pl.blkn[r_i] * sizeof(uint16_t));
		}
		for (r_i = 0; r_i < pl.n_threads; ++r_i)
		{
			memset(pl.visited_rid[r_i], 0, pl.n_total_read * sizeof(uint32_t));
		}
		
		kt_for(thread_n > 8 ? 8 : thread_n, estimate_coverage_of_seed_blocks, &pl, pl.n_total_read);
		kt_for(thread_n > 8 ? 8 : thread_n, estimate_coverage_of_read_blocks, &pl, pl.n_total_read);

		qsort(pl.ave_cov, pl.n_total_read, sizeof(double), compare_double);
		pl.ave_median = pl.ave_cov[(int)(pl.n_total_read*0.5)];
		fprintf(stderr, "[Coverage: %.3fs, %.3fGB] median coverage of all reads: %f...\n", realtime() - realtime0, peak_memory() / 1024.0 / 1024.0, pl.ave_median);
		is_next = choose_idx_read_for_next_iteration(&pl, pl.read_cov, pl.ave_median);

		// map reads without overlap to the last indexed reads
		if (pl.read_stat->iter_idx_n == 0) is_next = 0;
		if (is_next == 0)
		{
			pl.read_stat->iter_map_n = 0;
			for (r_i = 0; r_i < pl.n_total_read; ++r_i)
			{
				if (pl.read_stat->is_ove[r_i] == 1 && pl.ove_cl[r_i].n == 0)
				{
					extra_ove_num++;
					pl.read_stat->iter_map[pl.read_stat->iter_map_n++] = r_i;
				}
				if (pl.read_stat->is_ove[r_i] == 1)
				{
					extra_ove_num2++;
				}
			}
			fprintf(stderr, "[Extra Iteration: %.3fs, %.3fGB] %d reads are extract to extra iteration...\n", realtime() - realtime0, peak_memory() / 1024.0 / 1024.0, extra_ove_num);
			fprintf(stderr, "[Extra Iteration: %.3fs, %.3fGB] %d reads are marked as repetitive reads...\n", realtime() - realtime0, peak_memory() / 1024.0 / 1024.0, extra_ove_num2);
		}
	
#ifdef INTER_RESULT
		FILE *fp_iter = NULL;
		char temp_iter_dir[1024];
		memset(temp_iter_dir, 0, 1024);

		FILE *fp_batch = NULL;
		char temp_batch_dir[1024];
		memset(temp_batch_dir, 0, 1024);

		char iter[64] = {0};
		sprintf(iter, "%d", pl.iter);
		if (temp_file_perfix == NULL)
		{
			strcpy(temp_iter_dir, "./oves_iter_");
			strcat(temp_iter_dir, iter);
			strcat(temp_iter_dir, ".paf");
			strcpy(temp_batch_dir, "./oves_iter_");
			strcat(temp_batch_dir, iter);
			strcat(temp_batch_dir, ".info");
		}
		else
		{
			strcpy(temp_iter_dir, temp_file_perfix);
			strcat(temp_iter_dir, "_oves_iter_");
			strcat(temp_iter_dir, iter);
			strcat(temp_iter_dir, ".paf");
			strcpy(temp_batch_dir, temp_file_perfix);
			strcat(temp_batch_dir, "_oves_iter_");
			strcat(temp_batch_dir, iter);
			strcat(temp_batch_dir, ".info");
		}
		fprintf(stderr, "[To next Iteration: %.3fs, %.3fGB] The intermediate result for iteration %d output to %s\n", realtime() - realtime0, peak_memory() / 1024.0 / 1024.0, pl.iter, temp_iter_dir);


		fp_iter = fopen(temp_iter_dir, "w");
		if (fp_iter == NULL) exit(1);

		for (r_i = 0; r_i < pl.n_total_read; ++r_i)
		{
			if (pl.ove_cl[r_i].n == 0) continue;
			for (o_i = 0; o_i < pl.ove_cl[r_i].n; ++o_i)
			{
				fprintf(fp_iter, "%u\t%u\t%u\t%u\t%c\t", pl.ove_cl[r_i].ove[o_i].qid, pl.ove_cl[r_i].ove[o_i].ql, pl.ove_cl[r_i].ove[o_i].qs, pl.ove_cl[r_i].ove[o_i].qe, "+-"[pl.ove_cl[r_i].ove[o_i].rev]);
				fprintf(fp_iter, "%u\t%u\t%u\t%u\t", pl.ove_cl[r_i].ove[o_i].tid, pl.ove_cl[r_i].ove[o_i].tl, pl.ove_cl[r_i].ove[o_i].ts, pl.ove_cl[r_i].ove[o_i].te);
				fprintf(fp_iter, "%d\t%d\t%.2f\ti:%c\n", pl.ove_cl[r_i].ove[o_i].mbp, pl.ove_cl[r_i].ove[o_i].mln, (double)(int)pl.ove_cl[r_i].ove[o_i].score, "mi"[pl.read_stat->is_idx[r_i]]);
			}
		}

		fclose(fp_iter);

		// printf other batch information
		fp_batch = fopen(temp_batch_dir, "w");
		if (fp_batch == NULL) exit(1);

		fprintf(fp_batch,"%d %d %d %d\n",pl.iter,pl.index_n_sta,pl.index_idx_n_sta,pl.index_idx_n_end);
		fprintf(fp_batch,"%d\n", pl.read_stat->iter_idx_n);
		for (r_i = 0; r_i < pl.read_stat->iter_idx_n; ++r_i)	fprintf(fp_batch,"%d\t", pl.read_stat->iter_idx[r_i]);
		fprintf(fp_batch,"\n");
		fprintf(fp_batch,"%d\n", pl.n_total_read);
		for (r_i = 0; r_i < pl.n_total_read; ++r_i)	fprintf(fp_batch,"%d\t", pl.read_stat->is_idx[r_i]);
		fprintf(fp_batch,"\n");
		fprintf(fp_batch,"%d\n", pl.read_stat->iter_map_n);
		for (r_i = 0; r_i < pl.read_stat->iter_map_n; ++r_i)	fprintf(fp_batch,"%d\t", pl.read_stat->iter_map[r_i]);
		fprintf(fp_batch,"\n");
		fprintf(fp_batch,"%d\n", pl.n_total_read);
		for (r_i = 0; r_i < pl.n_total_read; ++r_i)	fprintf(fp_batch,"%d\t", pl.read_stat->is_ove[r_i]);
		fprintf(fp_batch,"\n");
		fprintf(fp_batch,"%ld %ld\n", pl.file_idx.n, pl.file_idx.m);
		for (r_i = 0; r_i < pl.file_idx.n; ++r_i)	fprintf(fp_batch,"%d\t", pl.file_idx.offs[r_i]);
		fprintf(fp_batch,"\n");
		fprintf(fp_batch,"%d %d %d %d\n", pl.file_add_p.last_rid, pl.file_add_p.upper_v, pl.file_add_p.add_n, pl.file_add_p.add_m);
		for (r_i = 0; r_i < pl.file_add_p.add_n; ++r_i)	fprintf(fp_batch,"%d\t", pl.file_add_p.add_p[r_i]);
		fprintf(fp_batch,"\n");

		fclose(fp_batch);
		fprintf(stderr, "[To next Iteration: %.3fs, %.3fGB] The intermediate result for batch (%d,%d) output to %s\n", realtime() - realtime0, peak_memory() / 1024.0 / 1024.0, pl.idx_batch, pl.map_batch, temp_batch_dir);

#endif
		for (i = 0; i < pl.n_total_read; i++)
		{
			if (pl.ove_cl[i].ove != NULL) {free(pl.ove_cl[i].ove); pl.ove_cl[i].ove = NULL;}
		}
		if (pl.ove_cl != NULL) {free(pl.ove_cl); pl.ove_cl = NULL;}
		if (pl.symm != NULL) {free(pl.symm); pl.symm = NULL;}
		if (pl.file_add_p.add_p != NULL) {free(pl.file_add_p.add_p); pl.file_add_p.add_p = NULL;}

		pl.iter++;
	}

	pl.ove_cl = (OVE_C_t *)calloc(pl.n_total_read, sizeof(OVE_C_t));
	if (pl.ove_cl == NULL)
	{
		fprintf(stderr, "[%s] calloc %ldGB pl.ove_cl memory error!\n", __func__, pl.n_total_read * sizeof(OVE_C_t) / 1024 / 1024 / 1024);
		exit(1);
	}
	for (i = 0; i < pl.n_total_read; i++)
	{
		pl.ove_cl[i].n = 0;
		pl.ove_cl[i].m  = pl.top_n;
		pl.ove_cl[i].ove = (OVE_t *)calloc(pl.ove_cl[i].m, sizeof(OVE_t));
	}
	pl.symm = (uint32_t *)calloc(pl.n_total_read, sizeof(uint32_t));

	if (extra_ove_num > 0)
	{
		pl.search_step = s_s*7.5;
		pl.file_add_p.last_rid = 0;
		pl.file_add_p.add_n = 1;
		pl.file_add_p.add_m = 2;
		pl.file_add_p.add_p = (uint32_t*)calloc(pl.file_add_p.add_m, sizeof(uint32_t));
		pl.index_n_sta = 0;
		pl.index_idx_n_sta = 0;
		pl.index_idx_n_end = 0;
		rep = thread_n;
		batch_size_base = (double)(memory - (double)pl.n_total_read  / 1024 / 1024 * (58 + 96 * top_n + 2 * (double)pl.x_len / split_len) / 1024 - (double)thread_n / 1024 * (20 + 16 * pl.x_len + 240 * pl.x_len * rep) / 1024 / 1024) * (double)pl.x_len / (double)(19 * pl.x_len + thread_n * 48 * top_n) * 1024 * 1024 * 1024;
		fprintf(stderr, "[Batch Size] Total %d reads, modify batch size %ld\n", pl.n_total_read, batch_size_base);
		pl.query_batch_size = batch_size_base;

		for(;;)
		{
			// step 2.1: indexing reads in read_stat->iter_idx
			mini_idx_t *mi = NULL;
			fp = bseq_open(index_fastq);
			pl.index_idx_n_sta = pl.index_idx_n_end;
			mi = indexing_input_read(fp, &pl.rep_n, &pl.index_n_sta, &pl.index_idx_n_end, pl.read_stat, &pl.file_idx, &pl.file_add_p);
			if (fp != NULL) bseq_close(fp);
			// fprintf(stderr, "index idx sta %d, end %d\n",pl.index_idx_n_sta,pl.index_idx_n_end);
			if (mi->mm.mi_n == 0)
			{
				if (mi->mm.mi_v != NULL) {free(mi->mm.mi_v); mi->mm.mi_v = NULL;}
				if (mi->mi_count != NULL) {free(mi->mi_count); mi->mi_count = NULL;}
				if (mi != NULL) {free(mi); mi = NULL;}
				break;
			}
			pl.idx_batch++;

			// step 2.2: mapping reads to indexed reads
			bf = bseq_open(read_fastq);
			pl.mi = mi;
			pl.fp = bf;
			pl.map_batch = 0;
			pl.n_processed = 0;
			pl.map_n_sta = 0;
			pl.map_idx_n_sta = 0;
			pl.map_idx_n_end = 0;

			kt_pipeline(thread_n < 3 ? thread_n : 3, mapping_pipeline, &pl, 3);

			if (bf != NULL)	bseq_close(bf);

			if (mi->mm.mi_v != NULL) {free(mi->mm.mi_v); mi->mm.mi_v = NULL;}
			if (mi->mi_count != NULL) {free(mi->mi_count); mi->mi_count = NULL;}
			if (mi != NULL) {free(mi); mi = NULL;}
		}
		if (pl.file_add_p.add_p != NULL) {free(pl.file_add_p.add_p); pl.file_add_p.add_p = NULL;}
	}

	if (pl.iter > 0)
	{
		char temp_iter_dir[1024];
		memset(temp_iter_dir, 0, 1024);

		char iter[64] = {0};
		sprintf(iter, "%d", pl.iter-1);
		if (temp_file_perfix == NULL)
		{
			strcpy(temp_iter_dir, "./oves_iter_");
			strcat(temp_iter_dir, iter);
			strcat(temp_iter_dir, ".paf");
		}
		else
		{
			strcpy(temp_iter_dir, temp_file_perfix);
			strcat(temp_iter_dir, "_oves_iter_");
			strcat(temp_iter_dir, iter);
			strcat(temp_iter_dir, ".paf");
		}
		fprintf(stderr, "[Current iteration: %.3fs, %.3fGB] The intermediate result loaded from %s\n", realtime() - realtime0, peak_memory() / 1024.0 / 1024.0, temp_iter_dir);
		load_paf_file(temp_iter_dir, pl.ove_cl, &pl.n_total_read);
		construct_symmetrical_graph(pl.n_total_read, pl.symm, pl.ove_cl, pl.read_stat);
	}

	// step 4: output
	FILE *fp_ove = NULL;
	char temp_ove_dir[1024];
	memset(temp_ove_dir, 0, 1024);
	if (temp_file_perfix != NULL)
	{
		strcpy(temp_ove_dir, temp_file_perfix);
		strcat(temp_ove_dir, "_oves.paf");
		fp_ove = freopen(temp_ove_dir, "w", stdout);
		// fp_ove = fopen(temp_ove_dir, "w");
		if (fp_ove == NULL)
		{
			fprintf(stderr, "[%s Wrong] Failed to open file %s!!!\n",  __func__, temp_ove_dir);
			exit(1);
		}
	}

	// gz file
	char file_type = 0;
	char *file_suff = NULL;
	int filename_len = strlen(read_fastq);
	char gz_type = read_fastq[filename_len-1];
	// fasta/q file
	if (gz_type == 'z') file_type = read_fastq[filename_len-4];
	else file_type = read_fastq[filename_len-1];
	if (file_type == 'q')	file_suff = "fastq";
	else file_suff = "fasta";

	FILE *fp_un_ove = NULL;
	char temp_un_ove_dir[1024];
	memset(temp_un_ove_dir, 0, 1024);
	if (temp_file_perfix == NULL)
	{
		strcpy(temp_un_ove_dir, "./un_oves.");
		strcat(temp_un_ove_dir, file_suff);
	}
	else
	{
		strcpy(temp_un_ove_dir, temp_file_perfix);
		strcat(temp_un_ove_dir, "_un_oves.");
		strcat(temp_un_ove_dir, file_suff);
	}

	// look up through the raw read file
	fp_un_ove = fopen(temp_un_ove_dir, "w");
	if (fp_un_ove == NULL)
	{
		fprintf(stderr, "[%s Wrong] Failed to open file %s!!!\n",  __func__, temp_un_ove_dir);
		exit(1);
	}
	fp = bseq_open(read_fastq);
	int un_ove_read_num = bseq_rewrite_read(fp, fp_un_ove, pl.symm, &pl.file_idx, file_type);
	fprintf(stderr, "[Result: %.3fs, %.3fGB] Rewrite %d reads without overlap to file %s...\n", realtime() - realtime0, peak_memory() / 1024.0 / 1024.0, un_ove_read_num, temp_un_ove_dir);
	if (fp != NULL) bseq_close(fp);
	fclose(fp_un_ove);

	uint32_t total_index_num = 0;
	for (r_i = 0; r_i < pl.n_total_read; ++r_i)
	{
		if (pl.read_stat->is_idx[r_i] == 1) total_index_num += 1;
		if (pl.ove_cl[r_i].n == 0) continue;
		for (o_i = 0; o_i < pl.ove_cl[r_i].n; ++o_i)
		{
			printf("%u\t%u\t%u\t%u\t%c\t", pl.ove_cl[r_i].ove[o_i].qid, pl.ove_cl[r_i].ove[o_i].ql, pl.ove_cl[r_i].ove[o_i].qs, pl.ove_cl[r_i].ove[o_i].qe, "+-"[pl.ove_cl[r_i].ove[o_i].rev]);
			printf("%u\t%u\t%u\t%u\t", pl.ove_cl[r_i].ove[o_i].tid, pl.ove_cl[r_i].ove[o_i].tl, pl.ove_cl[r_i].ove[o_i].ts, pl.ove_cl[r_i].ove[o_i].te);
			printf("%d\t%d\t%.2f\ti:%c\n", pl.ove_cl[r_i].ove[o_i].mbp, pl.ove_cl[r_i].ove[o_i].mln, pl.ove_cl[r_i].ove[o_i].score, "mi"[pl.read_stat->is_idx[r_i]]);
		}
	}
	if (temp_file_perfix != NULL) fclose(fp_ove);
	fprintf(stderr, "[Result: %.3fs, %.3fGB] Total indexed read num: %d...\n", realtime() - realtime0, peak_memory() / 1024.0 / 1024.0, total_index_num);

	if (trans_iter > 0)	generate_transitive_overlaps_file_core(temp_file_perfix, &pl);

	if (pl.symm != NULL) {free(pl.symm); pl.symm = NULL;}
	for (i = 0; i < pl.n_total_read; i++)
	{
		if (pl.ove_cl[i].ove != NULL) {free(pl.ove_cl[i].ove); pl.ove_cl[i].ove = NULL;}
	}
	if (pl.ove_cl != NULL) {free(pl.ove_cl); pl.ove_cl = NULL;}

	for (r_i = 0; r_i < pl.n_threads; r_i++)
	{
		if (pl.visited_rid[r_i] != NULL) {free(pl.visited_rid[r_i]); pl.visited_rid[r_i] = NULL;}
	}
	if (pl.visited_rid != NULL) {free(pl.visited_rid); pl.visited_rid = NULL;}

	for (r_i = 0; r_i < pl.n_total_read; r_i++)
	{
		if (pl.read_cov[r_i].cov != NULL) {free(pl.read_cov[r_i].cov); pl.read_cov[r_i].cov = NULL;}
	}
	if (pl.read_cov != NULL) {free(pl.read_cov); pl.read_cov = NULL;}
	if (pl.ave_cov != NULL) {free(pl.ave_cov); pl.ave_cov = NULL;}
	if (pl.blkn != NULL) {free(pl.blkn); pl.blkn = NULL;}

	if (pl.read_stat->rlen != NULL) {free(pl.read_stat->rlen); pl.read_stat->rlen = NULL;}
	if (pl.read_stat->is_idx != NULL) {free(pl.read_stat->is_idx); pl.read_stat->is_idx = NULL;}
	if (pl.read_stat->iter_idx != NULL) {free(pl.read_stat->iter_idx); pl.read_stat->iter_idx = NULL;}
	if (pl.read_stat->iter_map != NULL) {free(pl.read_stat->iter_map); pl.read_stat->iter_map = NULL;}
	if (pl.read_stat->is_ove != NULL) {free(pl.read_stat->is_ove); pl.read_stat->is_ove = NULL;}
	if (pl.read_stat != NULL) {free(pl.read_stat); pl.read_stat = NULL;}

	if (pl.file_idx.offs != NULL) {free(pl.file_idx.offs); pl.file_idx.offs = NULL;}

	return 0;
}
