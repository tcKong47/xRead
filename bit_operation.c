#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "bit_operation.h"
#include "overlapping.h"
// #include "assemblygraph.h"
// #include "linking.h"

uint8_t seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

int compare_rlen(const void *a, const void *b)
{
	len_t *l1 = (len_t *)a;
	len_t *l2 = (len_t *)b;

	if (l1->rlen < l2->rlen)
		return 1;
	if (l1->rlen > l2->rlen)
		return -1;
	else
		return 0;
}

int compare_mi(const void *a, const void *b)
{
	MINIMIZER_t *h1 = (MINIMIZER_t *)a;
	MINIMIZER_t *h2 = (MINIMIZER_t *)b;

	if (h1->x > h2->x)
		return 1;
	if (h1->x < h2->x)
		return -1;
	else
	{
		if (h1->y > h2->y)
			return 1;
		if (h1->y < h2->y)
			return -1;
		return 0;
	}
}

int compare_rep(const void *a, const void *b)
{
	uint32_t *h1 = (uint32_t *)a;
	uint32_t *h2 = (uint32_t *)b;

	if (*h1 < *h2)
		return 1;
	if (*h1 > *h2)
		return -1;
	else
		return 0;
}

int compare_tid(const void *a, const void *b)
{
	MR_t *v1 = (MR_t *)a;
	MR_t *v2 = (MR_t *)b;

	if (v1->t_id > v2->t_id)
		return 1;
	if (v1->t_id < v2->t_id)
		return -1;
	else
	{
		if (v1->qs > v2->qs)
			return 1;
		if (v1->qs < v2->qs)
			return -1;
		else
		{
			if (v1->ts > v2->ts)
				return 1;
			if (v1->ts < v2->ts)
				return -1;
			else
				return 0;
		}
	}
}

int compare_qs(const void *a, const void *b)
{
	MR_t *v1 = (MR_t *)a;
	MR_t *v2 = (MR_t *)b;

	if (v1->qs > v2->qs)
		return 1;
	if (v1->qs < v2->qs)
		return -1;
	else
	{
		if (v1->ts > v2->ts)
			return 1;
		if (v1->ts < v2->ts)
			return -1;
		else
			return 0;
	}
}

int compare_ove_len(const void *a, const void *b)
{
	OVE_t *h1 = (OVE_t *)a;
	OVE_t *h2 = (OVE_t *)b;

	if (h1->mln < h2->mln)
		return 1;
	if (h1->mln > h2->mln)
		return -1;
	else
		return 0;
}

int compare_ove_score(const void *a, const void *b)
{
	OVE_t *h1 = (OVE_t *)a;
	OVE_t *h2 = (OVE_t *)b;

	if (h1->score < h2->score)
		return 1;
	if (h1->score > h2->score)
		return -1;
	else
		return 0;

	// if ((double)h1->mbp/h1->mln < (double)h2->mbp/h2->mln)
	// 	return 1;
	// if ((double)h1->mbp/h1->mln > (double)h2->mbp/h2->mln)
	// 	return -1;
	// else
	// 	return 0;
}

// int compare_ave_cov(const void *a, const void *b)
// {
// 	read_cov_t *h1 = (read_cov_t *)a;
// 	read_cov_t *h2 = (read_cov_t *)b;

// 	if (h1->ave_cov > h2->ave_cov)
// 		return 1;
// 	if (h1->ave_cov < h2->ave_cov)
// 		return -1;
// 	else
// 		return 0;
// }

int compare_uint32_t(const void *a, const void *b)
{
	uint32_t *h1 = (uint32_t *)a;
	uint32_t *h2 = (uint32_t *)b;

	if (*h1 < *h2)
		return 1;
	if (*h1 > *h2)
		return -1;
	else
		return 0;
}

int compare_double(const void *a, const void *b)
{
	double *v1 = (double *)a;
	double *v2 = (double *)b;

	if (*v1 > *v2)
		return 1;
	if (*v1 < *v2)
		return -1;
	else
		return 0;
}

uint64_t hash64(uint64_t key, uint64_t mask)
{
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

