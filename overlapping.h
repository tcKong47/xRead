/*************************************************************************
	> File Name: overlapping.h
	> Author: Tangchao Kong
	> Mail: 
 ************************************************************************/

#ifndef OVERLAPPING_H
#define OVERLAPPING_H

#include "main.h"

typedef struct ArcNode
{
	uint32_t adjvex; //id of precursor vertex
	uint32_t weight; //weight of edge
	float penalty;	 //penalty of edge
	int dir;		 //direction of edge
} ANode, *pANode;

typedef struct VertexNode
{
	uint32_t anode_num; //number of adjacent node
	ANode *preedge;		//edge link between current node and precursor
} VNode;

typedef struct dpGraph
{
	VNode *vnode; //list of vertex
} Graph;

typedef struct
{
	uint8_t tstr, qstr;
	uint32_t t_id;
	int32_t cov;
	uint64_t ts, te;
	uint64_t qs, qe;
}MR_t;

typedef struct PATH
{
	int dir;
	float dist;
	int32_t pre_node;
}PATH_t;

typedef struct
{
	uint32_t idx_n;
	int32_t *idx;
} path_idx_t;

typedef struct
{
	uint32_t qid, tid;
	uint32_t ql, tl;
	uint32_t qs, ts;
	uint32_t qe, te;
	uint8_t rev;
	uint32_t mbp;
	uint32_t mln;
	float score;
}OVE_t;

typedef struct
{
	uint32_t n, m;
	OVE_t *ove;
}OVE_C_t;

typedef struct
{
	MR_t *vertex_mr;
	uint32_t v_n, v_m;
	uint32_t max_v_num;
	OVE_C_t *ove_cl;
	uint32_t *is_ove;
} thread_buf_t;

typedef struct
{
	uint32_t rid; // rid of last read
	uint32_t sta_pos, end_pos; // sta pos in seed read
	uint32_t ts, te; // sta pos in read[rid]
	uint8_t rev;
} cov_data_t;

typedef struct
{
	uint16_t *cov;
} read_cov_t;

typedef struct
{
	uint32_t *rlen;
	uint32_t *is_idx;
	uint32_t *iter_idx;
	uint32_t iter_idx_n;
	uint32_t *iter_map;
	uint32_t iter_map_n;
	uint32_t *is_ove;
}read_stat_t;

typedef struct
{
	uint32_t rlen;
	uint32_t nlen;
} len_t;

typedef struct
{
	uint32_t last_rid;
	uint32_t upper_v;
	uint32_t add_n, add_m;
    uint32_t *add_p;
} file_add_p_t;

uint32_t generate_transitive_overlaps_file(const char *trans_file);
uint32_t finding_overlapping(const char *read_fastq, const char *index_fastq, const char *temp_file_perfix);

#endif