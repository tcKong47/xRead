#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <malloc.h>
#include <string.h>
#include <math.h>

#include "main.h"
#include "ktime.h"
#include "overlapping.h"
#include "bit_operation.h"

#define _DEBUG

void init_graph(Graph *Dp_graph, uint32_t max_v_num)
{
	Dp_graph->vnode = (VNode *)calloc(max_v_num, sizeof(VNode));
	if (Dp_graph->vnode == NULL)
	{
		fprintf(stderr, "[%s] calloc %ldGB Dp_graph->vnode memory error!\n", __func__, max_v_num * sizeof(VNode) / 1024 / 1024 / 1024);
		exit(1);
	}
}

void free_graph(Graph *Dp_graph)
{
	if (Dp_graph->vnode != NULL) {free(Dp_graph->vnode);Dp_graph->vnode = NULL;}
}

static void find_longest_path(Graph *graph, uint32_t vi, PATH_t *dist_path, uint32_t vertex_num)
{
	float temp = 0;
    float penalty;
	float current_dist = 0;
	int32_t pre_node = -1;
    int weight, dir;
	uint32_t adjvex;

    //travel the processer of vertex vi
    pANode arcnode;
	int i;
	for (i = 0; i < graph->vnode[vi].anode_num; ++i)
	{
		arcnode = &graph->vnode[vi].preedge[i];
		weight = arcnode->weight;
		penalty = arcnode->penalty;
		dir = arcnode->dir;
		adjvex = arcnode->adjvex;
		temp = dist_path[adjvex].dist + weight - penalty;

		if (current_dist <= temp && dist_path[adjvex].dir == dir) //only walk to the same direction
        {
            current_dist = temp;
            pre_node = adjvex;
        }
	}

	dist_path[vi].dist = current_dist; // change the distance because there are overlaps between mems
	dist_path[vi].pre_node = pre_node; //the front node

	if (vi >= vertex_num) printf("[%s] use unlegal memory, %d-%d, something wrong with it...", __func__, vi, vertex_num);
	return;
}

static int sparse_dynamic_programming(Graph *graph, uint32_t vertex_num, MR_t *vertex_mr, PATH_t *dist_path, uint8_t *out_degree)
{
	uint32_t i;

	for (i = 0; i < vertex_num; ++i)
	{
		if (graph->vnode[i].anode_num > 0)
		{
			find_longest_path(graph, i, dist_path, vertex_num);
		}
	}

	#ifdef _DEBUG
	int j;
	printf("show the path");
	for (i = 0; i < vertex_num; ++i)
	{
		j = i;
		if ((out_degree[i] == 0))
		{
			printf("\nvertex: %d\tdist = %f\n", i, dist_path[i].dist);
			printf("path: ");
			printf("%d(%d,%ld,%ld,%ld,%ld)->", i, vertex_mr[i].t_id, vertex_mr[i].ts, vertex_mr[i].te, vertex_mr[i].qs, vertex_mr[i].qe);

			j = dist_path[j].pre_node;
			while (j != -1)
			{
				printf("%d(%d,%ld,%ld,%ld,%ld)->", j, vertex_mr[j].t_id, vertex_mr[j].ts, vertex_mr[j].te, vertex_mr[j].qs, vertex_mr[j].qe);
				j = dist_path[j].pre_node;
			}
		}
	}
	printf("\nend\n");
	#endif

	return 0;
}

uint32_t create_graph(uint32_t cur_s_s, MR_t *vertex_mr, uint32_t vertex_num, uint32_t *max_index, PATH_t *dist_path, Graph *graph)
{
	int32_t i, j, weight, gap;
	int32_t ove_q, ove_tt;
	uint32_t non_ioslated_point = 0;
	uint32_t thre_num;
	int search_step = (vertex_num < cur_s_s)? vertex_num : cur_s_s, t_id, index_n;
	uint32_t *max_dis = NULL;
	uint8_t *out_degree = NULL;

#ifdef _DEBUG
	int vi = 0;
	printf("qid\thid\tcov\tqstrand\tqs\tqe\ttid\ttstrand\tts\tte\tlen\n");
	for (vi = 0; vi < vertex_num; ++vi)
	{
		printf("%d\t%d\t%d\t%ld\t%ld\t%d\t%d\t%ld\t%ld\n", vi, vertex_mr[vi].cov, vertex_mr[vi].qstr, vertex_mr[vi].qs, vertex_mr[vi].qe, vertex_mr[vi].t_id, vertex_mr[vi].tstr, vertex_mr[vi].ts, vertex_mr[vi].te);
	}
	printf("\n");
#endif

	// re inital
	for (i = 0; i < vertex_num; ++i)
	{
		graph->vnode[i].anode_num = 0;
		graph->vnode[i].preedge = (ANode*)calloc(search_step, sizeof(ANode));
		if (graph->vnode[i].preedge == NULL)
		{
			fprintf(stderr, "[%s] calloc %ldGB graph->vnode[i].preedge memory error!\n", __func__, search_step * sizeof(ANode) / 1024 / 1024 / 1024);
			exit(1);
		}
		dist_path[i].dir = vertex_mr[i].qstr == vertex_mr[i].tstr ? 0 : 1;
		dist_path[i].dist = vertex_mr[i].cov;
		dist_path[i].pre_node = -1;
	}

	max_dis = (uint32_t *)calloc(vertex_num, sizeof(uint32_t));
	if (max_dis == NULL)
	{
		fprintf(stderr, "[%s] calloc %ldGB max_dis memory error!\n", __func__, vertex_num * sizeof(uint32_t) / 1024 / 1024 / 1024);
		exit(1);
	}

	out_degree = (uint8_t *)calloc(vertex_num, sizeof(uint8_t));
	if (out_degree == NULL)
	{
		fprintf(stderr, "[%s] calloc %ldGB out_degree memory error!\n", __func__, vertex_num * sizeof(uint8_t) / 1024 / 1024 / 1024);
		exit(1);
	}

	uint32_t edge_num;
	for (i = 0; i < vertex_num - 1; ++i)
	{
		edge_num = 0;
		thre_num = (vertex_num < (i + search_step))? vertex_num : (i + search_step);
		for (j = i + 1; j < thre_num; ++j)
		{
			if (vertex_mr[i].t_id != vertex_mr[j].t_id)	break;
			if ((vertex_mr[i].tstr ^ vertex_mr[j].tstr) != (vertex_mr[i].qstr ^ vertex_mr[j].qstr))	continue;

			ove_q = vertex_mr[j].qs - vertex_mr[i].qe - 1;
			ove_tt = vertex_mr[j].tstr == vertex_mr[j].qstr ? (int32_t)(vertex_mr[j].ts - vertex_mr[i].te - 1) : (int32_t)(vertex_mr[i].ts - vertex_mr[j].te - 1);
			gap = (int32_t)abs(ove_q - ove_tt);

			if (ove_q > -k && ove_tt > -k && (gap <= abs(ove_q)*0.3 || gap < k))
            {
            	graph->vnode[j].preedge[graph->vnode[j].anode_num].adjvex = i;
				weight = ove_q >= 0 ? vertex_mr[j].cov: vertex_mr[j].cov + ove_q;
				graph->vnode[j].preedge[graph->vnode[j].anode_num].weight = weight;
				graph->vnode[j].preedge[graph->vnode[j].anode_num].penalty = gap / (float)weight;
				graph->vnode[j].preedge[graph->vnode[j].anode_num++].dir = vertex_mr[j].tstr == vertex_mr[j].qstr ? 0 : 1;
				out_degree[i] = 1;

				edge_num++;
        		non_ioslated_point++;
            }
			if (edge_num > 10)	break;
		}
	}

	sparse_dynamic_programming(graph, vertex_num, vertex_mr, dist_path, out_degree);

	index_n = 0;
	t_id = vertex_mr[0].t_id;
	for (i = 0; i < vertex_num; ++i)
	{
		if (vertex_mr[i].t_id != t_id)
		{
			if (dist_path[max_index[index_n]].pre_node != -1 || max_dis[index_n] > k)	index_n++;
			t_id = vertex_mr[i].t_id;
			max_index[index_n] = i;
		}
		if (max_dis[index_n] < dist_path[i].dist)
		{
			max_dis[index_n] = dist_path[i].dist;
			max_index[index_n] = i;
		}
	}
	if (dist_path[max_index[index_n]].pre_node != -1 || max_dis[index_n] > k) index_n++; //deal with the last node

	for (i = 0; i < vertex_num; ++i)
	{
		if (graph->vnode[i].preedge != NULL) {free(graph->vnode[i].preedge);graph->vnode[i].preedge = NULL;}
	}
	if (max_dis != NULL) {free(max_dis);max_dis = NULL;}
	if (out_degree != NULL) {free(out_degree);out_degree = NULL;}

	return index_n;
}