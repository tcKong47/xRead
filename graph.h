#ifndef GRAPH_H
#define GRAPH_H

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

#include "overlapping.h"

void init_graph(Graph *graph, uint32_t max_v_num);
void free_graph(Graph *graph);
uint32_t create_graph(uint32_t cur_s_s, MR_t *vertex_mr, uint32_t vertex_num, uint32_t *max_index, PATH_t *dist_path, Graph *graph);

// int create_graph(step_query_t *step, int64_t i, uint32_t pid);

#endif