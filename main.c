#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include "ktime.h"
#include "main.h"
#include "overlapping.h"

//global variable
float x_read;
int k;
int w;
int l;
float r_n;
int s_s;
int m_b;
int min_ove;
int max_hang;
int top_n;
float cov_ratio;
int part_size;
int X_read;
int split_len;
int iter;
int trans_iter;
int trans_file;

int thread_n;
char *output_file;
int waitingLen;
uint64_t batch_size_base;
int memory;
double realtime0;

static void init_map_param(param_map *opt)
{
	opt->x_read = 3;
	opt->k = 15;
	opt->w = 5;
	opt->l = 11;
	opt->r_n = 0.0005;
	opt->s_s = 10;

	opt->m_b = 100;
	opt->min_ove = 500;
	opt->max_hang = 2000;
	opt->top_n = 2;

	opt->cov_ratio = 0.5;
	opt->part_size = 3;
	opt->X_read = 10;
	opt->split_len = 1000;
	opt->iter = 10;
	
	opt->trans_iter = 0;
	opt->trans_file = 0;

	opt->read_type = 1;
	opt->thread_n = 8;
	opt->memory = 16;
	opt->batch_size = opt->memory / (19 + 2.0/opt->thread_n) * 1024 *1024 * 1024;
	// 16 bytes for storing minimizers, 3 bytes for stroing map read bases
	// 1 bytes for storing other structs
}

static int help_usage(param_map *opt)
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:	xRead\n");

	fprintf(stderr, "Usage:		xRead <reads.fasta/fastq> [options]\n");
	fprintf(stderr, "\nProgram options:   \n");
    fprintf(stderr, "    -h --help                           Print help menu.\n");
    fprintf(stderr, "    -f --output-file		    [STR]    The path and name prefix of output file, if not set, using the stdout.\n");
	fprintf(stderr, "    -M --memory                [INT]    Maximum allowed memory. [%d]\n", opt->memory);
	fprintf(stderr, "    -t --thread-n              [INT]    Number of used threads. [%d]\n", opt->thread_n);
	fprintf(stderr, "    -p --read-type             [INT]    Specifiy the type of reads and set multiple paramters unless overriden. [%d]\n", opt->read_type);
	fprintf(stderr, "                                        	[1]: For eads with error rates ~12%%\n");
	fprintf(stderr, "                                        	[2]: For reads with error rates ~1%%\n");

	fprintf(stderr, "\nAlgorithm options:\n");
	fprintf(stderr, "    -x --x-longest-read        [FLOAT]  The longest x%% reads indexed for the first iteration. [%f]\n", opt->x_read);
	fprintf(stderr, "    -X --x-random-read         [FLOAT]  The longest X%% reads indexed for the second and later iteration. [%d]\n", opt->X_read);
	fprintf(stderr, "    -k --k-mer                 [INT]    The length of k-mer for sequences sketching. [%d]\n", opt->k);
	fprintf(stderr, "    -w --window-size           [INT]    Window size for sequences sketching. [%d]\n", opt->w);
	fprintf(stderr, "    -l --l-mer                 [INT]    The length of l-mer of the auxiliary index hash table. [%d]\n", opt->l);
	fprintf(stderr, "    -r --repeat-n              [FLOAT]  The proportion for filtering the most repetitive minimizer hits. [%.3f]\n", opt->r_n);
	fprintf(stderr, "    -e --search-step           [INT] 	 The number of search step of SDP graph construction. [%d]\n", opt->s_s);

	fprintf(stderr, "    -n --top-n                 [INT]    The number of overlaps each read retained for each iteration. [%d]\n", opt->top_n);
	fprintf(stderr, "    -b --matching-bases        [INT]    Minimum required matching bases for a single overlap between two reads. [%d]\n", opt->m_b);
	fprintf(stderr, "    -m --min-ove               [INT]    Minimum required length of overlap. [%d]\n", opt->min_ove);
	fprintf(stderr, "    -a --max-hang              [INT]    Maximum allowed total length of left and right overhangs. [%d]\n", opt->max_hang);

	fprintf(stderr, "    -L --split-len             [INT]    The length of split blocks for read coverage estimation. [%d]\n", opt->split_len);
	fprintf(stderr, "    -S --read-part-size        [INT]    Estimating the average coverage of every S consecutive blocks to mark the less-covered read portions. [%d]\n", opt->part_size);
	fprintf(stderr, "    -c --cov-ratio             [FLOAT]  The proportion of less-covered reads for next iteration. [%.3f]\n", opt->cov_ratio);
	fprintf(stderr, "    -I --iter-times            [INT]    Maximum allowed number of iterations. [%d]\n", opt->iter);

	fprintf(stderr, "    -R --trans-iter            [INT]    Enabling the calculating of transitive overlaps for R iterations. [%d]\n", opt->iter);
	fprintf(stderr, "    -N --trans-file            	     Skipping the overlapping discovery, only perfomred R transitive iterations to generate comprehensive graph.\n");
	
	fprintf(stderr, "\nProgram:	xRead\n");
	fprintf(stderr, "Usage:		xRead <overlaps.paf> -N [options]\n");
	fprintf(stderr, "    -R --trans-iter            [INT]    Enabling the calculating of transitive overlaps for R iterations. [%d]\n", opt->iter);
	
	return 1;
}

static const char *short_option = "f:p:t:x:k:w:l:r:e:b:m:a:n:c:S:X:L:I:R:NM:h";

static struct option long_option[] = {
	{"output-file", required_argument, NULL, 'f'},
	{"read-type", required_argument, NULL, 'p'},
	{"thread-n", required_argument, NULL, 't'},

	{"x-longest-read", required_argument, NULL, 'x'},
	{"k-mer", required_argument, NULL, 'k'},
	{"window-size", required_argument, NULL, 'w'},
	{"l-mer", required_argument, NULL, 'l'},
	{"repeat-n", required_argument, NULL, 'r'},
	{"search-step", required_argument, NULL, 'e'},

	{"matching-bases", required_argument, NULL, 'b'},
	{"min-ove", required_argument, NULL, 'm'},
	{"max-hang", required_argument, NULL, 'a'},
	{"top-n", required_argument, NULL, 'n'},

	{"cov-ratio", required_argument, NULL, 'c'},
	{"read-part-size", required_argument, NULL, 'S'},
	{"x-random-read", required_argument, NULL, 'X'},
	{"split-len", required_argument, NULL, 'L'},
	{"iter-times", required_argument, NULL, 'I'},

	{"trans-iter", required_argument, NULL, 'R'},
	{"trans-file", no_argument, NULL, 'N'},

	{"memory", required_argument, NULL, 'M'},
	{"help", no_argument, NULL, 'h'},
	{0, 0, 0, 0}};

int main(int argc, char *argv[])
{
	realtime0 = realtime();
	int c;
	float error, t, q, els;
	param_map *opt = (param_map* )calloc(1, sizeof(param_map));
	
	init_map_param(opt);

	while((c = getopt_long(argc, argv, short_option, long_option, NULL)) != -1)
	{
		switch(c)
		{
			case 'f': opt->output_file = strdup(optarg); break;
			case 'p': opt->read_type = atoi(optarg); break;
			case 't': opt->thread_n = atoi(optarg); break;

			case 'x': opt->x_read = atof(optarg); break;
			case 'k': opt->k = atoi(optarg); break;
			case 'w': opt->w = atoi(optarg); break;
			case 'l': opt->l = atoi(optarg); break;
			case 'r': opt->r_n = atof(optarg); break;
			case 'e': opt->s_s = atoi(optarg); break;

			case 'b': opt->m_b = atoi(optarg); break;
			case 'm': opt->min_ove = atoi(optarg); break;
			case 'a': opt->max_hang = atoi(optarg);break;
			case 'n': opt->top_n = atoi(optarg); break;

			case 'c': opt->cov_ratio = atof(optarg);break;
			case 'S': opt->part_size = atoi(optarg);break;
			case 'X': opt->X_read = atoi(optarg);break;
			case 'L': opt->split_len = atoi(optarg);break;
			case 'I': opt->iter = atoi(optarg);break;

			case 'R': opt->trans_iter = atoi(optarg);break;
			case 'N': opt->trans_file = 1; break;

			case 'M': opt->memory = atoi(optarg);break;
			case 'h': return help_usage(opt); break;
			default: return help_usage(opt); break;
		}
	} 

	if (argc - optind < 1)
		return help_usage(opt);

	if (opt->read_type > 2)
	{
		fprintf(stderr, "Input warning: read type should be 1 or 2, 1 default\n");
		opt->read_type = 1;
	}
	if (opt->thread_n < 1 || opt->thread_n > 128)
	{
		fprintf(stderr, "Input error: -t cannot be less than 1 or more than 128\n");
		exit(1);
	}
	if (opt->x_read <= 0 || opt->x_read > 100)
	{
		fprintf(stderr, "Input error: -x cannot be less than 0 or more than 100\n");
		exit(1);
	}
	if (opt->k < 10 || opt->k > 22)
	{
		fprintf(stderr, "Input error: -k cannot be less than 10 or more than 20\n");
		exit(1);
	}
	if (opt->w < 1 || opt->w > 64)
	{
		fprintf(stderr, "Input error: -w cannot be less than 1 or more than 10\n");
		exit(1);
	}
	if (opt->l < 6 || opt->l > 14)
	{
		fprintf(stderr, "Input error: -l cannot be less than 6 or more than 14\n");
		exit(1);
	}
	if (opt->l >= opt->k)
	{
		fprintf(stderr, "Input error: -l cannot be more than or equal to -k\n");
		exit(1);
	}
	if (opt->s_s < 1)
	{
		fprintf(stderr, "Input error: -e cannot be less than 1\n");
		exit(1);
	}
	if (opt->m_b < opt->k)
	{
		fprintf(stderr, "Input error: -b cannot be less than -k\n");
		exit(1);
	}
	if (opt->min_ove < opt->k)
	{
		fprintf(stderr, "Input error: -m cannot be less than -k\n");
		exit(1);
	}
	if (opt->max_hang < opt->k)
	{
		fprintf(stderr, "Input error: -a cannot be less than -k\n");
		exit(1);
	}
	if (opt->top_n < 1)
	{
		fprintf(stderr, "Input error: -n cannot be less than 1\n");
		exit(1);
	}
	if (opt->cov_ratio < 0.09 || opt->cov_ratio > 0.81)
	{
		fprintf(stderr, "Input error: -c cannot be less than 0.1 or more than 0.8\n");
		exit(1);
	}
	if (opt->part_size < 1 || opt->part_size > 10)
	{
		fprintf(stderr, "Input error: -S cannot be less than 1 or more than 10\n");
		exit(1);
	}
	if (opt->X_read < 1 || opt->X_read > 100)
	{
		fprintf(stderr, "Input error: -X cannot be less than 1 or more than 100\n");
		exit(1);
	}
	if (opt->split_len < 500 || opt->split_len > 2000)
	{
		fprintf(stderr, "Input error: -L cannot be less than 500 or more than 2000\n");
		exit(1);
	}
	if (opt->iter < 1 || opt->iter > 16)
	{
		fprintf(stderr, "Input error: -I cannot be less than 1 or more than 16\n");
		exit(1);
	}
	if (opt->memory < 2 || opt->memory > 128)
	{
		fprintf(stderr, "Input error: -M cannot be less than 4 or more than 128\n");
		exit(1);
	}

	char *read_fastq;
	read_fastq = strdup(argv[optind]);
	char *index_fastq;
	index_fastq = strdup(argv[optind]);

	els = 0.05;
	error = 0.2;
	if (opt->read_type == 1) //ont 88%
		error = 0.12;
	else if (opt->read_type == 2) //ont 99%
		error = 0.01;

	t = log10(els)/log10(1 - pow(1 - error, opt->k));
	q = 1/error - opt->k * pow(1 - error, opt->k)/(1 - pow(1 - error, opt->k));
	waitingLen = (int)(t * q);
	opt->batch_size = opt->memory / (19 + 2.0/opt->thread_n) * 1024 *1024 * 1024;
	fprintf(stderr, "[Initial] Configuring memory limitaion as %dGB, and batch size as %ld\n", opt->memory, opt->batch_size);

	if (opt->batch_size < 100000000 || opt->batch_size > 5000000000)
	{
		fprintf(stderr, "Input error: batch size cannot be less than 100000000 or more than 4000000000\n");
		exit(1);
	}

	x_read = opt->x_read;
	k = opt->k;
	w = opt->w;
	l = opt->l;
	r_n = opt->r_n;
	s_s = opt->s_s;

	m_b = opt->m_b;
	min_ove = opt->min_ove;
	max_hang = opt->max_hang;
	top_n = opt->top_n;

	cov_ratio = opt->cov_ratio;
	part_size = opt->part_size;
	X_read = opt->X_read;
	split_len = opt->split_len;
	iter = opt->iter;
	trans_iter = opt->trans_iter;
	trans_file = opt->trans_file;
	memory = opt->memory;
	batch_size_base = opt->batch_size;
	thread_n = opt->thread_n;
	output_file = opt->output_file;

	if (trans_iter > 0 && trans_file == 1)
	{
		fprintf(stderr, "[Initial] Start generateing transitive read overlaps...\n");
		generate_transitive_overlaps_file(read_fastq);
	}
	else
	{
		fprintf(stderr, "[Initial] Start discovering read overlaps...\n");
		finding_overlapping(index_fastq, read_fastq, output_file);
	}
	if (opt != NULL)	free(opt);
	if (read_fastq != NULL) {free(read_fastq); read_fastq = NULL;}
	if (index_fastq != NULL) {free(index_fastq); index_fastq = NULL;}

	fprintf(stderr, "[Result]\tReal time:%.3f sec\t%.3fh\tCPU:%.3f sec\tMemory peak:%.3f GB\n", realtime() - realtime0, (realtime() - realtime0) / 3600.0, cputime(), peak_memory() / 1024.0 / 1024.0);

	return 0;
}
