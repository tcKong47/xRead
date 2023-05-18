#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bseq.h"
#include "kseq.h"
#include "overlapping.h"
#include "bit_operation.h"

KSEQ_INIT(gzFile, gzread)

struct bseq_file_s {
	gzFile fp;
	kseq_t *ks;
};

bseq_file_t *bseq_open(const char *fn)
{
	bseq_file_t *fp;
	gzFile f;
	f = fn && strcmp(fn, "-")? gzopen(fn, "rb") : gzdopen(fileno(stdin), "rb");
	if (f == 0) return 0;
	fp = (bseq_file_t*)calloc(1, sizeof(bseq_file_t));
	fp->fp = f;
	fp->ks = kseq_init(fp->fp);
	return fp;
}

void bseq_seek_cur(bseq_file_t *fp, long int offset)
{
    offset = kseq_offset(fp->ks, offset);
    gzseek(fp->fp, offset, SEEK_CUR);
}

void bseq_close(bseq_file_t *fp)
{
	kseq_destroy(fp->ks);
	gzclose(fp->fp);
	free(fp);
}

len_t *bseq_read_l(bseq_file_t *fp, uint64_t batch, uint32_t *_n)
{
    len_t *len = NULL;
    uint32_t seqii = 0, m = 0;
    int64_t kr1 = 1;
    kseq_t *seq = fp->ks;
    for (seqii = 0; (seqii < batch) && ((kr1 = kseq_read(seq)) > 0); seqii++)
    {
        if (seqii >= m)
        {
            m = m == 0 ? 256 : m << 1;
            len = (len_t *)realloc(len, m * sizeof(len_t));
            if (len == NULL)
            {
                fprintf(stderr, "[%s] calloc %ldGB len memory error!\n", __func__, m * sizeof(len_t) / 1024 /1024 / 1024);
                exit(1);
            }
        }
        len[seqii].rlen = seq->seq.l;
        len[seqii].nlen = seq->name.l + seq->comment.l + seq->qual.l + seq->skip_base;
    }
    *_n = seqii;
    return len;
}

READ_t *bseq_read(bseq_file_t *fp, uint64_t batch, uint32_t *_n)
{
    READ_t *read_info = NULL;
    uint32_t idx = 0, m = 0;
    uint64_t base = 0;
    int64_t kr1 = 1;
    kseq_t *seq1 = fp->ks;

    for (; (base < batch) && ((kr1 = kseq_read(seq1)) > 0); )
    {
        if (idx >= m)
        {
            m = m == 0 ? 256 : m << 1;
            read_info = (READ_t *)realloc(read_info, m * sizeof(READ_t));
            if (read_info == NULL)
            {
                fprintf(stderr, "[%s] calloc %ldGB read_info memory error!\n", __func__, m * sizeof(READ_t) / 1024 /1024 / 1024);
                exit(1);
            }
        }
        read_info[idx].read_name = strdup(seq1->name.s);
        read_info[idx].read_length = seq1->seq.l;
        read_info[idx].read_seq = strdup(seq1->seq.s);
        base += read_info[idx++].read_length;
    }
    *_n = idx;
    return read_info;
}

int bseq_rewrite_read(bseq_file_t *fp, FILE *fp_un_ove, uint32_t *symn, file_idx_t *file_idx, char file_type)
{
    char c = file_type == 'q' ? '@' : '>';
    int read_num = 0;
    uint32_t seqi = 0;
    int64_t kr1 = 1;
    kseq_t *seq1 = fp->ks;
    for (; ; seqi++)
    {
        if (symn[seqi] == 0)
        {
            kr1 = kseq_read(seq1);
            if (kr1 <= 0) break;
            
            if (seq1->comment.l > 0) fprintf(fp_un_ove, "%c%s %s\n", c,(seq1->name.s), (seq1->comment.s));
            else fprintf(fp_un_ove, "%c%s\n", c, (seq1->name.s));
            fprintf(fp_un_ove, "%s\n", (seq1->seq.s));
            if (seq1->qual.l > 0)
            {
                if (seq1->comment.l > 0) fprintf(fp_un_ove, "+%s %s\n", (seq1->name.s), (seq1->comment.s));
                else fprintf(fp_un_ove, "+%s\n", (seq1->name.s));
                fprintf(fp_un_ove, "%s\n", (seq1->qual.s));
            }
            read_num++;
        }
        else
        {
            if (seqi >= file_idx->n) break;
            bseq_seek_cur(fp, file_idx->offs[seqi]);
        }
    }
    return read_num;
}


// the rid(seqii) is the same of all kind of reads
READ_t *bseq_read_idx(bseq_file_t *fp, uint64_t batch, uint32_t *_n, uint32_t *_sta, uint32_t *_sta_idx, uint32_t *iter_idx, file_idx_t *file_idx)
{
    READ_t *read_info = NULL;
    uint32_t idx = 0, m = 0, seqi = *_sta, seqidx = *_sta_idx;
    uint64_t base = 0;
    int64_t kr1 = 1;
    kseq_t *seq1 = fp->ks;

    for (; base < batch; seqi++)
    {
        if (seqi == iter_idx[seqidx])
        {
            kr1 = kseq_read(seq1);
            if (kr1 <= 0) break;
            if (idx >= m)
            {
                m = m == 0 ? 256 : m << 1;
                read_info = (READ_t *)realloc(read_info, m * sizeof(READ_t));
                if (read_info == NULL)
                {
                    fprintf(stderr, "[%s] calloc %ldGB read_info memory error!\n", __func__, m * sizeof(READ_t) / 1024 /1024 / 1024);
                    exit(1);
                }
            }
            read_info[idx].rid = seqi;
            read_info[idx].read_name = strdup(seq1->name.s);
            read_info[idx].read_length = seq1->seq.l;
            read_info[idx].read_seq = strdup(seq1->seq.s);
            base += read_info[idx++].read_length;
            seqidx++;
        }
        else
        {
            if (seqi >= file_idx->n) break;
            bseq_seek_cur(fp, file_idx->offs[seqi]);
        }
    }
    *_n = idx;
    *_sta = seqi;
    *_sta_idx = seqidx;
    return read_info;
}

READ_t *bseq_read_map(bseq_file_t *fp, uint64_t batch, uint32_t *_n, uint32_t *_sta, uint32_t *_sta_idx, uint32_t *iter_map, file_idx_t *file_idx)
{
    READ_t *read_info = NULL;
    uint32_t idx = 0, m = 0, seqi = *_sta, seqidx = *_sta_idx;
    uint64_t base = 0;
    int64_t kr1 = 1;
    kseq_t *seq1 = fp->ks;

    for (; base < batch; seqi++)
    {
        if (seqi == iter_map[seqidx])
        {
            kr1 = kseq_read(seq1);
            if (kr1 <= 0) break;
            if (idx >= m)
            {
                m = m == 0 ? 32767 : m << 1;
                read_info = (READ_t *)realloc(read_info, m * sizeof(READ_t));
                if (read_info == NULL)
                {
                    fprintf(stderr, "[%s] calloc %ldGB read_info memory error!\n", __func__, m * sizeof(READ_t) / 1024 /1024 / 1024);
                    exit(1);
                }
            }
            read_info[idx].rid = seqi;
            read_info[idx].read_name = strdup(seq1->name.s);
            read_info[idx].read_length = seq1->seq.l;
            read_info[idx].read_seq = strdup(seq1->seq.s);
            base += read_info[idx++].read_length;
            seqidx++;
        }
        else
        {
            if (seqi >= file_idx->n) break;
            bseq_seek_cur(fp, file_idx->offs[seqi]);
        }
    }
    *_n = idx;
    *_sta = seqi;
    *_sta_idx = seqidx;
    return read_info;
}