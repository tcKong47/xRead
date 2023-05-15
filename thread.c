#include <pthread.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

/************
 * kt_for() *
 ************/

#include "overlapping.h"

struct kt_for_t;

typedef struct {
	struct kt_for_t *t;
	uint64_t i;
} ktf_worker_t;

typedef struct kt_for_t {
	int n_threads;
	int64_t n;
    void (*func)(void *, int64_t, int);
    void *data;
    ktf_worker_t *workers;
} kt_for_t;

static inline int64_t steal_work(kt_for_t *t)
{
	int i, min_i = -1;
	int64_t k, min = UINT64_MAX;
	for (i = 0; i < t->n_threads; ++i)
	{
        if (min > t->workers[i].i)
        {
            min = t->workers[i].i;
            min_i = i;
        }
    }
	k = __sync_fetch_and_add(&t->workers[min_i].i, t->n_threads);
    k = k >= t->n ? -1 : k;
	return k;
}

static void *ktf_worker(void *data)
{
	ktf_worker_t *w = (ktf_worker_t*)data;
	int64_t i;
	for (;;) {
		i = __sync_fetch_and_add(&w->i, w->t->n_threads);
		if (i >= w->t->n) break;
		w->t->func(w->t->data, i, w - w->t->workers);
	}
	while ((i = steal_work(w->t)) >= 0)
	{
        w->t->func(w->t->data, i, w - w->t->workers);
    }
	// pthread_exit(0);
    return ((void*)0);
}

void kt_for(int n_threads, void (*func)(void *, int64_t, int), void *data, int64_t n)
{
    int i;
	kt_for_t aux;
	pthread_t *pid = NULL;
    pthread_attr_t attr;

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    aux.func = func;
    aux.data = data;
    aux.n_threads = n_threads;
    aux.n = n;
	aux.workers = (ktf_worker_t*)calloc(n_threads , sizeof(ktf_worker_t));
	for (i = 0; i < n_threads; ++i)
	{
        aux.workers[i].t = &aux;
        aux.workers[i].i = i;
    }

    pid = (pthread_t *)calloc(n_threads , sizeof(pthread_t));
    for (i = 0; i < n_threads; ++i) pthread_create(&pid[i], &attr, ktf_worker, &aux.workers[i]);
    for (i = 0; i < n_threads; ++i) pthread_join(pid[i], 0);

    if (aux.workers != NULL) {free(aux.workers); aux.workers = NULL;}
    if (pid != NULL) {free(pid); pid = NULL;}
}

/*****************
 * kt_pipeline() *
 *****************/

struct ktp_t;

typedef struct
{
    struct ktp_t *pl;
    int64_t index;
    int step;
    void *data;
} ktp_worker_t;

typedef struct ktp_t
{
    void *shared;
    void *(*func)(void *, int, void *);
    int64_t index;
    int n_workers, n_steps;
    ktp_worker_t *workers;
    pthread_mutex_t mutex;
    pthread_cond_t cv;
} ktp_t;

static void *ktp_worker(void *data)
{
    ktp_worker_t *w = (ktp_worker_t *)data;
    ktp_t *p = w->pl;
    while(w->step < p->n_steps)
    {
        pthread_mutex_lock(&p->mutex);
        for(;;)
        {
            int i;
            for (i = 0; i < p->n_workers; i++)
            {
                if (w == &p->workers[i]) continue;
                if (p->workers[i].step <= w->step && p->workers[i].index < w->index) break;
            }
            if (i == p->n_workers) break;
            // fprintf(stderr, "[%s]thread %ld, step %d, waiting\n", __func__, w->index, w->step);
            pthread_cond_wait(&p->cv, &p->mutex);
        }
        pthread_mutex_unlock(&p->mutex);

        // fprintf(stderr, "[%s]thread %ld, step %d\n", __func__, w->index, w->step);
        w->data = p->func(p->shared, w->step, w->step? w->data : 0);

        pthread_mutex_lock(&p->mutex);
        w->step = w->step == p->n_steps - 1 || w->data ? (w->step + 1) % p->n_steps : p->n_steps;
        if (w->step == 0) w->index = p->index++;
        // fprintf(stderr, "[%s]thread %ld, step %d, update\n", __func__, w->index, w->step);
        pthread_cond_broadcast(&p->cv);
        pthread_mutex_unlock(&p->mutex);
    }
    // pthread_exit(0);
    return ((void*)0);
}

void kt_pipeline(int n_threads, void *(*func)(void *, int, void *), void *shared_data, int n_steps)
{
    int i;
    ktp_t aux;
    pthread_t *pid = NULL;
    pthread_attr_t attr;

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    aux.n_workers = n_threads;
    aux.n_steps = n_steps;
    aux.func = func;
    aux.shared = shared_data;
    aux.index = 0;
    pthread_mutex_init(&aux.mutex, 0);
    pthread_cond_init(&aux.cv, 0);
    aux.workers = (ktp_worker_t *)calloc(n_threads, sizeof(ktp_worker_t));
    for(i = 0; i < n_threads; ++i)
    {
        aux.workers[i].step = 0;
        aux.workers[i].pl = &aux;
        aux.workers[i].data = 0;
        aux.workers[i].index = aux.index++;
    }

    pid = (pthread_t *)calloc(n_threads , sizeof(pthread_t));
    for(i = 0; i < n_threads; ++i)	pthread_create(&pid[i], &attr, ktp_worker, &aux.workers[i]);
    for(i = 0; i < n_threads; ++i)	pthread_join(pid[i], 0);

    pthread_mutex_destroy(&aux.mutex);
    pthread_cond_destroy(&aux.cv);

    if (aux.workers != NULL) {free(aux.workers); aux.workers = NULL;}
    if (pid != NULL) {free(pid); pid = NULL;}
}