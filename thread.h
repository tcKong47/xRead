#ifndef THREAD_H
#define THREAD_H

#include <stdint.h>

void kt_pipeline(int n_threads, void *(*func)(void *, int, void *), void *shared_data, int n_steps);
void kt_for(int n_threads, void (*func)(void *, int64_t, int), void *data, int64_t n);


#endif