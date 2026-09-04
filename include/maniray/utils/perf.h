#ifndef _MR_PERF_H
#define _MR_PERF_H

#include <time.h>

#ifndef MR_DISABLE_TIMER

#define MR_INIT_TIMER(start, end) struct timespec start, end;
#define MR_START_TIMER(start, end) clock_gettime(CLOCK_MONOTONIC, &start);
#define MR_STOP_TIMER(start, end, str) \
    do { \
        clock_gettime(CLOCK_MONOTONIC, &end); \
        long long __elapsed_us = (end.tv_sec - start.tv_sec) * 1000000LL + \
                            (end.tv_nsec - start.tv_nsec) / 1000; \
        printf("%s: %.2f ms\n", (str), (double)__elapsed_us / 1000.0); \
    } while (0)

#else

#define MR_INIT_TIMER(start, end)
#define MR_START_TIMER(start, end)
#define MR_STOP_TIMER(start, end, str)

#endif

#endif // _MR_PERF_H