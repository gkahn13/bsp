#ifndef TIMER_H_DEF
#define TIMER_H_DEF

#include <time.h>
#include <sys/sysinfo.h>

namespace util {

typedef struct Timer{
	struct timespec tic;
	struct timespec toc;
} Timer;

/* read current time */
inline void Timer_tic(Timer* t) {
	clock_gettime(CLOCK_MONOTONIC_RAW, &t->tic);
}

/* return time passed since last call to tic on this timer */
inline double Timer_toc(Timer* t)
{
	struct timespec temp;
	clock_gettime(CLOCK_MONOTONIC_RAW, &t->toc);

	if ((t->toc.tv_nsec - t->tic.tv_nsec)<0) {
		temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec-1;
		temp.tv_nsec = 1000000000+(t->toc.tv_nsec - t->tic.tv_nsec);
	} else {
		temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec;
		temp.tv_nsec = t->toc.tv_nsec - t->tic.tv_nsec;
	}

	return (double)temp.tv_sec + (double)temp.tv_nsec / 1000000000;
}

}
#endif // TIMER_H_DEF
