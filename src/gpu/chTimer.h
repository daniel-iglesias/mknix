
#ifndef __CUDAHANDBOOK_TIMER_H__
#define __CUDAHANDBOOK_TIMER_H__

#include <iostream>
#include <iomanip>
#include <stdlib.h>


#ifdef __cplusplus
extern "C" {
#endif

// On Linux and MacOS, use clock_gettime()
typedef struct timespec chTimerTimestamp;

struct cpuClock{
 chTimerTimestamp start, end;
 double elapsedMicroseconds;
};


inline void
chTimerGetTime( chTimerTimestamp *p )
{
    clock_gettime(CLOCK_MONOTONIC, p);
}

inline double
chTimerElapsedTime( chTimerTimestamp *pStart, chTimerTimestamp *pStop )
{
    return (double) (pStop->tv_sec - pStart->tv_sec) +
                    (pStop->tv_nsec - pStart->tv_nsec)/1e9;
}

/*void cpuTick(cpuClock *ck){
  chTimerGetTime(&ck->start);
}

void cpuTock(cpuClock *ck){
  chTimerGetTime(&ck->end);
  double nano = chTimerElapsedTime(&ck->start, &ck->end);
  std::cout << "Time Consumed: " << nano * 1e6 << " microseconds" <<std::endl;
}*/

//#endif
/*
inline double
chTimerBandwidth( chTimerTimestamp *pStart, chTimerTimestamp *pStop, double cBytes )
{
    double et = chTimerElapsedTime( pStart, pStop );
    return cBytes / et;
}*/

#ifdef __cplusplus
};
#endif

#endif
