
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
chTimerGetTime( chTimerTimestamp *p );

inline double
chTimerElapsedTime( chTimerTimestamp *pStart, chTimerTimestamp *pStop );

void cpuTick(cpuClock *ck);

void cpuTock(cpuClock *ck, std::string function_name=" ");

void saveClock(cpuClock *ck, std::string file_path ,std::string message=" ");

//#endif

inline double
chTimerBandwidth( chTimerTimestamp *pStart, chTimerTimestamp *pStop, double cBytes );

#ifdef __cplusplus
};
#endif

#endif//__CUDAHANDBOOK_TIMER_H__
