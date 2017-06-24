

#include "chTimer.h"
#include <fstream>
#include <ios>



inline void
chTimerGetTime( chTimerTimestamp *p )
{
    clock_gettime(CLOCK_MONOTONIC, p);
}

inline double
chTimerElapsedTime( chTimerTimestamp *pStart, chTimerTimestamp *pStop )
{ //returns in seconds up to the nanosecond resolution
    return (double) (pStop->tv_sec - pStart->tv_sec) +
                    (pStop->tv_nsec - pStart->tv_nsec)/1e9;
}

void cpuTick(cpuClock *ck){
  chTimerGetTime(&ck->start);
}

void cpuTock(cpuClock *ck, std::string function_name){
  chTimerGetTime(&ck->end);
  double nano = chTimerElapsedTime(&ck->start, &ck->end);
  ck->elapsedMicroseconds = nano * 1e6 ;
  std::cout << function_name <<" Time Consumed: " << ck->elapsedMicroseconds  << " microseconds" <<std::endl;
}

//#endif

inline double
chTimerBandwidth( chTimerTimestamp *pStart, chTimerTimestamp *pStop, double cBytes )
{
    double et = chTimerElapsedTime( pStart, pStop );
    return cBytes / et;
}

void saveClock(cpuClock *ck, std::string file_path ,std::string message){
    std::ofstream s_out(file_path.c_str());
    s_out << message << std::endl;
    s_out << message <<" Time Consumed: " << ck->elapsedMicroseconds  << " microseconds" <<std::endl;
    s_out.close();
}

void clockFullStats(std::vector<double> &clockResultsVector,
                    std::string function_name)
{
  double clockAvg = 0.0;
  int clockMeasures = clockResultsVector.size();
  for(int i = 0; i < clockMeasures; i++)
        clockAvg += clockResultsVector.at(i);
  clockAvg /= clockMeasures;

  std::cout << ">>  Clock measures for " << function_name << std::endl;
  std::cout << ">>  Number of measures " << clockMeasures << std::endl;
  std::cout << ">>  Average time " << clockAvg << " microseconds" << std::endl;
  std::cout << std::endl;
}
