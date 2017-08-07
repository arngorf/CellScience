#ifndef PMPH_TIMER
#define PMPH_TIMER

#include <sys/time.h>
#include <time.h>

namespace DebugTimer {

int timeval_subtract(struct timeval* result,
                     struct timeval* t2,
                     struct timeval* t1
                     )
{
    unsigned int resolution=1000000;
    long int diff = (t2->tv_usec + resolution * t2->tv_sec) -
    (t1->tv_usec + resolution * t1->tv_sec) ;
    result->tv_sec = diff / resolution;
    result->tv_usec = diff % resolution;
    return (diff<0);
};



unsigned long int tick()
{
    static struct timeval timeStart, timeEnd, timeDiff;

    gettimeofday(&timeEnd, NULL);

    timeval_subtract(&timeDiff, &timeEnd, &timeStart);

    gettimeofday(&timeStart, NULL);

    unsigned long int elapsed = (timeDiff.tv_sec*1e6+timeDiff.tv_usec);

    return elapsed;
};

}; // DebugTimer

#endif //PMPH_TIMER