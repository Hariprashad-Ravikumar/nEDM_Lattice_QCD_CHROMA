// -*- C++ -*-
/*! @file
 * @brief Timer support
 *
 * A stopwatch like timer.
 */

#ifndef ADAT_STOPWATCH_H
#define ADAT_STOPWATCH_H

#include<sys/time.h>

namespace Util 
{

//! Timer
  class StopWatch 
  {
  public:
    //! Constructor
    StopWatch();

    //! Destructor
    ~StopWatch();

    //! Reset the timer
    void reset();

    //! Start the timer
    void start();

    //! Stop the timer
    void stop();

    //! Get time in microseconds
    double getTimeInMicroseconds();

    //! Get time in seconds
    double getTimeInSeconds();

    //! Get split time in microseconds
    double getSplitTimeInMicroseconds();

    //! Get split time in seconds
    double getSplitTimeInSeconds();

  private:
    //! Get split time
    void splitTime();

    struct timeval t_split;

  private:
    long sec;
    long usec;
    bool startedP;
    bool stoppedP;

    struct timeval t_start;
    struct timeval t_end;
  };

} // namespace Util

#endif
