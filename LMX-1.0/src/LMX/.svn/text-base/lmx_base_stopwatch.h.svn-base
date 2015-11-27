/***************************************************************************
 *   Copyright (C) 2005 by Daniel Iglesias                                 *
 *   dani@localhost.localdomain                                            *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU Library General Public License as       *
 *   published by the Free Software Foundation; either version 2 of the    *
 *   License, or (at your option) any later version.                       *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU Library General Public     *
 *   License along with this program; if not, write to the                 *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef LMXSTOPWATCH_H
#define LMXSTOPWATCH_H

#include <ctime>
#include <sys/time.h>

//////////////////////////////////////////// Doxygen file documentation entry:
    /*!
      \file lmx_base_stopwatch.h

      \brief This file contains both the declaration and implementation for Stopwatch class.

      This file contains a tool for timing operations. See class lmx::Stopwatch for more details.

      \author Daniel Iglesias Ibáñez
      
    */
//////////////////////////////////////////// Doxygen file documentation (end)



namespace lmx {

  /**
   * \class Stopwatch
   * \brief Class Stopwatch
   *
   * This class is used to measure the time that it takes to do a group of operations.
   * Those operations must be consecutive. The usage is very simple:
   * just put the operations you want to measure inside a block (use braces) and create an empty Stopwatch object before any operation takes place.
   * 
   * WARNING: If using more than one core (processor) in execution, use class ExactStopWatch instead.
   * 
   * Inspired by Danny Kalev's recipe. More info in: http://www.informit.com/guides/content.asp?g=cplusplus&seqNum=156
   *
   * \author Daniel Iglesias
   * 
   */
class Stopwatch{

public:
 
/**
 * Empty constructor
 */
 Stopwatch() :
     start( std::clock() ), //start counting time
     vervose(1),
     new_out(0),
     out (&std::cout)
 {}

/**
 * Constructor with name of file
 */
 Stopwatch(char* file_name) :
     start( std::clock() ), //start counting time
     vervose(1),
     new_out(1)
 { out = new std::ofstream(file_name); }

  /**
   * Destructor
   */
 ~Stopwatch()
  {
    if(vervose){
      std::clock_t total_ticks = std::clock()-start; //get elapsed time
      *out << std::endl;
      *out << "---------------------timer info:--------------------" << std::endl;
      *out << "|                                                  |" << std::endl;
      *out << "| total of ticks for this activity: "<< total_ticks << std::endl;
      *out << "| in seconds: "<<double(total_ticks)/CLOCKS_PER_SEC<<std::endl;
      *out << "|                                                  |" << std::endl;
      *out << "-------------------end timer info.------------------" << std::endl;
      *out << std::endl;
    }
    if (new_out) {
      delete out;
      out = 0;
    }
  }

  /**
   * Time count access.
   *
   * @return value of the counter.
   */
  double getTime()
  { std::clock_t total = std::clock()-start;
    return double(total)/CLOCKS_PER_SEC;
  }
    
  /**
   * Turns off the default time info message.
   */
  void setQuiet()
  { vervose = 0;
  }

private:
 std::clock_t start;
 int vervose;
 bool new_out;
 std::ostream* out;
};



  /**
   * \class ExactStopwatch
   * \brief Class ExactStopwatch
   *
   * This class is better than class Stopwatch when using parallelized executables because it uses the system's time
   * instead of counting processor clicks. The usage is very simple: just put the operations you want to measure inside
   * a block (use braces) and create an empty Stopwatch object before any operation takes place.
   *
   * WARNING: Ther might be some portability issues.
   *
   * \author Daniel Iglesias
   * 
   */
class ExactStopwatch{

  public:
 
  /**
   * Empty constructor
   */
    ExactStopwatch() :
      start( std::clock() ), //start counting time
      vervose(1),
      new_out(0),
      out (&std::cout)
    { gettimeofday ( &t1, (struct timezone*)0 );
    }
 
  /**
  * Constructor with name of file
  */
    ExactStopwatch(char* file_name) :
      start( std::clock() ), //start counting time
      vervose(1),
      new_out(1)
    { out = new std::ofstream(file_name); }

  /**
   * Destructor
   */
    ~ExactStopwatch()
    {
      if(vervose){
        std::clock_t total_ticks = std::clock()-start; //get elapsed time
        gettimeofday ( &t2, (struct timezone*)0 );
        total_time = (t2.tv_sec-t1.tv_sec) + 1E-6*(double)(t2.tv_usec-t1.tv_usec);
        *out << std::endl;
        *out << "---------------------timer info:--------------------" << std::endl;
        *out << "|                                                  |" << std::endl;
        *out << "| total of ticks for this activity: "<< total_ticks << std::endl;
        *out << "| in seconds: "<< total_time <<std::endl;
        *out << "|                                                  |" << std::endl;
        *out << "-------------------end timer info.------------------" << std::endl;
        *out << std::endl;
      }
      if (new_out) {
        delete out;
        out = 0;
      }
    }

    /**
     * Time count access.
     * 
     * @return value of the counter.
     */
    double getTime()
    { gettimeofday ( &t2, (struct timezone*)0 );
      total_time = (t2.tv_sec-t1.tv_sec) + 1E-6*(double)(t2.tv_usec-t1.tv_usec);
      return total_time;
    }
    
    /**
     * Turns off the default time info message.
     */
    void setQuiet()
    { vervose = 0;
    }

  private:
    std::clock_t start;
    double total_time;
    int vervose;
    bool new_out;
    std::ostream* out;
    struct timeval t1, t2;
};

}

#endif
