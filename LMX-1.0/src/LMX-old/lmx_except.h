/***************************************************************************
 *   Copyright (C) 2006 by Daniel Iglesias                                 *
 *   http://code.google.com/p/lmx                                          *
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

    /* ******************************************************************* */
    /*                                                                     */
    /* Based on dal_except.h in gmm library:                               */
    /*                                                                     */
    /* ******************************************************************* */
    /*                                                                     */
    /* Library :  Dynamic Array Library (dal)                              */
    /* File    :  dal_except.h : Exceptions.                               */
    /*     									                                               */
    /*                                                                     */
    /* Date : September 01, 2002.                                          */
    /* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                  */
    /*          Julien Pommier, Julien.pommier@gmm.insa-tlse.fr            */
    /*                                                                     */
    /* Copyright (C) 2002  Yves Renard.                                    */
    /*                                                                     */
    /* ******************************************************************* */

#ifndef LMX_EXCEPT_H__
#define LMX_EXCEPT_H__

#include <iostream>
#include <exception>
#include <stdexcept>
#include <string>
#include <sstream>
#include <algorithm>

using std::endl; using std::cout; using std::cerr;
using std::ends; using std::cin;

namespace lmx {

//////////////////////////////////////////// Doxygen file documentation entry:
    /**
     * \file lmx_except.h
     *
     * \brief Exception handling classes, functions and macros.
     *
     * \author Addapted by Daniel Iglesias
     * 
     */
//////////////////////////////////////////// Doxygen file documentation (end)

/* *********************************************************************** */
/*	LMX generic errors.                                                    */
/* *********************************************************************** */

  /* errors definitions  */

//   using std::invalid_argument;

    /**
     * \class dimension_error
     * \brief Manages matrices and vectors dimension errors.
     */
  class dimension_error : public std::logic_error {
  public:
    /**
     * Constructor with one argument.
     * @param what_arg String containing error's information.
     */
    dimension_error(const std::string& what_arg): std::logic_error (what_arg)
      { }
  };

    /**
     * \class file_not_found_error
     * \brief Manages read/write methods' errors when checking file names.
     */
  class file_not_found_error : public std::logic_error {
  public:
    /**
     * Constructor with one argument.
     * @param what_arg String containing error's information.
     */
    file_not_found_error(const std::string& what_arg)
      : std::logic_error (what_arg) { }
  };

    /**
     * \class internal_error
     * \brief Manages errors that occur inside LMX methods and functions.
     */
  class internal_error : public std::logic_error {
  public:
    /**
     * Constructor with one argument.
     * @param what_arg String containing error's information.
     */
    internal_error(const std::string& what_arg) : std::logic_error (what_arg)
      { }
  };

    /**
     * \class failure_error
     * \brief Manages other errors not classified above.
     */
  class failure_error : public std::logic_error {
  public:
    /**
     * Constructor with one argument.
     * @param what_arg String containing error's information.
     */
    failure_error(const std::string& what_arg) : std::logic_error (what_arg)
      { }
  };

//   class not_linear_error : public std::logic_error {
//   public:
//     not_linear_error(const std::string& what_arg) : std::logic_error (what_arg)
//       { }
//   };

    /**
     * \class to_be_done_error
     * \brief Manages errors in methods witch implementation is not finished.
     */
  class to_be_done_error : public std::logic_error {
  public:
    /**
     * Constructor with one argument.
     * @param what_arg String containing error's information.
     */
    to_be_done_error(const std::string& what_arg) : std::logic_error (what_arg)
      { }
  };

  /** Function for catching any error of STL error's hierarchy. */ 
  #define LMX_STANDARD_CATCH_ERROR   catch(std::logic_error e) \
  { \
    cerr << "============================================\n";\
    cerr << "|      An error has been detected !!!      |\n";\
    cerr << "============================================\n";\
    cerr << e.what() << endl << endl;\
  }\
  catch(std::runtime_error e)\
  {\
    cerr << "============================================\n";\
    cerr << "|      An error has been detected !!!      |\n";\
    cerr << "============================================\n";\
    cerr << e.what() << endl << endl;\
    exit(1);\
  }\
  catch(std::bad_alloc) {\
    cerr << "============================================\n";\
    cerr << "|  A bad allocation has been detected !!!  |\n";\
    cerr << "============================================\n";\
    exit(1);\
  }\
  catch(std::bad_typeid) { \
    cerr << "============================================\n";\
    cerr << "|  A bad typeid     has been detected !!!  |\n";\
    cerr << "============================================\n";\
    exit(1);\
  } \
  catch(std::bad_exception) { \
    cerr << "============================================\n";\
    cerr << "|  A bad exception  has been detected !!!  |\n";\
    cerr << "============================================\n";\
    exit(1);\
  } \
  catch(std::bad_cast) { \
    cerr << "============================================\n";\
    cerr << "|    A bad cast  has been detected !!!     |\n";\
    cerr << "============================================\n";\
    exit(1);\
  } \
  catch(...) {\
    cerr << "============================================\n";\
    cerr << "|  An unknown error has been detected !!!  |\n";\
    cerr << "============================================\n";\
    exit(1);\
  }
//   catch(ios_base::failure) { 
//     cerr << "============================================\n";
//     cerr << "| A ios_base::failure has been detected !!!|\n";
//     cerr << "============================================\n";
//     exit(1);
//   } 

    /**
     * \struct exception_callback
     * \brief Callback handler for lmx exceptions.
     */
  struct exception_callback {
/// \cond EXCEPT
    virtual ~exception_callback(){}
    virtual void callback(const std::string&) = 0; //{};

    static exception_callback *which_except(exception_callback *p = 0) {
      static exception_callback *exc_cback = 0;
      if (p != 0) exc_cback = p;
      return exc_cback;
    }

    static void do_exception_callback(const std::string &msg)
      { if (which_except()) which_except()->callback(msg); }

    static void set_exception_callback(exception_callback *e)
      { which_except(e); }
/// \endcond

  };

    /**
     * \struct exception_callback_debug
     * \brief crashing callback for debug mode.
     */
  struct exception_callback_debug : public lmx::exception_callback  {
/// \cond EXCEPT
    /*virtual */void callback(const std::string& msg)
    { cerr << msg << endl; *(int *)(0) = 0; }
/// \endcond
  };

  /**
   * User's function for changing the default exception callback
   * @param e Pointer to the new exception callback.
   */
  inline void set_exception_callback(exception_callback *e)
  { exception_callback::which_except(e); }

#ifdef HAVE_PRETTY_FUNCTION
  /** Variable for storing the full function name where the error occurs. */ 
#  define LMX_PRETTY_FUNCTION __PRETTY_FUNCTION__
#else 
  /** Variable for storing the full function name where the error occurs. */
#  define LMX_PRETTY_FUNCTION ""
#endif

  /** standard function for throwing errors */ 
#define LMX_THROW(type, thestr) {                                    \
    std::stringstream msg;                                           \
    msg << "Error in " __FILE__ << ", line "                          \
        << __LINE__ << " " << LMX_PRETTY_FUNCTION << ": \n" << thestr << ends; \
    lmx::exception_callback::do_exception_callback(msg.str());       \
    throw (type)(msg.str());                                         \
  }

#ifdef DEBUG_MODE
#  define LMX_INTERNAL_ERROR(thestr) { \
  cerr << "Internal error: " << LMX_PRETTY_FUNCTION << " " << thestr << endl; \
   ::abort(); \
   }
#else
  /** standard function for throwing internal errors */ 
#  define LMX_INTERNAL_ERROR(thestr) LMX_THROW(lmx::internal_error, "Internal error: " << thestr)
#endif


    /**
     * \struct warning_level
     * \brief Manages importance level in warnings.
     */
  struct warning_level {
    /**
     * Function for setting warning level.
     * @param l Importance level.
     * @return level's value.
     */
    static int level(int l = -2)
    { static int level_ = 3;
      return (l != -2) ? (level_ = l) : level_;
    }
  };

  /** user function for changing the level warning. */
  inline void set_warning_level(int l) { warning_level::level(std::max(0,l)); }

  /** user function for throwing warning messages. */
#define LMX_WARNING(level_, thestr) {                                 \
    std::stringstream msg;                                            \
    msg << "Level " << level_ << " Warning in " __FILE__ << ", line "  \
        << __LINE__ << " " << LMX_PRETTY_FUNCTION << ": " << thestr << ends; \
    if ((level_) <= lmx::warning_level::level())                      \
       std::cerr << msg.str() << std::endl;                           \
  } 

  // Warning levels : 0 always printed
  //                  1 very important : specify a possible error in the code.
  //                  2 important : specify a default of optimization for inst.
  //                  3 remark
  //                  4 ignored by default.
  
} /* end of namespace lmx.                                                */



#endif /* LMX_EXCEPT_H__ */
