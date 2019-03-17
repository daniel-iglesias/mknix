/***************************************************************************
 *   Copyright (C) 2005 by Daniel Iglesias                                 *
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
 ***************************************************************************/#ifndef IOHB_H
#define IOHB_H

#include <stdio.h>
#include<complex>

//////////////////////////////////////////// Doxygen file documentation entry:
/*!
  \file lmx_base_iohb.h

  \brief This file contains functions for Harwell Boeing format file reading and writing.

  Addapted from library Harwell-Boeing File I/O in C, V. 1.0.

  \author Adapted by Daniel Iglesias

*/
//////////////////////////////////////////// Doxygen file documentation (end)

namespace lmx {

/*************************************************************************/
/*                                                                       */
/*  Functions to read and write Harwell Boeing format.                   */
/*                                                                       */
/*************************************************************************/

// Fri Aug 15 16:29:47 EDT 1997
//
//                      Harwell-Boeing File I/O in C
//                               V. 1.0
//
//          National Institute of Standards and Technology, MD.
//                            K.A. Remington
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//                                NOTICE
//
// Permission to use, copy, modify, and distribute this software and
// its documentation for any purpose and without fee is hereby granted
// provided that the above copyright notice appear in all copies and
// that both the copyright notice and this permission notice appear in
// supporting documentation.
//
// Neither the Author nor the Institution (National Institute of Standards
// and Technology) make any representations about the suitability of this
// software for any purpose. This software is provided "as is" without
// expressed or implied warranty.
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


/// \cond IOHB
inline void IOHBTerminate(const char * a) {LMX_THROW(lmx::failure_error, a); }

inline bool is_complex_double__(std::complex<double>) { return true; }

inline bool is_complex_double__(double) { return false; }

inline int ParseIfmt(const char * fmt, int * perline, int * width)
{
    if (sscanf(fmt, " (%dI%d)", perline, width) != 2) LMX_THROW(lmx::failure_error, "invalid HB I-format : " << fmt);
    return *width;
}

inline int ParseRfmt(const char * fmt, int * perline, int * width,
                     int * prec, int * flag)
{
    char p;
    *perline = *width = *flag = *prec = 0;
    if (sscanf(fmt, " (%d%c%d.%d)", perline, &p, width, prec) < 3 ||
        !strchr("PEDF", p)) LMX_THROW(lmx::failure_error, "invalid HB REAL format : " << fmt);
    *flag = p;
    return *width;
}
/// \endcond

/** matrix input/output for Harwell-Boeing format */
struct HarwellBoeing_IO
{
/// \cond IOHB
    int nrows() const { return Nrow; }

    int ncols() const { return Ncol; }

    int nnz() const { return Nnzero; }

    int is_complex() const { return Type[0] == 'C'; }

    int is_symmetric() const { return Type[1] == 'S'; }

    int is_hermitian() const { return Type[1] == 'H'; }

    HarwellBoeing_IO() { clear(); }

    HarwellBoeing_IO(const char * filename)
    {
        clear();
        open(filename);
    }

    ~HarwellBoeing_IO() { close(); }

    /* open filename and reads header */
    void open(const char * filename);

    template<typename T>
    void read(int& M, int& N, int& nonzeros, int *& colptr, int *& rowind, T *& val);
    /* read the opened file */
/*    template <typename T, int shift> void read(csc_matrix<T, shift>& A);
    template <typename MAT> void read(MAT &M);*/
    /* save the matrix */
//     template <typename T, int shift> static void write(const char *filename, const csc_matrix<T, shift>& A);
//     template <typename MAT> static void write(const char *filename, const MAT& A);
/// \endcond
private:
    FILE * f;
    char Title[73], Key[9], Rhstype[4], Type[4];
    int Nrow, Ncol, Nnzero, Nrhs;
    char Ptrfmt[17], Indfmt[17], Valfmt[21], Rhsfmt[21];
    int Ptrcrd, Indcrd, Valcrd, Rhscrd;
    int lcount;


    void close()
    {
        if (f) fclose(f);
        clear();
    }

    void clear()
    {
        Nrow = Ncol = Nnzero = Nrhs = 0;
        f = 0;
        lcount = 0;
        memset(Type, 0, sizeof Type);
        memset(Key, 0, sizeof Key);
        memset(Title, 0, sizeof Title);
    }

    char * getline(char * buf)
    {
        //char * junk = fgets(buf, BUFSIZ, f);
        ++lcount;
        if (sscanf(buf, "%*s") < 0) LMX_THROW(lmx::failure_error, "blank line in HB file at line " << lcount);
        return buf;
    }

    int substrtoi(const char * p, size_type len)
    {
        char s[100];
        len = std::min(len, sizeof s - 1);
        strncpy(s, p, len);
        s[len] = 0;
        return atoi(s);
    }

    double substrtod(const char * p, size_type len, int Valflag)
    {
        char s[100];
        len = std::min(len, sizeof s - 1);
        strncpy(s, p, len);
        s[len] = 0;
        if (Valflag != 'F' && !strchr(s, 'E')) {
            /* insert a char prefix for exp */
            auto last = strlen(s);
            for (auto j = last + 1; j >= 0; j--) {
                s[j] = s[j - 1];
                if (s[j] == '+' || s[j] == '-') {
                    s[j - 1] = (char)Valflag;
                    break;
                }
            }
        }
        return atof(s);
    }

    template<typename IND_TYPE>
    int readHB_data(IND_TYPE colptr[], IND_TYPE rowind[],
                    double val[])
    {
        /************************************************************************/
        /*  This function opens and reads the specified file, interpreting its  */
        /*  contents as a sparse matrix stored in the Harwell/Boeing standard   */
        /*  format and creating compressed column storage scheme vectors to hold*/
        /*  the index and nonzero value information.                            */
        /*                                                                      */
        /*    ----------                                                        */
        /*    **CAVEAT**                                                        */
        /*    ----------                                                        */
        /*  Parsing real formats from Fortran is tricky, and this file reader   */
        /*  does not claim to be foolproof.   It has been tested for cases when */
        /*  the real values are printed consistently and evenly spaced on each  */
        /*  line, with Fixed (F), and Exponential (E or D) formats.             */
        /*                                                                      */
        /*  **  If the input file does not adhere to the H/B format, the  **    */
        /*  **             results will be unpredictable.                 **    */
        /*                                                                      */
        /************************************************************************/
        int i, ind, col, offset, count;
        int Ptrperline, Ptrwidth, Indperline, Indwidth;
        int Valperline = 0, Valwidth = 0, Valprec, Nentries;
        int Valflag = 0;           /* Indicates 'E','D', or 'F' float format */
        char line[BUFSIZ];

        /*  Parse the array input formats from Line 3 of HB file  */
        ParseIfmt(Ptrfmt, &Ptrperline, &Ptrwidth);
        ParseIfmt(Indfmt, &Indperline, &Indwidth);
        if (Type[0] != 'P') {          /* Skip if pattern only  */
            ParseRfmt(Valfmt, &Valperline, &Valwidth, &Valprec, &Valflag);
        }

        /*  Read column pointer array:   */
        offset = 0;         /* if base 0 storage is declared (via macro def),  */
        /* then storage entries are offset by 1            */

        for (count = 0, i = 0; i < Ptrcrd; i++) {
            getline(line);
            for (col = 0, ind = 0; ind < Ptrperline; ind++) {
                if (count > Ncol) break;
                colptr[count] = substrtoi(line + col, (size_t)Ptrwidth) - offset;
                count++;
                col += Ptrwidth;
            }
        }

        /*  Read row index array:  */
        for (count = 0, i = 0; i < Indcrd; i++) {
            getline(line);
            for (col = 0, ind = 0; ind < Indperline; ind++) {
                if (count == Nnzero) break;
                rowind[count] = substrtoi(line + col, (size_t)Indwidth) - offset;
                count++;
                col += Indwidth;
            }
        }

        /*  Read array of values:  */
        if (Type[0] != 'P') {          /* Skip if pattern only  */
            if (Type[0] == 'C') {
                Nentries = 2 * Nnzero;
            } else { Nentries = Nnzero; }

            count = 0;
            for (i = 0; i < Valcrd; i++) {
                getline(line);
                if (Valflag == 'D') {
                    // const_cast Due to aCC excentricity
                    char * p;
                    while ((p = const_cast<char *>(strchr(line, 'D')))) *p = 'E';
                }
                for (col = 0, ind = 0; ind < Valperline; ind++) {
                    if (count == Nentries) break;
                    val[count] = substrtod(line + col, (size_t)Valwidth, Valflag);
                    count++;
                    col += Valwidth;
                }
            }
        }
        return 1;
    }
};

/// \cond IOHB
inline void HarwellBoeing_IO::open(const char * filename)
{
    int Totcrd, Neltvl, Nrhsix;
    char line[BUFSIZ];
    close();
    f = fopen(filename, "r");
    if (!f) {LMX_THROW(lmx::failure_error, "could not open " << filename); }
    /* First line: */
    sscanf(getline(line), "%72c%8s", Title, Key);
    Key[8] = Title[72] = 0;
    /* Second line: */
    Totcrd = Ptrcrd = Indcrd = Valcrd = Rhscrd = 0;
    sscanf(getline(line), "%d%d%d%d%d", &Totcrd, &Ptrcrd, &Indcrd, &Valcrd, &Rhscrd);

    /* Third line: */
    Nrow = Ncol = Nnzero = Neltvl = 0;
    if (sscanf(getline(line), "%3c%d%d%d%d", Type, &Nrow, &Ncol, &Nnzero, &Neltvl) < 1) {
        IOHBTerminate("Invalid Type info, line 3 of Harwell-Boeing file.\n");
    }
    std::for_each(Type, Type + 3, (int (*)(int)) toupper);
//    std::for_each(Type, Type+3, toupper);
    /*  Fourth line:  */
    if (sscanf(getline(line), "%16c%16c%20c%20c", Ptrfmt, Indfmt, Valfmt, Rhsfmt) < 3) {
        IOHBTerminate("Invalid format info, line 4 of Harwell-Boeing file.\n");
    }
    Ptrfmt[16] = Indfmt[16] = Valfmt[20] = Rhsfmt[20] = 0;

    /*  (Optional) Fifth line: */
    if (Rhscrd != 0) {
        Nrhs = Nrhsix = 0;
        if (sscanf(getline(line), "%3c%d%d", Rhstype, &Nrhs, &Nrhsix) != 1) {
            IOHBTerminate("Invalid RHS type information, line 5 of"
                                  " Harwell-Boeing file.\n");
        }
    }
}

/* only valid for double and complex<double> csc matrices */
template<typename T>
void
HarwellBoeing_IO::read(int& M, int& N, int& nonzeros, int *& colptr, int *& rowind, T *& val)
{

    typedef int IND_TYPE;

    if (!f) LMX_THROW(lmx::failure_error, "no file opened!");
    if (Type[0] == 'P') LMX_THROW(lmx::failure_error, "Bad HB matrix format (pattern matrices not supported)");
    if (is_complex_double__(T()) && Type[0] == 'R') LMX_THROW(lmx::failure_error,
                                                              "Bad HB matrix format (file contains a REAL matrix)");
    if (!is_complex_double__(T()) && Type[0] == 'C') LMX_THROW(lmx::failure_error,
                                                               "Bad HB matrix format (file contains a COMPLEX matrix)");

    N = ncols();
    M = nrows();
    nonzeros = nnz();
    val = 0;
    colptr = new IND_TYPE[ncols() + 1];
    rowind = new IND_TYPE[nnz()];
    val = new T[nnz()];
    readHB_data(colptr, rowind, (double *) val);
}

/*  template <typename MAT> void 
  HarwellBoeing_IO::read(MAT &M) {
    csc_matrix<typename gmm::linalg_traits<MAT>::value_type> csc;
    read(csc); 
    resize(M, mat_nrows(csc), mat_ncols(csc));
    copy(csc, M);
  }*/

//   template <typename IND_TYPE> 
//   inline int writeHB_mat_double(const char* filename, int M, int N, int nz,
// 				const IND_TYPE colptr[],
// 				const IND_TYPE rowind[], 
// 				const double val[], int Nrhs,
// 				const double /*rhs*/[], const double /*guess*/[],
// 				const double /*exact*/[], const char* Title,
// 				const char* Key, const char* Type, 
// 				const char* Ptrfmt, const char* Indfmt,
// 				const char* Valfmt, const char* Rhsfmt,
// 				const char* Rhstype, int shift) {
//     /************************************************************************/
//       /*  The writeHB function opens the named file and writes the specified  */
//       /*  matrix and optional right-hand-side(s) to that file in              */
//       /*  Harwell-Boeing format.                                              */
//       /*                                                                      */
//       /*  For a description of the Harwell Boeing standard, see:              */
//       /*            Duff, et al.,  ACM TOMS Vol.15, No.1, March 1989          */
//       /*                                                                      */
//       /************************************************************************/
//       FILE *out_file;
//       int i,entry,offset/* , j, acount, linemod */;
//       int totcrd, ptrcrd, indcrd, valcrd, rhscrd;
//       int nvalentries, nrhsentries;
//       int Ptrperline, Ptrwidth, Indperline, Indwidth;
//       int Rhsperline, Rhswidth, Rhsprec, Rhsflag;
//       int Valperline, Valwidth, Valprec;
//       int Valflag;           /* Indicates 'E','D', or 'F' float format */
//       char pformat[16],iformat[16],vformat[19],rformat[19];
//     
//       if ( Type[0] == 'C' )
// 	{ nvalentries = 2*nz; nrhsentries = 2*M; }
//       else
// 	{ nvalentries = nz; nrhsentries = M; }
//     
//       if ( filename != NULL ) {
// 	if ( (out_file = fopen( filename, "w")) == NULL )
// 	  DAL_THROW(gmm::failure_error,"Error: Cannot open file: " << filename);
//       } else out_file = stdout;
//     
//       if ( Ptrfmt == NULL ) Ptrfmt = "(8I10)";
//       ParseIfmt(Ptrfmt, &Ptrperline, &Ptrwidth);
//       sprintf(pformat,"%%%dd",Ptrwidth);
//       ptrcrd = (N+1)/Ptrperline;
//       if ( (N+1)%Ptrperline != 0) ptrcrd++;
//     
//       if ( Indfmt == NULL ) Indfmt =  Ptrfmt;
//       ParseIfmt(Indfmt, &Indperline, &Indwidth);
//       sprintf(iformat,"%%%dd",Indwidth);
//       indcrd = nz/Indperline;
//       if ( nz%Indperline != 0) indcrd++;
//     
//       if ( Type[0] != 'P' ) {          /* Skip if pattern only  */
// 	if ( Valfmt == NULL ) Valfmt = "(4E20.13)";
// 	ParseRfmt(Valfmt, &Valperline, &Valwidth, &Valprec, &Valflag);
// 	if (Valflag == 'D') *strchr(Valfmt,'D') = 'E';
// 	if (Valflag == 'F')
// 	  sprintf(vformat, "%% %d.%df", Valwidth, Valprec);
// 	else
// 	  sprintf(vformat, "%% %d.%dE", Valwidth, Valprec);
// 	valcrd = nvalentries/Valperline;
// 	if ( nvalentries%Valperline != 0) valcrd++;
//       } else valcrd = 0;
//     
//       if ( Nrhs > 0 ) {
// 	if ( Rhsfmt == NULL ) Rhsfmt = Valfmt;
// 	ParseRfmt(Rhsfmt,&Rhsperline,&Rhswidth,&Rhsprec, &Rhsflag);
// 	if (Rhsflag == 'F')
// 	  sprintf(rformat,"%% %d.%df",Rhswidth,Rhsprec);
// 	else
// 	  sprintf(rformat,"%% %d.%dE",Rhswidth,Rhsprec);
// 	if (Rhsflag == 'D') *strchr(Rhsfmt,'D') = 'E';
// 	rhscrd = nrhsentries/Rhsperline; 
// 	if ( nrhsentries%Rhsperline != 0) rhscrd++;
// 	if ( Rhstype[1] == 'G' ) rhscrd+=rhscrd;
// 	if ( Rhstype[2] == 'X' ) rhscrd+=rhscrd;
// 	rhscrd*=Nrhs;
//       } else rhscrd = 0;
//     
//       totcrd = 4+ptrcrd+indcrd+valcrd+rhscrd;
//     
//     
//       /*  Print header information:  */
//     
//       fprintf(out_file,"%-72s%-8s\n%14d%14d%14d%14d%14d\n",Title, Key, totcrd,
// 	      ptrcrd, indcrd, valcrd, rhscrd);
//       fprintf(out_file,"%3s%11s%14d%14d%14d\n",Type,"          ", M, N, nz);
//       fprintf(out_file,"%-16s%-16s%-20s", Ptrfmt, Indfmt, Valfmt);
//       //     if ( Nrhs != 0 ) {
//       //       /*    Print Rhsfmt on fourth line and                                 */
//       //       /*      optional fifth header line for auxillary vector information:  */
//       //       fprintf(out_file,"%-20s\n%-14s%d\n",Rhsfmt,Rhstype,Nrhs);
//       //     } else
//       fprintf(out_file,"\n");
//     
//       offset = 1 - shift;  /* if base 0 storage is declared (via macro def), */
//       /* then storage entries are offset by 1           */
//     
//       /*  Print column pointers:   */
//       for (i = 0; i < N+1; i++) {
// 	entry = colptr[i]+offset;
// 	fprintf(out_file,pformat,entry);
// 	if ( (i+1)%Ptrperline == 0 ) fprintf(out_file,"\n");
//       }
//     
//       if ( (N+1) % Ptrperline != 0 ) fprintf(out_file,"\n");
//     
//       /*  Print row indices:       */
//       for (i=0;i<nz;i++) {
// 	entry = rowind[i]+offset;
// 	fprintf(out_file,iformat,entry);
// 	if ( (i+1)%Indperline == 0 ) fprintf(out_file,"\n");
//       }
//     
//       if ( nz % Indperline != 0 ) fprintf(out_file,"\n");
//     
//       /*  Print values:            */
//     
//       if ( Type[0] != 'P' ) {          /* Skip if pattern only  */
// 	for (i=0;i<nvalentries;i++) {
// 	  fprintf(out_file,vformat,val[i]);
// 	  if ( (i+1)%Valperline == 0 ) fprintf(out_file,"\n");
// 	}
// 	if ( nvalentries % Valperline != 0 ) fprintf(out_file,"\n");
//       }
//     
//       if ( fclose(out_file) != 0) {
// 	DAL_THROW(gmm::failure_error,"Error closing file in writeHB_mat_double().");
//       } else return 1;
//     }
// 
//   template <typename T, int shift> void
//   HarwellBoeing_IO::write(const char *filename, const csc_matrix<T, shift>& A) {
//     const char *t = 0;    
//     if (is_complex_double__(T()))
//       if (mat_nrows(A) == mat_ncols(A)) t = "CUA"; else t = "CRA";
//     else
//       if (mat_nrows(A) == mat_ncols(A)) t = "RUA"; else t = "RRA";
//     writeHB_mat_double(filename, mat_nrows(A), mat_ncols(A),
// 		       A.jc[mat_ncols(A)], A.jc, A.ir,
// 		       (double *)A.pr,
// 		       0, 0, 0, 0, "GETFEM++ CSC MATRIX", "CSCMAT",
// 		       t, 0, 0, 0, 0, "F", shift);
//   }
// 
//   template <typename MAT> void
//   HarwellBoeing_IO::write(const char *filename, const MAT& A) {
//     gmm::csc_matrix<typename gmm::linalg_traits<MAT>::value_type> 
//       tmp(gmm::mat_nrows(A), gmm::mat_ncols(A));
//     gmm::copy(A,tmp); 
//     HarwellBoeing_IO::write(filename, tmp);
//   }
//   
// 
//   /** save a "double" or "std::complex<double>" matrix into a HarwellBoeing file */
//   template <typename T, int shift> inline void
//   Harwell_Boeing_save(const char *filename, const csc_matrix<T, shift>& A) {
//     HarwellBoeing_IO h; h.write(filename, A);
//   }
// 
//  /** load a "double" or "std::complex<double>" matrix from a HarwellBoeing file */
//   template <typename T, int shift> void
//   Harwell_Boeing_load(const char *filename, csc_matrix<T, shift>& A) {
//     HarwellBoeing_IO h(filename); h.read(A);
//   }
/// \endcond

}; // namespace lmx



#endif
