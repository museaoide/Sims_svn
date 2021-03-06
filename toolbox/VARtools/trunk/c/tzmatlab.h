/*********
 * _iv:   Pointer to TSivector (integer vector).
 * _im:   Pointer to TSimatrix (integer matrix).
 * _dv:   Pointer to TSdvector (double vector).
 * _dm:   Pointer to TSdmatrix (double matrix).
 * _dc:   Pointer to TSdcell (double cell -- pointer to pointer to a matrix).
 * _dcv:  Pointer to TSdcellvec (double cell -- pointer to pointer to a vector).
 * _d4:   Pointer to TSdfourth (double fourth dimension -- pointer to pointer to pointer to a matrix).
 * _dzv:  Pointer to TSdzvector (double complex vector).
 * _dzm:  Pointer to TSdzmatrix (double complex matrix).
 *
 * _s:  structure variable.
 * _ps: pointer to a structure.
 * _sv: an array of structures.
 *
 * ???????? OLD NOTATIONS ??????????
 * _v:  C row or column vector pointer.
 * _vint:  C row or column vector pointer to integer.
 * _m:  C matrix pointer.
 * _mint:  C matrix pointer to integer.
 * _m3: C 3-D matrix pointer.
 * _ppm:  C pointer to pointer to a matrix.
 * d_???_ppm:  the number of pointers that are pointed to by _ppm.
 * rv_???_ppm: a vector (with dimension d_???_ppm) pointer of numbers of rows,  each of the numbers coresponding to a pointed matrix.
 * cv_???_ppm: a vector (with dimension d_???_ppm) pointer of numbers of columns, each of the numbers coresponding to a pointed matrix.
 * d_???_v: dimension size.
 * r_???_m: numbers of rows.
 * c_???_m: number of columns.
 * r_???_m3: number of rows.
 * c_???_m3: number of columns.
 * t_???_m3: number of a third dimension.
*********/


#ifndef __TZMATLAB__
#define __TZMATLAB__
   #include <stdio.h>
   #include <stdlib.h>                      // For size_t, exit, malloc, free, qsort, EXIT_FAILURE.
   #include <memory.h>                      //For memcpy, etc.  Alternative: string.h

   //#include "mkl_cblas.h"                  // Intel MKL C Blas library.
   #include "mkl.h"   //This includes mkl_cblas.h and mkl_lapack.h.
   #include <imsl.h>       //IMSL math package.
   #include <imsls.h>    //IMSL statistical package.

   #define SWITCHTOIMSLCMATH             // define: use IMSL special functions like gammlog; undef: use my own default code if it exists.


   //-------Only one of the following for math library.--------
   #define INTELCMATHLIBRARY             // define: use Intel MKL LAPACK library; undef: use others.
   //#define IMSLCMATHLIBRARY         // define: use IMSL C Math library; undef: use others.
   //#define MATLABCMATHLIBRARY            // define: use Matlab C math library; undef: use others.

   //-------Only one of the following for math library.--------
   #define SWITCHTOINTELCMATH             // define: use Intel MKL LAPACK library; undef: use others.
   //#define SWITCHTOTZCMATH            // define: use my own C math library; undef: use others.

   //-------Only one of the following for optimization routines except that CG?_ and CSMINWEL_ can be chosen together.--------
   //#define IMSL_OPTIMIZATION             // IMSL optimization routines.
   #define CSMINWEL_OPTIMIZATION      //Sims's optimization routine.
   #define CGI_OPTIMIZATION             //Polak-Ribiere conjugate gradient method without using derivative information in performing the line minimization.
   //#define CGII_OPTIMIZATION          //NOT available yet! Pletcher-Reeves conjugate gradient method using derivative information in performing the line minimization.

   //-------Only one of the following for random number generating routines.--------
   #define IMSL_RANDOMNUMBERGENERATOR             // IMSL random number generator.
   //#define CASE2_RANDOMNUMBERGENERATOR   //Imported from the C recipe book -- case 2 and my own (Iskander) code for generating a gamma distribution.

   //-------Only one of the following statistical packages.--------
   #define IMSL_STATISTICSTOOLBOX             // IMSL statistical package.



   //-------If define: use matlab API interface; otherwise (undef), use C console.
   //#define WIN_MATLABAPI                   // define: use matlab API interface; undef: use C dos console.


   //---------------
   #ifdef MATLABCMATHLIBRARY
      #include "matlab.h"                     // For all mlf???? functions.
      #include "matrix.h"                     // For mxGetM, mxCreatDoubleMatrix, etc.
   #endif
   #ifdef WIN_MATLABAPI                    // define: use matlab API interface; undef: use C dos console.
      #include "mex.h"                       // For all mex??? calls.  Matlab API (application program interface or external interface).
      #define printf mexPrintf
      #define malloc mxMalloc
      #define calloc mxCalloc
      #define free mxFree
   #endif


   //-------------- Attributes for the real double matrix type TSdmatrix.                    --------------
   //-------------- Whenever a matrix is initialized, the default is M_GE, but nothing else. --------------
   #define M_UNDEF  0        //0 or NULL: No attribute will be given when memory is allocated but no values are initialized.
   #define M_GE     0x0001   //1:    A general matrix.
   #define M_SU     0x0002   //2:    A symmetric (must be square) matrix but only the upper triangular part is referenced.
   #define M_SL     0x0004   //4:    A symmetric (must be square) matrix but only the lower triangular part is referenced.
   #define M_UT     0x0008   //8:    A upper triangular (trapezoidal if nrows < ncols) matrix but only the upper triangular part is referenced.
   #define M_LT     0x0010   //16:   A lower triangular (trapezoidal if nrows > ncols) matrix but only the lower triangular part is referenced.
   #define M_CN     0x0020   //32:   A constant (CN) matrix (All elements are the same or no (N) change from one to another).
//   #define M_UTU    0x0040   //2^6:  An unit upper triangular matrix.
//   #define M_LTU    0x0080   //2^7:  An unit lower triangular matrix.
   //-------------- Attributes for the real double vector type TSdvector.  --------------
   #define V_UNDEF 0                   //Zero or NULL: No values have been assigned to the double vector.
   #define V_DEF   1                   //True: Values have been assigned to the double vector.


   //-------------- Other macro definitions.                    --------------
   #define BLOCKSIZE_FOR_INTEL_MKL 64   //A machine-dependent value (typically, 16 to 64) required for optimum performance of the blocked algorithm in Intel MKL.
   #define NEARINFINITY 1.0E+300
   #define BIGREALNUMBER 1.0E+30


   //---------------
   #define tzMalloc(elt_count, type)           (type *)m_alloc((elt_count)*sizeof(type))
   #define tzCalloc(elt_count, type)           (type *)c_alloc((elt_count), sizeof(type))
   #define tzDestroy(x)   {if ((x)) { \
                             free((x)); \
                             (x) = NULL; \
                          }}
   #define tzFclose(x)   {if ((x)) { \
                             fclose((x)); \
                             (x) = (FILE *)NULL; \
                          }}
   #define mos(i, j, nrows)  ((j)*(nrows)+(i))   //i: ith row; j: jth column; nrows: number of rows for the matrix.
            //Offset(os) for a matrix(m) in column major order and with base 0.  See Reek pp.241-242.
   #define square(x)        ((x) * (x))    //Must be careful to avoid using, say, square(tmpd=2) or square(++x).
   #define _max(a, b)  ((a)>(b) ? (a) : (b))  // Macro max or __max is already defined in stdlib.h in MS visual C++, but mex.h may overwrite the max macro so we use _max.
   #define _min(a, b)  ((a)>(b) ? (b) : (a))
   #define swap(a, b, stemp)  {(stemp)=(a); (a)=(b); (b)=(stemp);}


   //---------------
   typedef struct {
      int *v;  //v: vector.
      int n;
      int flag;    //flag: no legal values are assigned if 0 and legal values are assigned if 1.
   } TSivector;
	typedef struct {
      int *M;  //M: matrix.
      int nrows, ncols;
		int flag;  //flag: Refers to M_GE, M_SU, M_SL, M_UT, and M_LT in tzmatlab.h.
   } TSimatrix;
   typedef struct {
      TSivector **C;    //ncells-by-1 cells (C) and a ponter to vector in each cell.
      int ncells;  //Number of pointers (cells) to pointer.
   } TSicellvec;
   typedef struct {
      TSimatrix **C;    //ncells-by-1 cells (C) and a ponter to vector in each cell.
      int ncells;  //Number of pointers (cells) to pointer.
   } TSicell;
	//=== Real types.
	typedef struct {
      double *v;  //v: vector.
      int n;
      int flag;   //flag: no legal values are assigned if 0 and legal values are assigned if 1.
   } TSdvector;
   typedef struct {
      double *M;  //M: matrix.
      int nrows, ncols;
      int flag;   //flag: Refers to M_GE, M_SU, M_SL, M_UT, and M_LT in tzmatlab.h.
   } TSdmatrix;
   typedef struct {
      TSdmatrix **C;    //ncells-by-1 cells (C) and a pointer to matrix in each cell.
      int ncells;  //Number of pointers (cells) to pointer.
   } TSdcell;
   typedef struct {
      TSdvector **C;    //ncells-by-1 cells (C) and a ponter to vector in each cell.
      int ncells;  //Number of pointers (cells) to pointer.
   } TSdcellvec;
   typedef struct {
      TSdcell **F;    //ndims-by-1 fourths (F) and a pointer to cell in each fourth.
      int ndims;  //Number of pointers (fourth dimensions) to pointer.
   } TSdfourth;
   typedef struct {
      TSdcellvec **F;    //ndims-by-1 fourths (F) and a pointer to cellvec in each fourth.
      int ndims;  //Number of pointers (fourth dimensions) to pointer.
   } TSdfourthvec;
   //=== Complex types.
   typedef struct {
      TSdvector *real;     //Real part.
      TSdvector *imag;     //Imaginary part.
   } TSdzvector;
   typedef struct {
      TSdmatrix *real;   //Real part.
      TSdmatrix *imag;   //Imaginary part.
   } TSdzmatrix;




   //---------------
   void fn_DisplayError(char *msg_s);
   void *m_alloc(size_t size);
   void *c_alloc(size_t elt_count, size_t elt_size);
   //+
   TSivector *CreateVector_int(int _n);
   TSivector *DestroyVector_int(TSivector *x_iv);
   TSimatrix *CreateMatrix_int(int nrows, int ncols);
   TSimatrix *DestroyMatrix_int(TSimatrix *x_im);
	TSicellvec *CreateCellvec_int(TSivector *n_iv);
	TSicellvec *DestroyCellvec_int(TSicellvec *x_icv);
   TSicell *CreateCell_int(TSivector *row_iv, TSivector *col_iv);
   TSicell *DestroyCell_int(TSicell *x_ic);



   TSdvector *CreateVector_lf(int _n);
   TSdvector *DestroyVector_lf(TSdvector *x_iv);
   TSdmatrix *CreateMatrix_lf(int nrows, int ncols);
   TSdmatrix *DestroyMatrix_lf(TSdmatrix *x_im);
   TSdcell *CreateCell_lf(TSivector *row_iv, TSivector *col_iv);
   TSdcell *DestroyCell_lf(TSdcell *x_dc);
   TSdcellvec *CreateCellvec_lf(TSivector *n_iv);
   TSdcellvec *DestroyCellvec_lf(TSdcellvec *x_dcv);
   TSdfourth *CreateFourth_lf(int ndims, TSivector *row_iv, TSivector *col_iv);
   TSdfourth *DestroyFourth_lf(TSdfourth *x_d4);
   TSdfourthvec *CreateFourthvec_lf(int ndims, TSivector *n_iv);
   TSdfourthvec *DestroyFourthvec_lf(TSdfourthvec *x_d4v);

   TSdzvector *CreateVector_dz(int _n);
   void DestroyVector_dz(TSdzvector *x_dzv);
   TSdzmatrix *CreateMatrix_dz(int nrows, int ncols);
   void DestroyMatrix_dz(TSdzmatrix *x_dzm);

   //+
   TSdmatrix *CreateZeroMatrix_lf(const int nrows, const int ncols);
   TSdmatrix *CreateIdentityMatrix_lf(const int nrows, const int ncols);
   //TSdvector *CreateZerosVector_lf(int _n);
   TSivector *CreateConstantVector_int( const int _n, const int _k);
   TSimatrix *CreateConstantMatrix_int(const int nrows, const int ncols, const int _n);
   TSicellvec *CreateConstantCellvec_int(TSivector *n_iv, const int _n);
   TSicell *CreateConstantCell_int(TSivector *row_iv, TSivector *col_iv, const int _n);
	TSdvector *CreateConstantVector_lf(const int _n, const double _alpha);
   TSdmatrix *CreateConstantMatrix_lf(const int nrows, const int ncols, const double _alpha);
   TSdcellvec *CreateConstantCellvec_lf(TSivector *n_iv, const double _alpha);
   TSdcell *CreateConstantCell_lf(TSivector *row_iv, TSivector *col_iv, const double _alpha);
   //+
   void InitializeConstantVector_lf(TSdvector *x_dv, const double _alpha);
   void InitializeConstantVector_int(TSivector *x_iv, const int _k);
   void InitializeConstantMatrix_lf(TSdmatrix *x_dm, const double _alpha);
   void InitializeDiagonalMatrix_lf(TSdmatrix *x_dm, const double _alpha);
   void InitializeConstantMatrix_int(TSimatrix *x_dm, const int _alpha);
   void InitializeConstantCellvec_lf(TSdcellvec *x_dcv, const double _alpha);
   void InitializeConstantCell_lf(TSdcell *x_dc, const double _alpha);
   void InitializeConstantFourthvec_lf(TSdfourthvec *x_d4v, const double _alpha);
   void InitializeConstantFourth_lf(TSdfourth *x_d4, const double _alpha);


   void NegateColofMatrix_lf(TSdvector *y_dv, TSdmatrix *x_dm, int _j);
   void InitializeConstantColofMatrix_lf(TSdmatrix *X_dm, int jx, double _alpha);


   FILE *tzFopen(char *filename, char *mode);
#endif
