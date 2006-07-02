#ifndef __CSMINWEL_H__
#define __CSMINWEL_H__

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <float.h>

#include "tzmatlab.h"

void csminwel(double (*fcn)(double *x, int n, double **args, int *dims),
            double *x, int n, double *H, double *gh,
            int (*grad)(double *x, int n, double *g, double **args, int *dims),
            double *fh, double crit, int *itct, int nit,
            int *fcount, int *retcodeh, double **args, int *dims);
// Alternative but less clear way:  ... (double (*fcn)(double *, int, double **, int *), ...

void csminwel_SetPrintFile(char *filename);

#endif
