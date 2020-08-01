
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mex.h"


typedef char			si1;
typedef unsigned char		ui1;
typedef short			si2;
typedef unsigned short		ui2;
typedef int		        si4;
typedef unsigned int	        ui4;
typedef long    		si8;
typedef unsigned long 	        ui8;
typedef float			sf4;
typedef double			sf8;

#define ABS(x)  ( ((x) >= 0) ? (x) : -(x) )


void    lad_reg(sf8 *y, sf8 *ba, sf8 *ma, si8 n);
sf8     quantval(sf8 *x, si8 len, sf8 quantile, ui1 preserve_input, sf8 *buff);
sf8     kth_smallest(sf8 *a, si8 n, si8 k, sf8 lo_p);
void	mexFunction(si4 nlhs, mxArray *plhs[], si4 nrhs, const mxArray *prhs[]);


void    lad_reg(sf8 *y, sf8 *ba, sf8 *ma, si8 n)
{
        sf8     b, m, t, *yp, *buff, *bp, min_y, max_y, min_m, max_m, thresh, m_sum;
        sf8     d, m_eps, b_eps, lad_eps, test_m, lad, upper_m, lower_m, safe_eps;
        si8     i, mid_idx, n_bytes;
        
        
        // least absolute differences linear regression
        // assumes x to be 1:length(x)
        // fit: y = mx + b
        
        // allocate
        buff = (sf8 *) calloc((size_t) n, sizeof(sf8));
        if (buff == NULL) {
                mexErrMsgTxt("[lad_reg] could not allocate enough memory");
                return;
        }
        
        // setup
        yp = y;
        min_y = max_y = *yp;
        for (i = n; --i;) {
                if (*++yp > max_y)
                        max_y = *yp;
                else if (*yp < min_y)
                        min_y = *yp;
        }
        lower_m = min_m = (min_y - max_y) / (sf8) (n - 1);
        upper_m = max_m  = -min_m;
        safe_eps = DBL_EPSILON * (sf8) 1000.0;
        thresh = safe_eps * (sf8) 10.0;
        d = max_m - min_m;
        mid_idx = (n - 1) >> 1;
        n_bytes = n * sizeof(sf8);
        
        // search
        while (d > thresh) {
                m = (upper_m + lower_m) / (sf8) 2.0;
                bp = buff; yp = y;
                m_sum = (sf8) 0.0;
                for (i = n; i--;)
                        *bp++ = *yp++ - (m_sum += m);
                b = quantval(buff, n, 0.5, 0, NULL);  // median
                bp = buff; lad = (sf8) 0.0;
                for (i = n; i--;) {
                        t = *bp++ - b;
                        lad += ABS(t);
                }
                m_eps = m + safe_eps;
                bp = buff; yp = y;
                m_sum = (sf8) 0.0;
                for (i = n; i--;)
                        *bp++ = *yp++ - (m_sum += m_eps);
                b_eps = quantval(buff, n, 0.5, 0, NULL);  // median
                bp = buff; lad_eps = (sf8) 0.0;
                for (i = n; i--;) {
                        t = *bp++ - b_eps;
                        lad_eps += ABS(t);
                }
                test_m = lad_eps - lad;
                if (test_m > (sf8) 0.0)
                        upper_m = m;
                else if (test_m < (sf8) 0.0)
                        lower_m = m;
                else
                        break;
                d = upper_m - lower_m;
        }
        *ba = b;
        *ma = m;
        
        // clean up
        free(buff);
        
        return;
}


sf8     quantval(sf8 *x, si8 len, sf8 quantile, ui1 preserve_input, sf8 *buff)
{
        ui1     free_buff;
        sf8     q, fk, lo_p;
        si8     lo_k;
        
        
        if (len == 1)
                return(*x);
        
        free_buff = 0;
        if (preserve_input) {
                if (buff == NULL) {
                        buff = (sf8 *) malloc((size_t) len * sizeof(sf8));
                        if (buff == NULL) {
                                fprintf(stderr, "%s(): Not enough memory => exiting (line %d)\n", __FUNCTION__, __LINE__);
                                exit(1);
                        }
                        free_buff = 1;
                }
                memcpy(buff, x, len * sizeof(sf8));
        } else {
                buff = x;
        }
        
        if (quantile == (sf8) 1.0) {
                lo_k = len - 2;
                lo_p = (sf8) 0.0;
        } else {
                fk = quantile * (sf8) (len - 1);
                lo_k = (si8) fk;
                lo_p = (sf8) 1.0 - (fk - (sf8) lo_k);
        }
        
        if (len == 2) {
                if (x[0] <= x[1])
                        return((x[0] * lo_p) + (x[1] * (1.0 - lo_p)));
                return((x[1] * lo_p) + (x[0] * (1.0 - lo_p)));
        }

        q = kth_smallest(buff, len, lo_k, lo_p);

        if (free_buff)
                free(buff);
        
        return(q);
}


// kth_smallest()
// Algorithm from Niklaus Wirth's book: "Algorithms + data structures = programs".
// Code here is derived from code by Nicolas Devillard. Public domain.
sf8     kth_smallest(sf8 *x, si8 n, si8 lo_k, sf8 lo_p)
{
        sf8             lo_v, *lp, *mp, *last_mp, *lo_kp, *hi_kp;
        register sf8    v, t, *xip, *xjp;
        
                
        lp = x;
        last_mp = mp = x + (n - 1);
        lo_kp = x + lo_k;
        hi_kp = lo_kp + 1;
        while (lp < mp) {
                v = *lo_kp;
                xip = lp;
                xjp = mp;
                do {
                        for (; *xip < v; ++xip);
                        for (; v < *xjp; --xjp);
                        if (xip <= xjp) {
                                t = *xip;
                                *xip++ = *xjp;
                                *xjp-- = t;
                        }
                } while (xip <= xjp);
                
                if (xjp < lo_kp)
                        lp = xip;
                if (hi_kp < xip)
                        last_mp = mp;
                if (lo_kp < xip)
                        mp = xjp;
        }
        lo_v = *lo_kp;
     
        lp = lo_kp; mp = last_mp;
        while (lp < mp) {
                v = *hi_kp;
                xip = lp;
                xjp = mp;
                do {
                        for (; *xip < v; ++xip);
                        for (; v < *xjp; --xjp);
                        if (xip <= xjp) {
                                t = *xip;
                                *xip++ = *xjp;
                                *xjp-- = t;
                        }
                } while (xip <= xjp);

                if (xjp < hi_kp)
                        lp = xip;
                if (hi_kp < xip)
                        mp = xjp;
        }
        
        return((lo_v * lo_p) + (*hi_kp * ((sf8) 1.0 - lo_p)));
}


// Mex gateway routine
void    mexFunction(si4 nlhs, mxArray *plhs[], si4 nrhs, const mxArray *prhs[])
{
        sf8     *y, *b, *m;
        si8     n, cols;
        ui8     dims[2];
        
        
        //  Check for proper number of arguments
        if (nrhs != 1) {
                mexErrMsgTxt("[lad_reg] one input required: input_vector");
                return;
        }
        if (nlhs != 2) {
                mexErrMsgTxt("[lad_reg] two oubputs required: b & m");
                return;
        }
        
        // Check to make sure the input argument is a double vector
        if (mxIsDouble(prhs[0]) != 1) {
                mexErrMsgTxt("[lad_reg] input_vector must be of type double");
                return;
        }
        // Check to make sure the first input argument is a column vector
        cols = (si8) mxGetN(prhs[0]);
        if (cols > 1) {
                mexErrMsgTxt("[lad_reg] input_vector must be a column vector");
                return;
        }
        // Get the rows (length) of the input vector
        n = (si8) mxGetM(prhs[0]);
        // Get a c pointer to the input vector
        y = (sf8 *) mxGetPr(prhs[0]);
        
        // Allocate the first oubput scalar (b)
        dims[0] = 1; dims[1] = 1;
        plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
        
        // Create a c pointer to first oubput scalar (b)
        b = (sf8 *) mxGetPr(plhs[0]);
        
        // Allocate the second oubput scalar (m)
        dims[0] = 1; dims[1] = 1;
        plhs[1] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
        
        // Create a c pointer to second oubput scalar (m)
        m = (sf8 *) mxGetPr(plhs[1]);
        
        // Call the c subroutine
        lad_reg(y, b, m, n);
        
        return;
}

