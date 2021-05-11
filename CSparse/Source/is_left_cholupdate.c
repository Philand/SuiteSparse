#include "cs.h"

/* --- is_left_cholupdate */

csi *is_pre_update (csi *A_col_modified, csi n , const iss *S)
{
    // I_0 = A_col_modified | I_1 = L_col_modified
    csi k, i ;
    csi *parent, *L_col_modified ;
    parent = S->parent ;
    L_col_modified = cs_malloc (n, sizeof(csi)) ;

    for (k = 0 ; k < n ; k++)
    {
        if (A_col_modified == 1)
        {
            L_col_modified [k] = 1 ;
            for (i = parent [k] ; i != -1 ; i = parent [i])
                L_col_modified [k] = 1 ;
        }
    }
    return L_col_modified ;
}

csn *is_left_cholupdate (const cs *A, const iss *S, csn *N, csi *L_col_modified)
{
    cs *L ;
    csi k, n, i, j ;
    csi *A_colptr, *A_rowind, *Lp, *Li, *L_rowptr, *L_colind, *Lpk ;
    double lkj, lkk ;
    double *Ax, *a, *Lx ;

    n = A->n ;
    L = N->L ;
    Lp = L->p ; Li = L->i ; Lx = L->x ;

    for (k = 0 ; k < n ; k++)
    {
        if (Lp [k+1] - Lp [k] > 1)
            Lpk [k] = Lp [k] + 1 ;
        else
            Lpk [k] = -1 ;
    }

    for (k = 0 ; k < n ; k++)
    {
        if (L_col_modified [k] == 1)
        {
            /* a (k:n) = A (k:n,k) */
            for (i = A_colptr [k] ; i < A_colptr [k+1] ; i++)
            {
               a [A_rowind [i]] = Ax [i] ;
            }

            /* for j = find (L (k,;)) */
            for (i = L_rowptr [k] ; i < L_rowptr [k+1] ; i++)
            {
                j = L_colind [i] ;
                lkj = Lx [Lpk [j]] ;
              
                for (int p = Lpk [j] ; p < Lp [j+1] ; p++)
                    a [Li[p]] -= Lx [p]*lkj;
              
                Lpk [j] ++ ;
            }

            /* L (k,k) = sqrt (a (k)) */ 
            Lx [ Lp [k]] = lkk = sqrt (a [Li [Lp [k]]]) ;
            a [Li [Lp [k]]] = 0.0 ;

            /* L (k+1:n,k) = a (k+1:n) / L (k,k) */
            for (j = Lp [k] + 1 ; j < Lp [k+1] ; j++)
            {
                Lx [j] = a [Li [j]] / lkk ;
                a [Li [j]] = 0.0 ;
            }
        }
    }
    return N ;
}