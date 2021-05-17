#include "cs.h"

/* --- is_left_cholupdate */

csi *is_pre_update (csi *I0, csi I0_size, csi *I1, csi *I1_size, const iss *S)
{
    csi k, i, j, count, I1_max_size, ok ;
    csi *parent ;
    parent = S->parent ;
    I1_max_size = I0_size ;
    ok = 1 ;

    i = count = 0 ;
    for (k = I0 [i] ; i < I0_size ; k = I0 [i++] )
    {
        for (j = k ; j != -1 ; j = parent [j])
        {
            I1 [count++] = j ;
            if (count == I1_max_size)
            {
                I1_max_size *= 2 ;
                I1 = realloc (I1, I1_max_size * sizeof (csi)) ;
            }
        }
    }
    *I1_size = count ;

    return I1 ;
}


csn *is_left_cholupdate (const cs *A, const iss *S, csn *N, csi *I1, csi I1_size)
{
    cs *L ;
    csi k, n, i, j, o ;
    csi *A_colptr, *A_rowind, *Lp, *Li, *L_rowptr, *L_colind, *Lpk ;
    double lkj, lkk ;
    double *Ax, *a, *Lx ;

    n = A->n ;
    L = N->L ;
    Lp = L->p ; Li = L->i ; Lx = L->x ;
    o = 0 ;

    L_colind = S->L_colind ;
    L_rowptr =  S->L_rowptr ;
    a = cs_malloc(n, sizeof (double)) ;                 /* get csi workspace */
    Ax = A->x ;
    A_colptr = A->p ;
    A_rowind = A->i ;
    Lpk = cs_malloc(n, sizeof (csi)) ;             

    for (k = 0 ; k < n ; k++)
    {
        if (Lp [k+1] - Lp [k] > 1)
            Lpk [k] = Lp [k] + 1 ;
        else
            Lpk [k] = -1 ;
    }

    for (k = I1 [o] ; o < I1_size ; k = I1 [o++] )
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

    cs_free (Lpk) ;
    cs_free (a) ;

    return N ;
}