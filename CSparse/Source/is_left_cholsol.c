#include "cs.h"

/* --- is_left_cholsol.c --- */

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

int csiComparator ( const void * first, const void * second)
{
    csi firstCsi = * (const csi *) first ;
    csi secondCsi = * (const csi *) second ;
    return (int) firstCsi - secondCsi ;
}

/* build of a column, returning the numbers of elements */
void is_build_Cx_column (double * Cx, double *a, csi * Ci, csi nb_nz_col)
{
    csi k ;
    for (k = 0 ; k < nb_nz_col ; k++)
    {
        Cx [k] = a [Ci [k]] ;
        a [Ci [k]] = 0.0 ;
    }
}

cs *order_column (cs *C)
{
    csi k, j, n, nz, *Cp, *Ci ;
    double *a, *Cx ;
    n = C->n ; Cp = C->p ; Ci = C->i ; Cx = C->x ;
    a = cs_malloc (n, sizeof (double)) ;
    for (k = 0 ; k < n ; k++)
    {
        for (j = Cp [k] ; j < Cp [k+1] ; j++)
            a [Ci [j]] = Cx [j] ;
        nz = Cp [k+1] - Cp [k] ;
        qsort (&Ci [Cp [k]], nz, sizeof (csi), csiComparator) ;
        is_build_Cx_column (&Cx [Cp [k]], a, &Ci [Cp [k]], nz) ; 
    }
    free (a) ;
    return (cs_done (C, NULL, NULL, 1)) ;
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

/* x=A\b where A is symmetric positive definite; b overwritten with solution */
csi is_left_cholsol (csi order, const cs *A, double *b)
{
    double *x ;
    iss *S ;
    csn *N ;
    cs *C ;
    csi n, ok, *Perm, *pinv ;
    if (!CS_CSC (A) || !b) return (0) ;     /* check inputs */
    n = A->n ;

    Perm = cs_amd (order, A) ;     /* P = amd(A+A'), or natural */
    pinv = cs_pinv (Perm, n) ;     /* find inverse permutation */
    cs_free (Perm) ;
    if (order && !pinv) return (0) ;

    C = is_symperm (A, pinv, 1) ;
    //C = order_column (C) ;

    S = is_left_schol (order, C) ;     /* ordering and symbolic analysis */
    N = is_left_chol (C, S) ;          /* numeric Cholesky factorization */
    x = cs_malloc (n, sizeof (double)) ;    /* get workspace */
    print_fill_in (C, N) ;
    ok = (S && N && x) ;
    if (ok)
    {
        cs_ipvec (pinv, b, x, n) ;   /* x = P*b */
        cs_lsolve (N->L, x) ;        /* x = L\x */
        cs_ltsolve (N->L, x) ;       /* x = L'\x */
        cs_pvec (pinv, x, b, n) ;    /* b = P'*x */
    }
    cs_free (x) ;
    cs_free (pinv) ;
    is_sfree (S) ;
    cs_nfree (N) ;
    cs_spfree (C) ;
    return (ok) ;
}