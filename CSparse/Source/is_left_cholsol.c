#include "cs.h"

/* --- is_left_cholsol.c --- */

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

/* C = A(p,p) where A and C are symmetric the upper part stored; pinv not p */
cs *is_symperm (const cs *A, const csi *pinv, csi values)
{
    csi i, j, p, q, i2, j2, n, *Ap, *Ai, *Cp, *Ci, *w ;
    double *Cx, *Ax ;
    cs *C ;
    if (!CS_CSC (A)) return (NULL) ;                    /* check inputs */
    n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    C = cs_spalloc (n, n, Ap [n], values && (Ax != NULL), 0) ; /* alloc result*/
    w = cs_calloc (n, sizeof (csi)) ;                   /* get workspace */
    if (!C || !w) return (cs_done (C, w, NULL, 0)) ;    /* out of memory */
    Cp = C->p ; Ci = C->i ; Cx = C->x ;

    for (j = 0 ; j < n ; j++)           /* count entries in each column of C */
    {
        j2 = pinv ? pinv [j] : j ;      /* column j of A is column j2 of C */
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            i = Ai [p] ; 
            i2 = pinv ? pinv [i] : i ;  /* row i of A is row i2 of C */
            if (i > j)
                w [CS_MIN (i2, j2)]++ ; /* column count of C */
            else
                w [CS_MAX (i2, j2)]++ ; /* column count of C */
        }
    }
    cs_cumsum (Cp, w, n) ;              /* compute column pointers of C */
    for (j = 0 ; j < n ; j++)
    {
        j2 = pinv ? pinv [j] : j ;      /* column j of A is column j2 of C */
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            i = Ai [p] ;
            i2 = pinv ? pinv [i] : i ;  /* row i of A is row i2 of C */
            if (i > j)
                Ci [q = w [CS_MIN (i2, j2)]++] = CS_MAX (i2, j2) ;
            else
                Ci [q = w [CS_MAX (i2, j2)]++] = CS_MIN (i2, j2) ;
            if (Cx) Cx [q] = Ax [p] ;
        }
    }
    return (cs_done (C, w, NULL, 1)) ;  /* success; free workspace, return C */
}

csi *cs_pinv_identity (csi const *p, csi n)
{
    csi k, *pinv ;
    if (!p) return (NULL) ;                     /* p = NULL denotes identity */
    pinv = cs_malloc (n, sizeof (csi)) ;        /* allocate result */
    if (!pinv) return (NULL) ;                  /* out of memory */
    for (k = 0 ; k < n ; k++) pinv [k] = k ;    /* invert the permutation */
    return (pinv) ;                             /* return result */
}

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
    C = order_column (C) ;

    S = is_left_schol (order, C) ;     /* ordering and symbolic analysis */
    N = is_left_chol (C, S) ;          /* numeric Cholesky factorization */
    x = cs_malloc (n, sizeof (double)) ;    /* get workspace */
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