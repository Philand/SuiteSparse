#include "cs.h"

csi *cs_pinv_identity (csi const *p, csi n)
{
    csi k, *pinv ;
    if (!p) return (NULL) ;                     /* p = NULL denotes identity */
    pinv = cs_malloc (n, sizeof (csi)) ;        /* allocate result */
    if (!pinv) return (NULL) ;                  /* out of memory */
    for (k = 0 ; k < n ; k++) pinv [k] = k ;    /* invert the permutation */
    return (pinv) ;                             /* return result */
}

/* x=A\b where A is symmetric positive definite; b overwritten with solution */
csi is_left_cholsol (csi order, const cs *A, double *b)
{
    double *x ;
    iss *S ;
    csn *N ;
    cs *A_perm ;
    csi n, ok, *Perm, *pinv ;
    if (!CS_CSC (A) || !b) return (0) ;     /* check inputs */
    n = A->n ;

    Perm = cs_amd (order, A) ;              /* P = amd(A+A'), or natural */
    pinv = cs_pinv (Perm, n) ;     /* find inverse permutation */
    cs_free (Perm) ;
    // if (order && !pinv) return (is_sfree (S)) ;
    // A_perm = cs_symperm (A, pinv, 0) ;

    S = is_left_schol (order, A) ;     /* ordering and symbolic analysis */
    N = is_left_chol (A, S) ;          /* numeric Cholesky factorization */
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
    return (ok) ;
}
