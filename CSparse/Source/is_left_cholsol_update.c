#include "cs.h"
#include <time.h>

static double tic (void) { return (clock () / (double) CLOCKS_PER_SEC) ; }
static double toc (double t) { double s = tic () ; return (CS_MAX (0, s-t)) ; }

/* x=A\b where A is symmetric positive definite; b overwritten with solution */
csi is_left_cholsol_update (csi order, cs *A, double *b, FILE * filePtr)
{
    double *x, t ;
    iss *S ;
    csn *N ;
    cs *C ;
    csi n, ok, *Perm, *pinv ;
    csi *I0, *I1 ;
    csi I0_size, I1_size ;
    if (!CS_CSC (A) || !b) return (0) ;     /* check inputs */
    n = A->n ;

    Perm = cs_amd (order, A) ;     /* P = amd(A+A'), or natural */
    pinv = cs_pinv (Perm, n) ;     /* find inverse permutation */
    cs_free (Perm) ;
    if (order && !pinv) return (0) ;

    C = is_symperm (A, pinv, 1) ;      /* permuting matrix */
    S = is_left_schol (order, C) ;     /* ordering and symbolic analysis */
    N = is_left_chol (C, S) ;          /* numeric Cholesky factorization */
    x = cs_malloc (n, sizeof (double)) ;    /* get workspace */
    print_fill_in (C, N) ;
    //printf ("L = \n") ; cs_print (N->L, 0) ;

    I0_size = I1_size = 0 ;
    I0 = cs_malloc (2, sizeof (csi)) ;
    I0 = is_load_update_matrix (filePtr, A, I0, &I0_size) ;

    t = tic () ;
    I1 = cs_malloc (I0_size, sizeof (csi)) ;
    I1 = is_pre_update2 (I0, I0_size, I1, &I1_size, S, n) ;

    C = is_symperm (A, pinv, 1) ;
    N = is_left_cholupdate (C, S, N, I1, I1_size) ;

    cs_free (I0) ;
    cs_free (I1) ;

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
    printf ("%8.2f s  | ", toc (t)) ;
    return (ok) ;
}