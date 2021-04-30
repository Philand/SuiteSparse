#include "cs.h"

/* --- is_left_cholsol.c --- */

csi *cs_pinv_identity (csi const *p, csi n)
{
    csi k, *pinv ;
    if (!p) return (NULL) ;                     /* p = NULL denotes identity */
    pinv = cs_malloc (n, sizeof (csi)) ;        /* allocate result */
    if (!pinv) return (NULL) ;                  /* out of memory */
    for (k = 0 ; k < n ; k++) pinv [k] = k ;    /* invert the permutation */
    return (pinv) ;                             /* return result */
}

/* quicksort sur Ci */
/* vecteur dense pour stocker les valeurs de Cx */
/* puis on les stocke dans Cx selon l'ordre d'indices de Ci */

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

/* C = A(p,q) where p and q are permutations of 0..m-1 and 0..n-1. */
cs *is_permute_so (const cs *A, const csi *pinv, csi values)
{
    csi t, j, k, nz = 0, m, n, *Ap, *Ai, *Cp, *Ci ;
    double *Cx, *Ax, *a ;
    cs *C ;
    if (!CS_CSC (A)) return (NULL) ;    /* check inputs */
    m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    C = cs_spalloc (m, n, Ap [n], values && Ax != NULL, 0) ;  /* alloc result */
    if (!C) return (cs_done (C, NULL, NULL, 0)) ;   /* out of memory */
    a = cs_malloc (n, sizeof (double)) ;
    for (k = 0 ; k < n ; k++)
        a [k] = 0.0 ;
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    for (k = 0 ; k < n ; k++)
    {
        Cp [k] = nz ;
        j = k ;
        for (t = Ap [j] ; t < Ap [j+1] ; t++)
        {
            a [pinv [Ai [t]]] =  Ax [t] ;
            Ci [nz++] = pinv ? (pinv [Ai [t]]) : Ai [t] ;
        }
        qsort (&Ci [Cp [k]], nz - Cp [k], sizeof (csi), csiComparator) ;
        is_build_Cx_column (&Cx [Cp [k]], a, &Ci [Cp [k]], nz - Cp [k]) ; 
    }
    Cp [n] = nz ;                       /* finalize the last column of C */
    return (cs_done (C, NULL, NULL, 1)) ;
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

    Perm = cs_amd (order, A) ;     /* P = amd(A+A'), or natural */
    pinv = cs_pinv (Perm, n) ;     /* find inverse permutation */
    cs_free (Perm) ;
    if (order && !pinv) return (0) ;
    if (order == 0)
        A_perm = cs_permute (A, pinv, NULL, 1) ;
    else if (order == 1)
    {
        A_perm = is_permute_so (A, pinv, 1) ;
        printf (" MATRIX A : \n") ;
        cs_print (A, 0) ;
        printf ("pinv = [ ") ;
        for (csi k = 0 ; k < n ; k++)
            printf ("%td ", pinv [k]) ;
        printf ("]\n") ;
        printf (" MATRIX A permuted : \n") ;
        cs_print (A_perm, 0) ;
    }

    S = is_left_schol (order, A_perm) ;     /* ordering and symbolic analysis */
    N = is_left_chol (A_perm, S) ;          /* numeric Cholesky factorization */
    printf ("L:\n") ; cs_print (N->L, 0) ;
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
    cs_spfree (A_perm) ;
    return (ok) ;
}