#include "cs.h"
#include "stdbool.h"
/* fichier de demo IS pour la compréhension de CSparse */

int main (void)
{

    /* ------------------------------ */
    /* --- Parameters of the demo --- */
    /*                                */

    bool symbolic_test, uplooking_test ;

    symbolic_test = false ;
    uplooking_test = true ;

    /*                                */
    /* ------------------------------ */

    /* ---------------------------------------------------------------------- */

    cs *T, *A ;
    csi k, order, mn, m, n ;
    double *b ;

    /* ---------------------------------------------------------------------- */
    
    T = cs_load (stdin) ;               /* load triplet matrix T from stdin */
    printf ("T:\n") ; cs_print (T, 0) ; /* print T */
    A = cs_compress (T) ;               /* A = compressed-column form of T */
    printf ("A:\n") ; cs_print (A, 0) ; /* print A */
    cs_spfree (T) ;                     /* clear T */

    /* ---------------------------------------------------------------------- */

    printf ("Execution of the cs_cholsol function on the matrix A : \n") ;

    m = A->m ;
    n = A-> n ;
    mn = CS_MAX (m,n) ;
    b = cs_malloc (mn, sizeof (double)) ;
    for (k = 0 ; k < n ; k++)
        b [k] = 1 + ((double) k) / m ;

    order = 1 ;
    cs_cholsol (order, A, b) ; /* solve Ax=b with Cholesky */

    printf ("Successful execution !") ;

    /* ---------------------------------------------------------------------- */
    // ---------------------------------------------------------------------- //
    /* ---------------------------------------------------------------------- */

    printf ("Decomposition of the function cs_cholsol to analyse the symbolic \
             analysis and the up-looking cholesky factorization : \n") ;

    // ---------------------------------------------------------------------- //

    css *S ; /* (cs_symbolic) object initialization */
    if (symbolic_test == true)
    {
        /* --- function cs_schol in details --- */

        csi n, *c, *post, *P ;
        cs *C ;
        css *S ;
        if (!CS_CSC (A)) return (NULL) ;        /* check inputs */
        n = A->n ;
        S = cs_calloc (1, sizeof (css)) ;       /* allocate result S */
        if (!S) return (NULL) ;                 /* out of memory */
        P = cs_amd (order, A) ;                 /* P = amd(A+A'), or natural */
        S->pinv = cs_pinv (P, n) ;              /* find inverse permutation */
        cs_free (P) ;
        if (order && !S->pinv) return (cs_sfree (S)) ;
        C = cs_symperm (A, S->pinv, 0) ;        /* C = spones(triu(A(P,P))) */
        S->parent = cs_etree (C, 0) ;           /* find etree of C */
        post = cs_post (S->parent, n) ;         /* postorder the etree */
        c = cs_counts (C, S->parent, post, 0) ; /* find column counts of chol(C) */
        cs_free (post) ;
        cs_spfree (C) ;
        S->cp = cs_malloc (n+1, sizeof (csi)) ; /* allocate result S->cp */
        S->unz = S->lnz = cs_cumsum (S->cp, c, n) ; /* find column pointers for L */
        cs_free (c) ;
        return ((S->lnz >= 0) ? S : cs_sfree (S)) ;

    }
    else
    {
        /* --- we just apply the cs_schol function to then apply cs_chol --- */

        S = cs_schol (order, A) ;    /* ordering and symbolic analysis */
    }

    // ---------------------------------------------------------------------- //

    csn *N ; /* (cs_numeric) object initialization */
    if (uplooking_test == true)
    {
        /* L = chol (A, [pinv parent cp]), pinv is optional */
        /* --- csn *cs_chol (const cs *A, const css *S) --- */

        double d, lki, *Lx, *x, *Cx ;
        csi top, i, p, k, n, *Li, *Lp, *cp, *pinv, *s, *c, *parent, *Cp, *Ci ;
        cs *L, *C, *E ;
        csn *N ;
        if (!CS_CSC (A) || !S || !S->cp || !S->parent) return (NULL) ;
        n = A->n ;
        N = cs_calloc (1, sizeof (csn)) ;       /* allocate result */
        c = cs_malloc (2*n, sizeof (csi)) ;     /* get csi workspace */
        x = cs_malloc (n, sizeof (double)) ;    /* get double workspace */
        cp = S->cp ; pinv = S->pinv ; parent = S->parent ;
        C = pinv ? cs_symperm (A, pinv, 1) : ((cs *) A) ;
        E = pinv ? C : NULL ;           /* E is alias for A, or a copy E=A(p,p) */
        if (!N || !c || !x || !C) return (cs_ndone (N, E, c, x, 0)) ;
        s = c + n ;
        Cp = C->p ; Ci = C->i ; Cx = C->x ;
        N->L = L = cs_spalloc (n, n, cp [n], 1, 0) ;    /* allocate result */
        if (!L) return (cs_ndone (N, E, c, x, 0)) ;
        Lp = L->p ; Li = L->i ; Lx = L->x ;
        for (k = 0 ; k < n ; k++) Lp [k] = c [k] = cp [k] ;
        for (k = 0 ; k < n ; k++)       /* compute L(k,:) for L*L' = C */
        {
            /* --- Nonzero pattern of L(k,:) ------------------------------------ */
            top = cs_ereach (C, k, parent, s, c) ;      /* find pattern of L(k,:) */
            x [k] = 0 ;                                 /* x (0:k) is now zero */
            for (p = Cp [k] ; p < Cp [k+1] ; p++)       /* x = full(triu(C(:,k))) */
            {
                if (Ci [p] <= k) x [Ci [p]] = Cx [p] ;
            }
            d = x [k] ;                     /* d = C(k,k) */
            x [k] = 0 ;                     /* clear x for k+1st iteration */
            /* --- Triangular solve --------------------------------------------- */
            for ( ; top < n ; top++)    /* solve L(0:k-1,0:k-1) * x = C(:,k) */
            {
                i = s [top] ;               /* s [top..n-1] is pattern of L(k,:) */
                lki = x [i] / Lx [Lp [i]] ; /* L(k,i) = x (i) / L(i,i) */
                x [i] = 0 ;                 /* clear x for k+1st iteration */
                for (p = Lp [i] + 1 ; p < c [i] ; p++)
                {
                    x [Li [p]] -= Lx [p] * lki ;
                }
                d -= lki * lki ;            /* d = d - L(k,i)*L(k,i) */
                p = c [i]++ ;
                Li [p] = k ;                /* store L(k,i) in column i */
                Lx [p] = lki ;
            }
            /* --- Compute L(k,k) ----------------------------------------------- */
            if (d <= 0) return (cs_ndone (N, E, c, x, 0)) ; /* not pos def */
            p = c [k]++ ;
            Li [p] = k ;                /* store L(k,k) = sqrt (d) in column k */
            Lx [p] = sqrt (d) ;
        }   
        Lp [n] = cp [n] ;               /* finalize L */
        N = cs_ndone (N, E, c, x, 1) ; /* success: free E,s,x; return N */
    }
    else
    {
        /* --- we just apply the cs_chol function --- */

        N = cs_chol (A, S) ;  /* numeric Cholesky factorization ( up-looking ) */
    }


    return (0) ;
}
