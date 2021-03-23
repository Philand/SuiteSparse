#include "cs.h"
/* fichier de demo IS pour la compréhension de CSparse */

int main (void)
{

    /* ------------------------------ */
    /* --- Parameters of the demo --- */
    /*                                */

    int symbolic_test, uplooking_test ;

    symbolic_test = 0 ;
    uplooking_test = 1 ;

    /*                                */
    /* ------------------------------ */
    /*     display of parameters      */

    printf ("==================================================== \n" \
            "==================================================== \n" \
            "              Beginning of execution                 \n" \
            "==================================================== \n" \
            "==================================================== \n") ;

    printf ("display of parameters : \n") ;
    printf ("symbolic_test = %d \n", symbolic_test) ;
    printf ("uplooking_test = %d \n", uplooking_test) ;

    /*                                */
    /* ------------------------------ */
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

    printf ("Successful execution ! \n") ;

    /* ---------------------------------------------------------------------- */
    // ---------------------------------------------------------------------- //
    /* ---------------------------------------------------------------------- */

    printf ("Decomposition of the function cs_cholsol to analyse the symbolic" \
             "analysis and the up-looking cholesky factorization : \n") ;

    // ---------------------------------------------------------------------- //

    css *S_ ; /* (cs_symbolic) object initialization */
    if (symbolic_test == 1)
    {
        /* --- function cs_schol in details --- */

        /* declaration of parameters */
        csi n, *c, *post, *P ;
        cs *C ;
        css *S ;

        /* check inputs */
        if (!CS_CSC (A)) return (NULL) ;

        /* initialization of parameters */
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
        S_ = ((S->lnz >= 0) ? S : cs_sfree (S)) ;

    }
    else
    {
        /* --- we just apply the cs_schol function to then apply cs_chol --- */

        S_ = cs_schol (order, A) ;    /* ordering and symbolic analysis */
    }

    // ---------------------------------------------------------------------- //

    csn *N_ ; /* (cs_numeric) object initialization */
    if (uplooking_test == 1)
    {

        printf ("================================ \n" \
                "Beginning of the up-looking test \n" \
                "================================ \n") ;

        /* L = chol (A, [pinv parent cp]), pinv is optional */
        /* --- csn *cs_chol (const cs *A, const css *S) --- */

        /* declaration of parameters */
        double d, lki, *Lx, *x, *Cx ;
        csi top, i, p, k, n, *Li, *Lp, *cp, *pinv, *s, *c, *parent, *Cp, *Ci ;
        cs *L, *C, *E ;
        csn *N ;

        /* check if we have the good objects */
        if (!CS_CSC (A) || !S_ || !S_->cp || !S_->parent) 
            return (NULL) ;

        /* initialization of parameters */
        n = A->n ;
        N = cs_calloc (1, sizeof (csn)) ;       /* allocate result */
        c = cs_malloc (2*n, sizeof (csi)) ;     /* get csi workspace */
        x = cs_malloc (n, sizeof (double)) ;    /* get double workspace */
        cp = S_->cp ; pinv = S_->pinv ; parent = S_->parent ;
        C = pinv ? cs_symperm (A, pinv, 1) : ((cs *) A) ;
        E = pinv ? C : NULL ;           /* E is alias for A, or a copy E=A(p,p) */
        if (!N || !c || !x || !C) 
            return (cs_ndone (N, E, c, x, 0)) ;
        s = c + n ;
        Cp = C->p ; Ci = C->i ; Cx = C->x ;
        N->L = L = cs_spalloc (n, n, cp [n], 1, 0) ;    /* allocate result */
        if (!L) 
            return (cs_ndone (N, E, c, x, 0)) ;
        Lp = L->p ; Li = L->i ; Lx = L->x ;
        for (k = 0 ; k < n ; k++) 
            Lp [k] = c [k] = cp [k] ;

        printf ("================================ \n" \
                " Beginning of the main for loop  \n" \
                "================================ \n") ;

        printf("column number : n = %td\n", n, \
               "------------------------------\n") ;

        for (k = 0 ; k < n ; k++)       /* compute L(k,:) for L*L' = C */
        {
            printf ("-------------------------------------------------------" \
                    "k = %td\n", k) ;

            printf ("--- Nonzero pattern of L(%td,:) --- \n", k) ;
            /* --- Nonzero pattern of L(k,:) ------------------------------------ */
            top = cs_ereach (C, k, parent, s, c) ;      /* find pattern of L(k,:) */
            printf ("top = %td\n", top) ;
            x [k] = 0 ;                                 /* x (0:k) is now zero */
            printf ("Start of the intern loop for (p) from Cp[%td] = %td to Cp[%td] = %td \n", k, Cp [k], k, Cp [k+1]) ;
            for (p = Cp [k] ; p < Cp [k+1] ; p++)       /* x = full(triu(C(:,k))) */
            {
                
                if (Ci [p] <= k)
                {
                    x [Ci [p]] = Cx [p] ;
                    printf ("x [Ci [p]] = Cx [p] = %f \n", x [Ci [p]]) ;
                }
            }
            d = x [k] ;                     /* d = C(k,k) */
            printf ("d = %f\n", d, "d = C(k,k)") ;
            x [k] = 0 ;                     /* clear x for k+1st iteration */
            /* --- Triangular solve --------------------------------------------- */
            printf ("--- Triangular solve --- \n") ;
            printf ("Start of the intern loop for (top) from %td to n = %td \n", top, n) ;
            for ( ; top < n ; top++)    /* solve L(0:k-1,0:k-1) * x = C(:,k) */
            {
                i = s [top] ;               /* s [top..n-1] is pattern of L(k,:) */
                printf ("i = %td \n", i) ;
                lki = x [i] / Lx [Lp [i]] ; /* L(k,i) = x (i) / L(i,i) */
                printf ("lki = %f \n", lki) ;
                x [i] = 0 ;                 /* clear x for k+1st iteration */

                printf ("Start of the intern loop for from Lp[%td] = %td to c[%td] = %td \n", i, Lp [i] + 1, i, c [i]) ;
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
            if (d <= 0)
            {
                printf (" matrix not definite positive ! \n") ;
                return (cs_ndone (N, E, c, x, 0)) ; /* not pos def */
            }
            p = c [k]++ ;
            Li [p] = k ;                /* store L(k,k) = sqrt (d) in column k */
            Lx [p] = sqrt (d) ;
        }   
        Lp [n] = cp [n] ;               /* finalize L */
        N_ = cs_ndone (N, E, c, x, 1) ; /* success: free E,s,x; return N */
    }
    else
    {
        /* --- we just apply the cs_chol function --- */

        N_ = cs_chol (A, S_) ;  /* numeric Cholesky factorization ( up-looking ) */
    }


    return (0) ;
}
