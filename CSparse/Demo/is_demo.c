#include "cs.h"
/* fichier de demo IS pour la compréhension de CSparse */

csi *cs_pinv_identity (csi const *p, csi n)
{
    csi k, *pinv ;
    if (!p) return (NULL) ;                     /* p = NULL denotes identity */
    pinv = cs_malloc (n, sizeof (csi)) ;        /* allocate result */
    if (!pinv) return (NULL) ;                  /* out of memory */
    for (k = 0 ; k < n ; k++) pinv [k] = k ;    /* invert the permutation */
    return (pinv) ;                             /* return result */
}

/* ---------------------------------------------------------------------- */

int main (void)
{
    printf ("=================================================================> IS_DEMO \n") ;

    /* ------------------------------ */
    /* --- Parameters of the demo --- */
    /*                                */
    /*
    int symbolic_test, uplooking_test, leftlooking_test ;

    symbolic_test = 1 ;
    uplooking_test = 1 ;
    leftlooking_test = 1 ;

    cs *T, *A ;
    csi k, order, mn, m, n ;
    double *b ; */

    /* ---------------------------------------------------------------------- */
    /*
    T = cs_load (stdin) ;         
    A = cs_compress (T) ;               
    cs_spfree (T) ;                 


    m = A->m ;
    n = A-> n ;
    mn = CS_MAX (m,n) ;
    b = cs_malloc (mn, sizeof (double)) ;
    for (k = 0 ; k < n ; k++)
        b [k] = 1 + ((double) k) / m ;

    order = 1 ;
    cs_cholsol (order, A, b) ;

    css *S_ ;
    if (symbolic_test == 1)
    {
        csi n, *c, *post, *P ;
        cs *C ;
        css *S ;

        if (!CS_CSC (A)) return (NULL) ;

        n = A->n ;
        S = cs_calloc (1, sizeof (css)) ;       
        if (!S) return (NULL) ;                 
        P = cs_amd (order, A) ;                 
        S->pinv = cs_pinv_identity (P, n) ;             
        cs_free (P) ;
        if (order && !S->pinv) return (cs_sfree (S)) ;
        C = cs_symperm (A, S->pinv, 0) ;      
        S->parent = cs_etree (C, 0) ;       
        post = cs_post (S->parent, n) ;         
        c = cs_counts (C, S->parent, post, 0) ; 
        cs_free (post) ;
        cs_spfree (C) ;
        S->cp = cs_malloc (n+1, sizeof (csi)) ; 
        S->unz = S->lnz = cs_cumsum (S->cp, c, n) ;
        cs_free (c) ;
        S_ = ((S->lnz >= 0) ? S : cs_sfree (S)) ;

    }
    else
    {
        S_ = cs_schol (order, A) ; 
    }*/

    // ---------------------------------------------------------------------- //
    /*
    csn *N_ ; 
    if (uplooking_test == 1)
    {
        double d, lki, *Lx, *x, *Cx ;
        csi top, i, p, k, n, *Li, *Lp, *cp, *pinv, *s, *c, *parent, *Cp, *Ci ;
        cs *L, *C, *E ;
        csn *N ;

        if (!CS_CSC (A) || !S_ || !S_->cp || !S_->parent) 
            return (NULL) ;

        n = A->n ;
        N = cs_calloc (1, sizeof (csn)) ;    
        c = cs_malloc (2*n, sizeof (csi)) ;     
        x = cs_malloc (n, sizeof (double)) ;    
        cp = S_->cp ; pinv = S_->pinv ; parent = S_->parent ;
        C = pinv ? cs_symperm (A, pinv, 1) : ((cs *) A) ;
        E = pinv ? C : NULL ;           
        if (!N || !c || !x || !C) 
            return (cs_ndone (N, E, c, x, 0)) ;
        s = c + n ;
        Cp = C->p ; Ci = C->i ; Cx = C->x ;
        N->L = L = cs_spalloc (n, n, cp [n], 1, 0) ;   
        if (!L) 
            return (cs_ndone (N, E, c, x, 0)) ;
        Lp = L->p ; Li = L->i ; Lx = L->x ;
        for (k = 0 ; k < n ; k++) 
            Lp [k] = c [k] = cp [k] ;

        for (k = 0 ; k < n ; k++)      
        {

            top = cs_ereach (C, k, parent, s, c) ;    
            x [k] = 0 ;                          
            for (p = Cp [k] ; p < Cp [k+1] ; p++)      
            {
                
                if (Ci [p] <= k)
                {
                    x [Ci [p]] = Cx [p] ;
                }
            }
            d = x [k] ;                    
            x [k] = 0 ;                     
            for ( ; top < n ; top++)   
            {
                i = s [top] ;          
                lki = x [i] / Lx [Lp [i]] ; 
                x [i] = 0 ;                

                for (p = Lp [i] + 1 ; p < c [i] ; p++)
                {
                    x [Li [p]] -= Lx [p] * lki ;
                }
                d -= lki * lki ;       
                p = c [i]++ ;
                Li [p] = k ;              
                Lx [p] = lki ;
            }
            if (d <= 0)
            {
                return (cs_ndone (N, E, c, x, 0)) ; 
            }
            p = c [k]++ ;
            Li [p] = k ;             
            Lx [p] = sqrt (d) ;
        }   

        Lp [n] = cp [n] ;               
        N_ = cs_ndone (N, E, c, x, 1) ;
    }
    else
    {
        N_ = cs_chol (A, S_) ; 
    }

    printf ("------------------------------------------- \n" \
            "--------- Let's display the result -------- \n") ;

    cs *Lf ;
    double *Lfx ;
    csi *Lfp, *Lfi ;

    Lf = N_->L ;
    Lfp = Lf->p ; Lfi = Lf->i ; Lfx = Lf->x ;

    printf ("L:\n") ; cs_print (Lf, 0) ; */

    return (0) ;
}