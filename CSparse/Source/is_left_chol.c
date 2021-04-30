#include "cs.h"

/* --- is_left_chol.c --- */

csn *is_left_chol (const cs *A, const iss *S)
{
    if (!CS_CSC (A) || !S) return (NULL) ;

    // printf ("\n \n") ;
    // cs_print (A, 0) ;

    csn *N ;
    cs *L ;
    csi *Lp, *Li, *Lpk ;
    csi *L_colind, *L_colptr, *L_rowind, *L_rowptr ;
    csi *A_colptr, *A_rowind ;
    csi n, k, i, j ;
    double *a, *Ax, *Lx, lkj, lkk ;

    L_colind = S->L_colind ;
    L_colptr = S->L_colptr ;
    L_rowind = S->L_rowind ;
    L_rowptr =  S->L_rowptr ;
    n = A->n ;
    N = cs_calloc (1, sizeof (csn)) ;                   /* allocate result */
    a = cs_malloc(n, sizeof (double)) ;                 /* get csi workspace */
    Ax = A->x ;
    A_colptr = A->p ;
    A_rowind = A->i ;
    N->L = L = cs_spalloc (n, n, L_colptr [n], 1, 0) ;  /* allocate result */
    Lp = L->p ; Li = L->i ; Lx = L->x ;
    Lpk = cs_malloc(n, sizeof (csi)) ;                  /* current pointer for column k of L */

    /* on connait le nombre de valeurs dans chaque colonne */
    for (k = 0 ; k < n+1 ; k++) 
        Lp [k] = L_colptr [k] ;

    for (k = 0 ; k < n ; k++)
    {
        // const int nnz = Lp[k+1]-Lp[k]
        // assert(nnz >= 1)
        // if( nnz > 1) Lpk[k] = 1; /// (Lp[k] + 1) 
        if (Lp [k+1] - Lp [k] > 1)
            Lpk [k] = Lp [k] + 1 ;
        else
            Lpk [k] = -1 ;
    }

    /*=====================================================*/
    /* algorithme (boucle calculant L colonne par colonne) */
    /*=====================================================*/
    
    for (k = 0 ; k < n ; k++)
    {
        printf ("\n--------------------------------------------- \n") ;
        printf ("--------------------------------------  k = %td \n \n", k) ;

        /* a (k:n) = A (k:n,k) */
        for (i = A_colptr [k] ; i < A_colptr [k+1] ; i++)
        {
            a [A_rowind [i]] = Ax [i] ;
        }

        printf ("1) Récupérons les valeurs de Ax dans a : \n \n") ;
        printf ("a = [ ") ;
        for (j = 0 ; j < n ; j++)
            printf ("%f ", a [j]) ;
        printf ("]\n") ;

        /* for j = find (L (k,;)) */
        for (i = L_rowptr [k] ; i < L_rowptr [k+1] ; i++)
        {
            j = L_colind [i] ;
            lkj = Lx [Lpk [j]] ;
              
            for (int p = Lpk [j] ; p < Lp [j+1] ; p++)
                  a [L_rowind[p]] -= Lx [p]*lkj;
              
            Lpk [j] ++ ;
        }

        printf ("\n2) Valeurs de a après substitutions des colonnes à gauche : \n \n") ;
        printf ("a = [ ") ;
        for (j = 0 ; j < n ; j++)
            printf ("%f ", a [j]) ;
        printf ("]\n") ;

        /* L (k,k) = sqrt (a (k)) */ 
        Lx [ Lp [k]] = lkk = sqrt (a [L_rowind [Lp [k]]]) ;
        // printf ("lkk = %f \n", lkk) ;
        a [L_rowind [Lp [k]]] = 0.0 ;
        Li [ Lp [k]] = L_rowind [Lp [k]];

        /* L (k+1:n,k) = a (k+1:n) / L (k,k) */
        for (j = Lp [k] + 1 ; j < Lp [k+1] ; j++)
        {
            Lx [j] = a [L_rowind [j]] / lkk ;
            a [L_rowind [j]] = 0.0 ;
            Li [j] = L_rowind [j];
        }

        printf ("\n3) Après refresh : \n \n") ;
        printf ("a = [ ") ;
        for (j = 0 ; j < n ; j++)
            printf ("%f ", a [j]) ;
        printf ("]\n") ;
    }

    cs_free (Lpk) ;
    cs_free (a) ;

    return N ;
}
