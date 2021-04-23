#include "cs.h"
/* fichier de demo 3 (InSimo) pour la compréhension de CSparse */
/* deuxième implémentation de la left-looking                  */

void is_bool_union (csi *P, csi *col_rowind, csi size, csi i)
{
    csi row ;
    for (csi k = 0 ; k < size ; k++)
    {   
        row = col_rowind [k] ;
        if (P [row] == 0 && row != i)
            P [row] = 1 ;
    }
}

// transformer en int --> renvoyer le nombre de non zéros de la colonne
csi is_build_column (csi * rowind, csi *P, csi n, csi k)
{
    csi j = 0 ;
    csi nb_nz_col = 0 ;
    for (csi i = k ; i < n ; i++)
    {
        if (P [i] == 1)
        {
            rowind [j++] = i ;
            P [i] = 0 ;
            nb_nz_col ++ ;
        }
    }
    return nb_nz_col ;
}

/* fonction is_write */
void is_write (csi *L_colind, csi *stack, csi top, csi n)
{
    csi imax = n - top ;
    for (csi i = 0 ; i < imax  ; i++)
        L_colind [i] = stack [i] ;
}

int main (void)
{
    printf ("=================================================================> IS_DEMO 3 \n") ;

    /* --- Récupération de la matrice --------------------------------------- */
    cs *T, *A ;
    T = cs_load (stdin) ;               /* load triplet matrix T from stdin */
    A = cs_compress (T) ;               /* A = compressed-column form of T */
    cs_spfree (T) ;                     /* clear T */

    /* ###################################################################### */

    /* ---------------------------------------------------------------------- */
    /* --- Analyse Symbolique de la matrice A ------------------------------- */
    /* ---------------------------------------------------------------------- */

    printf (" \n ------- Analyse Symbolique de la matrice A ------- \n") ;

    /*===========================*/
    /* déclaration des variables */
    /*===========================*/

    csi k, i, j, p; /* variables incrémentales */
    csi len, top ;
    csi n ;         /* taille de la matrice A */
    csi *A_colptr ; /* tab qui compte le nb d'élts par colonne de A */
    csi *A_rowind ; /* tab qui stocke les ind des lignes pour chaque col de A */

    csi *L_colptr ; /* tab qui compte le nb d'élts par colonne de L */
    csi *L_rowind ; /* tab qui stocke les ind des lignes pour chaque col de L */

    csi *L_rowptr ; /* tab qui compte le nb d'élts par ligne de L */
    csi *L_colind ; /* tab qui stocke les ind des col pour chaque ligne e L*/

    csi *parent ;   /* etree */
    csi *flag ;     /* tab pour marquer les noeuds visités (pour 1 itération) */
    csi *stack ;    /* pile pour ajouter les nonzeros de l'itération k */
    csi *P ;        /* boolean array pour l'opérateur d'union */

    csi nb_nz_col ;

    /*============================================*/
    /* initialisation des variables (allocations) */
    /*============================================*/

    n = A->n ; A_colptr = A->p ; A_rowind = A->i ;
    L_colptr = cs_malloc (n+1, sizeof(csi)) ;
    L_colptr [0] = 0 ;
    L_rowptr = cs_malloc (n+1, sizeof(csi)) ;
    L_rowptr [0] = 0 ;
    L_rowind = cs_malloc(A_colptr [n], sizeof(csi)) ;
    L_colind = cs_malloc(A_colptr [n], sizeof(csi)) ;
    parent = cs_malloc (n, sizeof(csi)) ;
    flag = cs_malloc (n, sizeof(csi)) ;
    stack = cs_malloc (n, sizeof(csi)) ;
    P = cs_malloc (n, sizeof(csi)) ;
    for (k = 0 ; k < n ; k++)
        P [k] = 0 ;

    /*==========================================================*/
    /* boucle calculant la structure des lignes et des colonnes */
    /*==========================================================*/

    for (k = 0 ; k < n ; k++)
    {
        parent [k] = -1 ;
        flag [k] = k ;
        top = n ;

        /* structure de la ligne k */
        for (p = A_colptr [k] ; p < A_colptr [k+1] ; p++)
        {
            i = A_rowind[p] ;

            if (i < k)
            {
                len = 0 ;

                for ( ; flag [i] != k ; i = parent [i])
                {
                    if (parent [i] == -1)
                        parent [i] = k ;
                    flag [i] = k ;
                    stack [len++] = i ;
                }

                while (len > 0)
                    stack [--top] = stack [--len] ;
            }    
        }
        
        is_write (&L_colind [L_rowptr [k]], &stack [top], top, n) ;
        L_rowptr [k+1] = L_rowptr [k] + n - top ;

        /* structure de la colonne k */
        
        /* L_k = A_k */
        for (i = A_colptr [k] ; i < A_colptr [k+1] ; i++)
            P [A_rowind [i]] = 1 ;

        /* for all i such that pi [k] = i */
        for (i = L_rowptr [k] ; i < L_rowptr [k+1] ; i++)
        {
            j = L_colind [i] ;
            if (parent [j] == k)
                is_bool_union(&P [0], &L_rowind [L_colptr [j]], L_colptr [j+1] - L_colptr [j], j) ;
        }
        nb_nz_col = is_build_column (&L_rowind [L_colptr [k]], &P [0], n, k) ;
        L_colptr [k+1] = L_colptr [k] + nb_nz_col ;                        /* upd nb_nz */
    }

    /* ---------------------------------------------------------------------- */
    /* --- Vérification de la phase symbolique ------------------------------ */
    /* ---------------------------------------------------------------------- */
    printf (" \n --- Vérification de l'indice pointeur des colonnes : --- \n") ;
    for (k = 0 ; k < n + 1 ; k++)
    {
        printf ("%td  ;  ", L_colptr [k]) ;
    }
    printf ("\n") ;

    printf (" \n --- Vérification de la structure des non-zéros : --- \n") ;
    for (k = 0 ; k < L_colptr [n] ; k++)
    {
        printf ("%td  ;  ", L_rowind [k]) ;
    }
    printf ("\n") ;

        printf (" \n --- Vérification de l'indice pointeur des lignes : --- \n") ;
    for (k = 0 ; k < n + 1 ; k++)
    {
        printf ("%td  ;  ", L_rowptr [k]) ;
    }
    printf ("\n") ;

    printf (" \n --- Vérification de la structure des non-zéros : --- \n") ;
    for (k = 0 ; k < L_rowptr [n] ; k++)
    {
        printf ("%td  ;  ", L_colind [k]) ;
    }
    printf ("\n") ;

    printf ("\n --- Tableau de parent à la fin de l'exécution : --- \n") ;
    for (k = 0 ; k < n ; k++)
    {
        printf ("%td  ;  ", parent [k]) ;
    }
    printf ("\n") ;

    /* ###################################################################### */

    /* ---------------------------------------------------------------------- */
    /* --- Factorisation Numérique de la matrice A -------------------------- */
    /* ---------------------------------------------------------------------- */

    printf (" \n ------- Factorisation Numérique de la matrice A ------- \n") ;

    /*===========================*/
    /* déclaration des variables */
    /*===========================*/
    
    csn *N ;        /* (cs_numeric) object initialization */
    cs *L ;
    csi *Lp, *Li, *Lpk ;
    double *a, *Ax, *Lx, lkj, lkk ;

    /*============================================*/
    /* initialisation des variables (allocations) */
    /*============================================*/

    N = cs_calloc (1, sizeof (csn)) ;                   /* allocate result */
    a = cs_malloc(n, sizeof (double)) ;                 /* get csi workspace */
    Ax = A->x ;
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
                  a [L_rowind[p]] -= Lx [p]*lkj;
              
            Lpk [j] ++ ;
        }

        /* L (k,k) = sqrt (a (k)) */
        Lx [ Lp [k]] = lkk = sqrt (a [L_rowind [Lp [k]]]) ;
        a [L_rowind [Lp [k]]] = 0.0 ;
        Li [ Lp [k]] = L_rowind [Lp [k]];

        /* L (k+1:n,k) = a (k+1:n) / L (k,k) */
        for (j = Lp [k] + 1 ; j < Lp [k+1] ; j++)
        {
            Lx [j] = a [L_rowind [j]] / lkk ;
            a [L_rowind [j]] = 0.0 ;
            Li [j] = L_rowind [j];
        }
    }

    /* --------------------------------------------------------------------- */
    /* --- Vérification de la phase numérique ------------------------------ */
    /* --------------------------------------------------------------------- */
    printf ("L:\n") ; cs_print (L, 0) ;

    /* --------------------------------------------------------------------- */
    /* --- Free ------------------------------------------------------------ */
    /* --------------------------------------------------------------------- */
    cs_free (L_rowptr) ;
    cs_free (L_colptr) ;
    cs_free (L_rowind) ;
    cs_free (L_colind) ;
    cs_free (parent) ;
    cs_free (stack) ;
    cs_free (P) ;
    cs_free (flag) ;
    cs_free (Lpk) ;
    cs_free (a) ;
    cs_spfree (A) ;
    cs_nfree (N) ;

    return 0 ;
}