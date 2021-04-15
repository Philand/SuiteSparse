#include "cs.h"
/* fichier de demo 3 (InSimo) pour la compréhension de CSparse */
/* deuxième implémentation de la left-looking                  */

csi is_bool_union (csi *P, csi *col_rowind, csi size, csi i)
{
    csi row, nb_new ;
    for (csi k = 0 ; k < size ; k++)
    {   
        row = col_rowind [k] ;
        if (P [row] == 0 && row != i)
        {
            P [row] = 1 ;
            nb_new ++ ;
        }
    }
    return nb_new ;
}

void is_build_column (csi * rowind, csi *P, csi n)
{
    csi k = 0 ;
    for (csi i = 0 ; i < n ; i++)
    {
        if (P [i] == 1)
        {
            rowind [k++] = i ;
            P [i] = 0 ;
        }
    }
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

    csi count1, count2, nb_new, start ;
    csi nb_row_col, nb_sup_values ;

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

    /*==========================================================*/
    /* boucle calculant la structure des lignes et des colonnes */
    /*==========================================================*/

    L_rowptr [0] = count1 = 0 ;
    for (k = 0 ; k < n ; k++)
    {
        parent [k] = -1 ;
        L_rowptr [k] = count1 ;
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

        count1 += n - top ;
        is_write (&L_colind [L_rowptr [k]], &stack [top], top, n) ;
        L_rowptr [k+1] = count1 ;

        /* structure de la colonne k */
        
        /* L_k = A_k */
        count2 = nb_row_col = A_colptr [k+1] - A_colptr [k] ;
        start = A_colptr [k] ;
        nb_sup_values = 0 ;

        for (i = 0 ; i < nb_row_col ; i++)
        {
            /* on filtre les valeurs de la triangulaire inférieure de A */
            if (A_rowind[start + i] >= k)
                P [A_rowind [start +i]] = 1 ;
            else
                count2 -- ;
        }

        /* for all i such that pi [k] = i */
        for (i = L_rowptr [k] ; i < L_rowptr [k+1] ; i++)
        {
            j = L_colind [i] ;
            if (parent [j] == k)
                count2 += is_bool_union(&P [0], &L_rowind [L_colptr [j]], L_colptr [j+1] - L_colptr [j], j) ;
        }
        is_build_column (&L_rowind [L_colptr [k]], &P [0], n) ;
        L_colptr [k+1] = L_colptr [k] + count2 ;                        /* upd nb_nz */
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
    csi *Lp, *Li, q, *count_update ;
    double *a, *Ax, *Lx, lkj, lkk ;

    /*============================================*/
    /* initialisation des variables (allocations) */
    /*============================================*/

    N = cs_calloc (1, sizeof (csn)) ;                   /* allocate result */
    a = cs_malloc(n, sizeof (double)) ;                 /* get csi workspace */
    Ax = A->x ;
    N->L = L = cs_spalloc (n, n, L_colptr [n], 1, 0) ;    /* allocate result */
    Lp = L->p ; Li = L->i ; Lx = L->x ;
    count_update = cs_malloc(n, sizeof (csi)) ;

    /* on connait le nombre de valeurs dans chaque colonne */
    for (k = 0 ; k < n+1 ; k++) 
        Lp [k] = L_colptr [k] ;

    /* on initialise les compteurs d'accès aux colonnes */
    for (k = 0 ; k < n ; k++)
    {
        if (Lp [k+1] - Lp [k] > 1)
            count_update [k] = 1 ;
        else
            count_update [k] = 0 ;
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

            p = count_update [j] ; 
            q = 0 ;
            lkj = Lx [Lp [j] + p] ;
              
            for(int q = Lp [j] + p ; q < Lp [j+1] ; q++)
                  a [L_rowind[q]] -= Lx [q]*lkj;
              
            count_update [j] ++ ;
        }

        /* L (k,k) = sqrt (a (k)) */
        Lx [ Lp [k]] = lkk = sqrt (a [L_rowind [Lp [k]]]) ;
        a [L_rowind [Lp [k]]] = 0 ;
        Li [ Lp [k]] = L_rowind [Lp [k]];

        /* L (k+1:n,k) = a (k+1:n) / L (k,k) */
        for (j = Lp [k] + 1 ; j < Lp [k+1] ; j++)
        {
            Lx [j] = a [L_rowind [j]] / lkk ;
            a [L_rowind [j]] = 0 ;
            Li [j] = L_rowind [j];
        }
    }

    /* --------------------------------------------------------------------- */
    /* --- Vérification de la phase numérique ------------------------------ */
    /* --------------------------------------------------------------------- */
    printf ("L:\n") ; cs_print (L, 0) ;

}