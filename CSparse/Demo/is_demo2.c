#include "cs.h"
/* fichier de demo 2 (InSimo) pour la compréhension de CSparse */
/* première implémentation de la left-looking                  */

/* fonction is_min : calcule le minimum d'une liste (nonzeros) d'une taille
   fixé (size) en excluant l'indice de la colonne (i)                       */
csi is_min (csi *nonzeros, csi size, csi i)
{
    csi min = 1000000000 ;
    for (csi k = 0 ; k < size ; k++)
    {
        if (nonzeros[k] < min && nonzeros[k] != i)
        {
            min = nonzeros[k] ;
        }

    }
    if (min == 1000000000)
        min = -1 ;

    return min ;
}

/* fonction is_union : calcule l'union de deux listes */
csi is_union (csi *list1, csi size1, csi *list2, csi size2, csi i)
{
    int already_in ;
    csi row, nb_new = 0 ;

    for (csi k = 0 ; k < size2 ; k++)
    {
        row = list2 [k] ;
        already_in = 0 ;
        for (csi j = 0 ; j < size1 ; j++)
        {
            if ( list1[j] == row)
            {
                already_in = 1 ;
            }
        }
        if (already_in == 0 && row != i)
        {
            list1 [size1 + nb_new] = row ;
            nb_new ++ ;
        }
    }

    return nb_new ;
}

int main (void)
{
    printf ("=================================================================> IS_DEMO 2 \n") ;

    /* --- Récupération de la matrice --------------------------------------- */
    cs *T, *A ;
    T = cs_load (stdin) ;               /* load triplet matrix T from stdin */
    A = cs_compress (T) ;               /* A = compressed-column form of T */
    cs_spfree (T) ;                     /* clear T */

    /* --- Analyse Symbolique de la matrice A ------------------------------- */
    printf (" \n ------- Analyse Symbolique de la matrice A ------- \n") ;
    csi k, i, n, *Ap, *Ai, *pi, *nonzeros, *indptr, count, nb_new, start ;
    csi nb_row_col, nb_sup_values, *start_update ; 

    n = A->n ; Ap = A->p ; Ai = A->i ;
    indptr = cs_malloc (n+1, sizeof(csi)) ;
    indptr [0] = 0 ;
    pi = cs_malloc (n, sizeof(csi)) ;
    nonzeros = cs_malloc(n*(n+1)/2, sizeof(csi)) ;
    start_update = cs_malloc(n, sizeof (csi)) ;

    for (k = 0 ; k < n ; k++)
    {
        /* L_k = A_k */
        count = nb_row_col = Ap [k+1] - Ap [k] ;
        start = Ap [k] ;
        nb_sup_values = 0 ;

        for (i = 0 ; i < nb_row_col ; i++)
        {
            /* on filtre les valeurs de la triangulaire inférieure de A */
            if (Ai[start + i] >= k)
            {
                nonzeros [indptr [k] + i - nb_sup_values] = Ai [start + i] ; 
            }
            else
            {
                nb_sup_values += 1 ; 
                count -- ;
            }
        }

        /* for all i such that pi [k] = i */
        for (i = 0 ; i < k ; i++)
        {
            if (pi [i] == k)
            {
                count += is_union(&nonzeros [indptr [k]], count, &nonzeros[indptr [i]], indptr [i+1] - indptr [i], i) ;
            }
        }
        pi [k] = is_min (&nonzeros [indptr [k]], count, k) ;    /* upd etree */
        indptr [k+1] = indptr [k] + count ;                     /* upd nb_nz */
        if (indptr [k+1] - indptr [k] > 1)
            start_update [k] = nonzeros [indptr [k] + 1] ;      /* upd st_cl */
        else
            start_update [k] = -1 ;
    }

    /* --- Vérification de la phase symbolique ------------------------------ */
    printf (" \n --- Vérification de l'indice pointeur des colonnes : --- \n") ;
    for (k = 0 ; k < n + 1 ; k++)
    {
        printf ("%td  ;  ", indptr [k]) ;
    }
    printf ("\n") ;

    printf (" \n --- Vérification de la structure des non-zéros : --- \n") ;
    for (k = 0 ; k < indptr [n] ; k++)
    {
        printf ("%td  ;  ", nonzeros [k]) ;
    }
    printf ("\n") ;

    printf ("\n --- Tableau de pi à la fin de l'exécution : --- \n") ;
    for (k = 0 ; k < n ; k++)
    {
        printf ("%td  ;  ", pi [k]) ;
    }
    printf ("\n") ;

    /* --- Factorisation Numérique de la matrice A -------------------------- */
    printf (" \n ------- Factorisation Numérique de la matrice A ------- \n") ;
    csn *N ;        /* (cs_numeric) object initialization */
    cs *L ;
    csi nb_nonzeros, id, *Lp, *Li, j, p, q, *count_update ;
    double *a, *Ax, *Lx, lkj, lkk ;

    N = cs_calloc (1, sizeof (csn)) ;                   /* allocate result */
    a = cs_malloc(n, sizeof (double)) ;                 /* get csi workspace */
    Ax = A->x ; Ap = A->p ; Ai = A->i ;
    N->L = L = cs_spalloc (n, n, indptr [n], 1, 0) ;    /* allocate result */
    Lp = L->p ; Li = L->i ; Lx = L->x ;
    count_update = cs_malloc(n, sizeof (csi)) ;

    /* on connait le nombre de valeurs dans chaque colonne */
    for (k = 0 ; k < n+1 ; k++) 
        Lp [k] = indptr [k] ;

    /* on initialise les compteurs d'accès aux colonnes */
    for (k = 0 ; k < n ; k++)
    {
        if (Lp [k+1] - Lp [k] > 1)
            count_update [k] = 1 ;
        else
            count_update [k] = 0 ;
    }
    
    for (k = 0 ; k < n ; k++)
    {
        /* a (k:n) = A (k:n,k) */
        for (i = Ap [k] ; i < Ap [k+1] ; i++)
        {
            a [Ai [i]] = Ax [i] ;
        }

        /* for j = find (L (k,;)) */
        for (j = 0 ; j < k ; j++)
        {
            /* on traite colonne par colonne */
            if (start_update [j] == k)
            {
                p = count_update [j] ; 
                q = 0 ;
                lkj = Lx [ Lp [j] + p ] ;

                while ((q != (Lp [k+1] - Lp [k] + 1)) && (p != (Lp [j+1] - Lp [j])))
                {
                    if (Li [Lp[j] + p] == nonzeros [Lp [k] + q])
                    {
                        a [ nonzeros [Lp [k] + q]] -= Lx [Lp [j] + p] * lkj ;
                        p ++ ; 
                        q ++ ;
                    }
                    else if (Li [Lp[j] + p] < nonzeros [Lp [k] + q])
                    {
                        p ++ ;
                    }
                    else if (Li [Lp[j] + p] > nonzeros [Lp [k] + q])
                    {
                        q ++ ;
                    }
                }
                count_update [j] ++ ;
                start_update [j] = nonzeros [indptr [j] + count_update [j]] ;
            }
        }

        /* L (k,k) = sqrt (a (k)) */
        Lx [ Lp [k]] = lkk = sqrt (a [nonzeros [Lp [k]]]) ;
        Li [ Lp [k]] = nonzeros [Lp [k]];

        /* L (k+1:n,k) = a (k+1:n) / L (k,k) */
        for (j = Lp [k] + 1 ; j < Lp [k+1] ; j++)
        {
            Lx [j] = a [nonzeros [j]] / lkk ;
            Li [j] = nonzeros [j];
        }
    }

    /* --- Vérification de la phase numérique ------------------------------ */
    printf ("L:\n") ; cs_print (L, 0) ;

}