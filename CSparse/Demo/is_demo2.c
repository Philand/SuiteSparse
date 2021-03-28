#include "cs.h"
/* fichier de demo 2 IS pour la compréhension de CSparse */

csi is_min (csi *nonzeros, csi size, csi i)
{
    printf("\n ----- fonction is_min ----- \n") ;
    printf ("size = %td \n", size) ;
    csi min = 1000000000 ;
    printf ("min = nonzeros [0] = %td \n", min) ;
    for (csi k = 0 ; k < size ; k++)
    {
        printf ("k = %td  et  nonzeros [%td] = %td \n", k, k, nonzeros [k]) ;
        if (nonzeros[k] < min && nonzeros[k] != i)
        {
            min = nonzeros[k] ;
            printf ("new min = %td \n", nonzeros [k]) ;
        }

    }
    if (min == 1000000000)
        min = -1 ;
    printf("final min = %td \n", min) ;
    printf("--- fin fonction is_min --- \n \n") ;
    return min ;
}

/*

P est de taille n

union(csi* Pattern, const csi* L_k, int size )
{
    for(int i=0;i<size;++i)
    {
        const int index = *L_k[i];
        Pattern[index] = 1.0;
    }
}

on traite la colonne k :

Pattern = malloc(n, -1);

union(Pattern, A_k);

for all j tel que Parent[j] = k dans l'arbre d'élmination
    union(Pattern, L_j, size(L_j));


il faut faire croitre ses tableaux pour stocker les indices de la colonne k

int L_k_size = 0;
for int i=0;i<n;++i
    if(P[i] == 1) ++L_k_size;

non_zero_begin_k = la première cas autorisée dans non_zero pour stocker les indices de la colonne k
for int i=0;i<n;++i
    if(Pattern[i] == 1 ) non_zero[non_zero_begin_k++] = i;

// pas comme ça dans le code plus bas
old_size = L_size;
new_size = resize(L,size(L)+L_k_size);
ptr = old_size
for int i=0;i<n;++i
    if(Pattern[i] == 1 ) L_index[ptr++] = i;

*/

csi is_union (csi *list1, csi size1, csi *list2, csi size2, csi i)
{
    printf("\n ----- fonction is_union ----- \n") ;
    int already_in ;
    csi row, nb_new ;

    nb_new = 0 ;

    for (csi k = 0 ; k < size2 ; k++)
    {
        row = list2 [k] ;
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
    printf ("nb_new = %td \n", nb_new) ;
    printf("\n --- fin fonction is_union --- \n") ;
    return nb_new ;
}

int main (void)
{

    /* --- Récupération de la matrice --- */

    cs *T, *A ;
    T = cs_load (stdin) ;               /* load triplet matrix T from stdin */
    A = cs_compress (T) ;               /* A = compressed-column form of T */
    cs_spfree (T) ;                     /* clear T */

    /* --- Création de l'object pour notre factorisation --- */

    csn *N ;                             /* (cs_numeric) object initialization */
    N = cs_calloc (1, sizeof (csn)) ;    /* allocate result */

    /* --- Analyse Symbolique de la matrice A --- */
    csi k, i, n, *Ap, *Ai, *pi, *nonzeros, *indptr, count, nb_new, start , nb_row_col ;

    n = A->n ; Ap = A->p ; Ai = A->i ;
    indptr = cs_malloc (n+1, sizeof(csi)) ;
    indptr [0] = 0 ;
    pi = cs_malloc (n, sizeof(csi)) ;
    /*for (k = 0 ; k < n ; k++)
    {
        pi [k] = -1 ;
    }*/
    nonzeros = cs_malloc(n*(n+1)/2, sizeof(csi)) ;

    printf ("Tableau de pi à l'initialisation : \n") ;
    for (k = 0 ; k < n ; k++)
    {
        printf ("%td \n", pi [k]) ;
    }
    
    printf ("n = %td \n", n) ;

    for (k = 0 ; k < n ; k++)
    {
        printf ("----------------------------------------- k = %td \n", k) ;
        /* L_k = A_k */
        count = nb_row_col = Ap [k+1] - Ap [k] ;
        printf ("starting count = %td \n", count) ;
        start = Ap [k] ;
        printf ("start = Ap [%td] = %td \n", k, start) ;
        for (i = 0 ; i < nb_row_col ; i++)
        {
            if (Ai[start + i] >= k)
                nonzeros [indptr [k] + i] = Ai [start + i] ; 
            else
            {
                count -- ;
                continue ;
            }
        }

        /* for all i such that pi [k] = i */
        for (i = 0 ; i < k ; i++)
        {
            if (pi [i] == k)
            {
                printf ("!!! %td a pour parent %td !!! \n", i, k) ;
                count += is_union(&nonzeros [indptr [k]], count, &nonzeros[indptr [i]], indptr [i+1] - indptr [i], i) ;
            }
        }
        pi [k] = is_min (&nonzeros [indptr [k]], count, k) ;
        indptr [k+1] = indptr [k] + count ;
        // printf ("%td \n", indptr [k]) ;
        // printf ("%td \n", indptr [k+1]) ;
    }

    /* --- Vérification de la phase symbolique --- */
    printf (" \n --- Vérification de l'indice pointeur des colonnes : --- \n") ;
    for (k = 0 ; k < n + 1 ; k++)
    {
        printf ("%td \n", indptr [k]) ;
    }    


    printf (" \n --- Vérification de la structure des non-zéros : --- \n") ;
    for (k = 0 ; k < indptr [n] ; k++)
    {
        printf ("%td \n", nonzeros [k]) ;
    }

    printf ("\n --- Tableau de pi à la fin de l'exécution : --- \n") ;
    for (k = 0 ; k < n ; k++)
    {
        printf ("%td \n", pi [k]) ;
    }


    /* --- Factorisation Numérique de la matrice A --- */
    /*
    for (k = 0 ; k < n ; k++)
    {

    }

    */

}