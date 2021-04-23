#include "cs.h"

/* union of 2 columns of their nonzeros */
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

/* build of a column, returning the numbers of elements */
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

/* build of a line (left columns) */
void is_write (csi *L_colind, csi *stack, csi top, csi n)
{
    csi imax = n - top ;
    for (csi i = 0 ; i < imax  ; i++)
        L_colind [i] = stack [i] ;
}

/* symbolic analysis for a left-looking Cholesky factorization */
iss *is_left_schol (csi order, const cs *A)
{
    if (!CS_CSC (A)) return (NULL) ;    /* check inputs */

    iss *S ;
    // csi *Perm ;
    csi k, i, j, p; /* variables incrémentales */
    csi len, top ;
    csi nb_nz_col ;
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

    S = is_symalloc (A) ;   /* allocate result S */
    if (!S) return (NULL) ; /* out of memory */

    n = A->n ; A_colptr = A->p ; A_rowind = A->i ;
    // L_colptr = cs_malloc (n+1, sizeof(csi)) ;
    L_colptr = S->L_colptr ;
    L_colptr [0] = 0 ;
    // L_rowptr = cs_malloc (n+1, sizeof(csi)) ;
    L_rowptr = S->L_rowptr ;
    L_rowptr [0] = 0 ;
    // L_rowind = cs_malloc(A_colptr [n], sizeof(csi)) ;
    L_rowind = S-> L_rowind ;
    //L_colind = cs_malloc(A_colptr [n], sizeof(csi)) ;
    L_colind = S-> L_colind ;
    //parent = cs_malloc (n, sizeof(csi)) ;
    parent = S->parent ;
    flag = cs_malloc (n, sizeof(csi)) ;
    stack = cs_malloc (n, sizeof(csi)) ;
    P = cs_malloc (n, sizeof(csi)) ;

    for (k = 0 ; k < n ; k++)
        P [k] = 0 ;

    for (k = 0 ; k < n ; k++)
    {
        parent [k] = -1 ;
        flag [k] = k ;
        top = n ;

        /* --- structure of line k ------------------------------------------ */
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

        /* --- structure of column k ---------------------------------------- */
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
        L_colptr [k+1] = L_colptr [k] + nb_nz_col ;
    }

    cs_free (stack) ;
    cs_free (P) ;
    cs_free (flag) ;

    // S->L_colind = L_colind ;
    // S->L_colptr = L_colptr ;
    // S->L_rowind = L_rowind ;
    // S->L_rowptr = L_rowptr ;
    // S->parent = parent ;

    // cs_free (L_rowptr) ;
    // cs_free (L_colptr) ;
    // cs_free (L_rowind) ;
    // cs_free (L_colind) ;
    // cs_free (parent) ;

    return S ;
}