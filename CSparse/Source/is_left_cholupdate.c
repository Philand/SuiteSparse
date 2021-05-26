#include "cs.h"

/* --- is_left_cholupdate */

csi is_ind_in_set (csi i, csi *set, csi set_size)
{
    csi k = 0 ;
    while (set [k] != i && k < set_size)
        k++ ;

    if (k == set_size)
        return 0 ;
    else 
        return 1 ;
}

csi is_maj_Lpk (csi Lp_j, csi Lp_j1, csi *Li, csi *I1, csi I1_size)
{
    csi value, value2 ;
    csi result ; /* Lpk [j] */
    value = Lp_j1 ;
    if (value - Lp_j > 1)
    {
        result = Lp_j + 1 ;
        value2 = value - Lp_j ;
        while (!is_ind_in_set(Li [result], I1, I1_size) && (value2 > 1))
        {
            result++ ;
            value2-- ;
        }
        if (value2 == 1)
            result = -1 ;
    }
    else
        result = -1 ;
        
    return result ;
}

csi is_unique (csi *vec, csi i, csi *size)
{
    csi k ;
    for (k = 0 ; k < (*size) ; k++)
    {
        if (vec [k] == i)
            return (0) ;
    }
    return (1) ;
}

csi is_add_entry (csi *vec, csi i, csi *size, csi *max_size)
{
    csi ok = 1 ;
    if (!vec || i < 0) return (0) ; /* check inputs */
    if (!is_unique (vec, i, size)) return (1) ; /* check if i is already in I0 */
    if (*size >= (*max_size))
    {
        (*max_size) *= 2 ;
        vec = cs_realloc (vec, (*max_size), sizeof (csi), &ok) ;
    }
    vec [(*size)++] = i ;
    return (1) ;
}

/* load a set I0 and changed values of cs matrix A from a file */
csi *is_load_update_matrix (FILE *f, cs *A, csi *I0, csi *I0_size)
{
    csi i, j, max_size, ok = 1 ;
    double x ;
    if (!f) return (NULL) ; /* check inputs */
    max_size = 2 ;
    while (fscanf (f, "%td %td %lf\n", &i, &j, &x) == 3)
    {
        A->x [(csi) i] = x ;
        if (!I0 || j < 0) return cs_free (I0) ; /* check inputs */
        if (!is_unique (I0, j, I0_size)) continue ; /* check if i is already in I0 */
        if (*I0_size >= max_size)
        {
            max_size *= 2 ;
            I0 = cs_realloc (I0, max_size, sizeof (csi), &ok) ;
        }
        I0 [(*I0_size)++] = j ;
    }
    qsort (I0, *I0_size, sizeof (csi), csiComparator) ;
    return (I0) ;
}

/* fonction construisant l'ensemble I1 à partie de l'ensemble I0 */
csi *is_pre_update (csi *I0, csi I0_size, csi *I1, csi *I1_size, const iss *S)
{
    csi k, i, j, count, I1_max_size, ok = 1 ;
    csi *parent ;
    parent = S->parent ;
    I1_max_size = I0_size ;
    ok = 1 ;

    i = count = 0 ;
    for (k = I0 [i] ; i < I0_size ; k = I0 [++i] )
    {
        for (j = k ; j != -1 ; j = parent [j])
        {
            if (!is_unique (I1, j, &count)) continue ;
            I1 [count++] = j ;
            if (count == I1_max_size)
            {
                I1_max_size *= 2 ;
                I1 = cs_realloc (I1, I1_max_size, sizeof (csi), &ok) ;
            }
        }
    }
    *I1_size = count ;
    qsort (I1, *I1_size, sizeof (csi), csiComparator) ;
    return (I1) ;
}

/* fonction effectuant la mise à jour partielle de Cholesky left-looking */
csn *is_left_cholupdate (const cs *A, const iss *S, csn *N, csi *I1, csi I1_size)
{
    cs *L ;
    csi k, n, i, j, o ;
    csi *A_colptr, *A_rowind, *Lp, *Li, *L_rowptr, *L_colind, *Lpk ;
    double lkj, lkk ;
    double *Ax, *a, *Lx ;

    n = A->n ;
    L = N->L ;
    Lp = L->p ; Li = L->i ; Lx = L->x ;
    o = 0 ;

    L_colind = S->L_colind ;
    L_rowptr =  S->L_rowptr ;
    a = cs_malloc (n, sizeof (double)) ; /* get csi workspace */
    Ax = A->x ;
    A_colptr = A->p ;
    A_rowind = A->i ;
    Lpk = cs_malloc (n, sizeof (csi)) ;     

    for (k = 0 ; k < n ; k++)
        Lpk [k] = is_maj_Lpk (Lp [k], Lp [k+1], Li, I1, I1_size) ;

    for (k = I1 [o] ; o < I1_size ; k = I1 [++o] )
    {
        /* a (k:n) = A (k:n,k) */
        for (i = A_colptr [k] ; i < A_colptr [k+1] ; i++)
            a [A_rowind [i]] = Ax [i] ;

        /* for j = find (L (k,;)) */
        for (i = L_rowptr [k] ; i < L_rowptr [k+1] ; i++)
        {
            j = L_colind [i] ;
            lkj = Lx [Lpk [j]] ;
              
            for (int p = Lpk [j] ; p < Lp [j+1] ; p++)
                a [Li[p]] -= Lx [p]*lkj;
              
            // mise à jour de Lpk [j]
            Lpk [j] = is_maj_Lpk (Lpk [j], Lp [j+1], Li, I1, I1_size) ;
        }

        /* L (k,k) = sqrt (a (k)) */ 
        Lx [ Lp [k]] = lkk = sqrt (a [Li [Lp [k]]]) ;
        a [Li [Lp [k]]] = 0.0 ;

        /* L (k+1:n,k) = a (k+1:n) / L (k,k) */
        for (j = Lp [k] + 1 ; j < Lp [k+1] ; j++)
        {
            Lx [j] = a [Li [j]] / lkk ;
            a [Li [j]] = 0.0 ;
        }
    }

    cs_free (Lpk) ;
    cs_free (a) ;

    return N ;
}