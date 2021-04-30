#include "cs_demo.h"
#include <time.h>

/* 1 if A is square & upper tri., -1 if square & lower tri., 0 otherwise */
static csi is_sym (cs *A)
{
    csi is_upper, is_lower, j, p, n = A->n, m = A->m, *Ap = A->p, *Ai = A->i ;
    if (m != n) return (0) ;
    is_upper = 1 ;
    is_lower = 1 ;
    for (j = 0 ; j < n ; j++)
    {
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            if (Ai [p] > j) is_upper = 0 ;
            if (Ai [p] < j) is_lower = 0 ;
        }
    }
    return (is_upper ? 1 : (is_lower ? -1 : 0)) ;
}

/* true for off-diagonal entries */
static csi dropdiag (csi i, csi j, double aij, void *other) { return (i != j) ;}

/* C = A + triu(A,1)' */
static cs *make_sym (cs *A)
{
    cs *AT, *C ;
    AT = cs_transpose (A, 1) ;          /* AT = A' */
    cs_fkeep (AT, &dropdiag, NULL) ;    /* drop diagonal entries from AT */
    C = cs_add (A, AT, 1, 1) ;          /* C = A+AT */
    cs_spfree (AT) ;
    return (C) ;
}

/* create a right-hand side */
static void rhs (double *x, double *b, csi m)
{
    csi i ;
    for (i = 0 ; i < m ; i++) b [i] = 1 + ((double) i) / m ;
    for (i = 0 ; i < m ; i++) x [i] = b [i] ;
}

/* infinity-norm of x */
static double norm (double *x, csi n)
{
    csi i ;
    double normx = 0 ;
    for (i = 0 ; i < n ; i++) normx = CS_MAX (normx, fabs (x [i])) ;
    return (normx) ;
}

/* compute residual, norm(A*x-b,inf) / (norm(A,1)*norm(x,inf) + norm(b,inf)) */
static void print_resid (csi ok, cs *A, double *x, double *b, double *resid)
{
    csi i, m, n ;
    if (!ok) { printf ("    (failed)\n") ; return ; }
    m = A->m ; n = A->n ;
    for (i = 0 ; i < m ; i++) resid [i] = -b [i] ;  /* resid = -b */
    cs_gaxpy (A, x, resid) ;                        /* resid = resid + A*x  */
    printf ("resid: %8.2e\n", norm (resid,m) / ((n == 0) ? 1 :
        (cs_norm (A) * norm (x,n) + norm (b,m)))) ;
}

static double tic (void) { return (clock () / (double) CLOCKS_PER_SEC) ; }
static double toc (double t) { double s = tic () ; return (CS_MAX (0, s-t)) ; }

static void print_order (csi order)
{
    switch (order)
    {
        case 0: printf ("natural    ") ; break ;
        case 1: printf ("amd(A+A')  ") ; break ;
        case 2: printf ("amd(S'*S)  ") ; break ;
        case 3: printf ("amd(A'*A)  ") ; break ;
    }
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

iss *is_symalloc (cs *A)
{
    if (!A) return (NULL) ;
    csi n = A->n ;
    csi *A_colptr = A->p ;
    iss *S = cs_calloc (1, sizeof (iss)) ;
    S->L_colptr = cs_malloc (n+1, sizeof(csi)) ;
    S->L_rowptr = cs_malloc (n+1, sizeof(csi)) ;
    S->L_rowind = cs_malloc(A_colptr [n], sizeof(csi)) ;
    S->L_colind = cs_malloc(A_colptr [n], sizeof(csi)) ;
    S->parent = cs_malloc (n, sizeof(csi)) ;
    return S ;
}

/* free a symbolic left-looking Cholesky factorization */
iss *is_sfree (iss *S)
{
    if (!S) return (NULL) ;     /* do nothing if S already NULL */
    cs_free (S->parent) ;
    cs_free (S->L_colptr) ;
    cs_free (S->L_rowind) ;
    cs_free (S->L_rowptr) ;
    cs_free (S->L_colind) ;
    return ((iss *) cs_free (S)) ;  /* free the iss struct and return NULL */
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

/* read a problem from a file; use %g for integers to avoid csi conflicts */
problem *get_problem (FILE *f, double tol)
{
    cs *T, *A, *C ;
    csi sym, m, n, mn, nz1, nz2 ;
    problem *Prob ;
    Prob = cs_calloc (1, sizeof (problem)) ;
    if (!Prob) return (NULL) ;
    T = cs_load (f) ;                   /* load triplet matrix T from a file */
    Prob->A = A = cs_compress (T) ;     /* A = compressed-column form of T */
    cs_spfree (T) ;                     /* clear T */
    if (!cs_dupl (A)) return (free_problem (Prob)) ; /* sum up duplicates */
    Prob->sym = sym = is_sym (A) ;      /* determine if A is symmetric */
    m = A->m ; n = A->n ;
    mn = CS_MAX (m,n) ;
    nz1 = A->p [n] ;
    cs_dropzeros (A) ;                  /* drop zero entries */
    nz2 = A->p [n] ;
    if (tol > 0) cs_droptol (A, tol) ;  /* drop tiny entries (just to test) */
    Prob->C = C = sym ? make_sym (A) : A ;  /* C = A + triu(A,1)', or C=A */
    if (!C) return (free_problem (Prob)) ;
    printf ("\n--- Matrix: %g-by-%g, nnz: %g (sym: %g: nnz %g), norm: %8.2e\n",
            (double) m, (double) n, (double) (A->p [n]), (double) sym,
            (double) (sym ? C->p [n] : 0), cs_norm (C)) ;
    if (nz1 != nz2) printf ("zero entries dropped: %g\n", (double) (nz1 - nz2));
    if (nz2 != A->p [n]) printf ("tiny entries dropped: %g\n",
            (double) (nz2 - A->p [n])) ;
    Prob->b = cs_malloc (mn, sizeof (double)) ;
    Prob->x = cs_malloc (mn, sizeof (double)) ;
    Prob->resid = cs_malloc (mn, sizeof (double)) ;
    return ((!Prob->b || !Prob->x || !Prob->resid) ? free_problem (Prob) : Prob) ;
}

/* free a problem */
problem *free_problem (problem *Prob)
{
    if (!Prob) return (NULL) ;
    cs_spfree (Prob->A) ;
    if (Prob->sym) cs_spfree (Prob->C) ;
    cs_free (Prob->b) ;
    cs_free (Prob->x) ;
    cs_free (Prob->resid) ;
    return (cs_free (Prob)) ;
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

/* --- is_left_schol.c --- */

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

    L_colptr = S->L_colptr ;
    L_colptr [0] = 0 ;
    L_rowptr = S->L_rowptr ;
    L_rowptr [0] = 0 ;
    L_colind = S->L_colind ;
    L_rowind = S->L_rowind ;
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

    cs_free (stack) ;
    cs_free (P) ;
    cs_free (flag) ;

    return S ;
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

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

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

/* --- is_left_cholsol.c --- */

csi *cs_pinv_identity (csi const *p, csi n)
{
    csi k, *pinv ;
    if (!p) return (NULL) ;                     /* p = NULL denotes identity */
    pinv = cs_malloc (n, sizeof (csi)) ;        /* allocate result */
    if (!pinv) return (NULL) ;                  /* out of memory */
    for (k = 0 ; k < n ; k++) pinv [k] = k ;    /* invert the permutation */
    return (pinv) ;                             /* return result */
}

/* quicksort sur Ci */
/* vecteur dense pour stocker les valeurs de Cx */
/* puis on les stocke dans Cx selon l'ordre d'indices de Ci */

int csiComparator ( const void * first, const void * second)
{
    csi firstCsi = * (const csi *) first ;
    csi secondCsi = * (const csi *) second ;
    return (int) firstCsi - secondCsi ;
}

/* build of a column, returning the numbers of elements */
void is_build_Cx_column (double * Cx, double *a, csi * Ci, csi nb_nz_col)
{
    csi k ;
    for (k = 0 ; k < nb_nz_col ; k++)
    {
        Cx [k] = a [Ci [k]] ;
        a [Ci [k]] = 0.0 ;
    }
}

/* C = A(p,q) where p and q are permutations of 0..m-1 and 0..n-1. */
cs *is_permute_so (const cs *A, const csi *pinv, csi values)
{
    csi t, j, k, nz = 0, m, n, *Ap, *Ai, *Cp, *Ci ;
    double *Cx, *Ax, *a ;
    cs *C ;
    if (!CS_CSC (A)) return (NULL) ;    /* check inputs */
    m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    C = cs_spalloc (m, n, Ap [n], values && Ax != NULL, 0) ;  /* alloc result */
    if (!C) return (cs_done (C, NULL, NULL, 0)) ;   /* out of memory */
    a = cs_malloc (n, sizeof (double)) ;
    for (k = 0 ; k < n ; k++)
        a [k] = 0.0 ;
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    for (k = 0 ; k < n ; k++)
    {
        Cp [k] = nz ;
        j = k ;
        for (t = Ap [j] ; t < Ap [j+1] ; t++)
        {
            a [pinv [Ai [t]]] =  Ax [t] ;
            Ci [nz++] = pinv ? (pinv [Ai [t]]) : Ai [t] ;
        }
        qsort (&Ci [Cp [k]], nz - Cp [k], sizeof (csi), csiComparator) ;
        is_build_Cx_column (&Cx [Cp [k]], a, &Ci [Cp [k]], nz - Cp [k]) ; 
    }
    Cp [n] = nz ;                       /* finalize the last column of C */
    return (cs_done (C, NULL, NULL, 1)) ;
}

/* x=A\b where A is symmetric positive definite; b overwritten with solution */
csi is_left_cholsol (csi order, const cs *A, double *b)
{
    double *x ;
    iss *S ;
    csn *N ;
    cs *A_perm ;
    csi n, ok, *Perm, *pinv ;
    if (!CS_CSC (A) || !b) return (0) ;     /* check inputs */
    n = A->n ;

    Perm = cs_amd (order, A) ;     /* P = amd(A+A'), or natural */
    pinv = cs_pinv (Perm, n) ;     /* find inverse permutation */
    cs_free (Perm) ;
    if (order && !pinv) return (0) ;
    if (order == 0)
        A_perm = cs_permute (A, pinv, NULL, 1) ;
    else if (order == 1)
    {
        A_perm = is_permute_so (A, pinv, 1) ;
        printf (" MATRIX A : \n") ;
        cs_print (A, 0) ;
        printf ("pinv = [ ") ;
        for (csi k = 0 ; k < n ; k++)
            printf ("%td ", pinv [k]) ;
        printf ("]\n") ;
        printf (" MATRIX A permuted : \n") ;
        cs_print (A_perm, 0) ;
    }

    S = is_left_schol (order, A_perm) ;     /* ordering and symbolic analysis */
    N = is_left_chol (A_perm, S) ;          /* numeric Cholesky factorization */
    printf ("L:\n") ; cs_print (N->L, 0) ;
    x = cs_malloc (n, sizeof (double)) ;    /* get workspace */
    ok = (S && N && x) ;
    if (ok)
    {
        cs_ipvec (pinv, b, x, n) ;   /* x = P*b */
        cs_lsolve (N->L, x) ;        /* x = L\x */
        cs_ltsolve (N->L, x) ;       /* x = L'\x */
        cs_pvec (pinv, x, b, n) ;    /* b = P'*x */
    }
    cs_free (x) ;
    cs_free (pinv) ;
    is_sfree (S) ;
    cs_nfree (N) ;
    cs_spfree (A_perm) ;
    return (ok) ;
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

csi test1 (problem *Prob)
{
    cs *A ;
    csi order, m, ok ;
    double *b, *x, *resid, t ;
    if (!Prob) return (0) ;
    A = Prob->A ; b = Prob->b ; x = Prob->x ; resid = Prob->resid ;
    m = A->m ;

    // if (!Prob->sym) return (1) ;
    for (order = 0 ; order <= 1 ; order++)      /* natural and amd(A+A') */
    {
        if (!order && m > 1000) continue ;
        printf ("Cholesky left-looking ") ;
        print_order (order) ;
        rhs (x, b, m) ;                         /* compute right-hand side */
        t = tic () ;
        ok = is_left_cholsol (order, A, x) ;    /* solve Ax=b with Cholesky */
        printf ("time: %8.2f ", toc (t)) ;
        print_resid (ok, A, x, b, resid) ;      /* print residual */
    }
    return (1) ;
}

/*
	- ./cs_demo1 < ../Matrix/t1
	- ./cs_demo2 < ../Matrix/t1
	- ./cs_demo2 < ../Matrix/ash219
	- ./cs_demo2 < ../Matrix/bcsstk01
	- ./cs_demo2 < ../Matrix/fs_183_1
	- ./cs_demo2 < ../Matrix/mbeacxc
	- ./cs_demo2 < ../Matrix/west0067
	- ./cs_demo2 < ../Matrix/lp_afiro
	- ./cs_demo2 < ../Matrix/bcsstk16
	- ./cs_demo3 < ../Matrix/bcsstk01
	- ./cs_demo3 < ../Matrix/bcsstk16
	- ./is_left_demo < ../Matrix/is_matrix1
	- ./is_left_demo < ../Matrix/is_matrix2
	- ./is_left_demo < ../Matrix/is_matrix0
*/

int main (void)
{
    problem *Prob = get_problem (stdin, 1e-14) ;
    test1 (Prob) ;
    free_problem (Prob) ;
    return (0) ;
}