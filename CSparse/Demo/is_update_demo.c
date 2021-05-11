#include "cs_demo.h"
#include <time.h>

/* -------------------------------------------------------------------------- */
/* ------------------------------- FUNCTIONS USEFULL ------------------------ */
/* -------------------------------------------------------------------------- */

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
    printf ("   %8.2e\n", norm (resid,m) / ((n == 0) ? 1 :
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
/* ------------------------------- UPDATE DEMO ------------------------------ */
/* -------------------------------------------------------------------------- */

csi *is_pre_update (csi *A_col_modified, csi n , const iss *S)
{
    // I_0 = A_col_modified | I_1 = L_col_modified
    csi k, i ;
    csi *parent, *L_col_modified ;
    parent = S->parent ;
    L_col_modified = cs_malloc (n, sizeof(csi)) ;
    for (k = 0 ; k < n ; k++)
        L_col_modified [k] = 0 ;

    // faire directement une boucle sur le tableau.

    for (k = 0 ; k < n ; k++)
    {
        if (A_col_modified [k] == 1)
        {
            L_col_modified [k] = 1 ;
            for (i = parent [k] ; i != -1 ; i = parent [i])
                L_col_modified [i] = 1 ;
        }
    }
    return L_col_modified ;
}

csn *is_left_cholupdate (const cs *A, const iss *S, csn *N, csi *L_col_modified)
{
    cs *L ;
    csi k, n, i, j ;
    csi *A_colptr, *A_rowind, *Lp, *Li, *L_rowptr, *L_colind, *Lpk ;
    double lkj, lkk ;
    double *Ax, *a, *Lx ;

    n = A->n ;
    L = N->L ;
    Lp = L->p ; Li = L->i ; Lx = L->x ;

    L_colind = S->L_colind ;
    L_rowptr =  S->L_rowptr ;
    a = cs_malloc(n, sizeof (double)) ;                 /* get csi workspace */
    Ax = A->x ;
    A_colptr = A->p ;
    A_rowind = A->i ;
    Lpk = cs_malloc(n, sizeof (csi)) ;             

    for (k = 0 ; k < n ; k++)
    {
        if (Lp [k+1] - Lp [k] > 1)
            Lpk [k] = Lp [k] + 1 ;
        else
            Lpk [k] = -1 ;
    }

    for (k = 0 ; k < n ; k++)
    {
        if (L_col_modified [k] == 1)
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
                    a [Li[p]] -= Lx [p]*lkj;
              
                Lpk [j] ++ ;
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
    }

    cs_free (Lpk) ;
    cs_free (a) ;

    return N ;
}

/* solve a linear system using Cholesky and after update factorization */
csi update_demo (problem *Prob)
{
    cs *A, *C ;
    double *b, *x, *resid,  t, tol ;
    csi k, m, n, ok, order, nb, ns, *r, *s, *rr, sprank, sym ;
    csd *D ;
    if (!Prob) return (0) ;
    A = Prob->A ; C = Prob->C ; b = Prob->b ; x = Prob->x ; resid = Prob->resid;
    m = A->m ; n = A->n ; sym = Prob->sym ;
    D = cs_dmperm (C, 1) ;                      /* randomized dmperm analysis */
    if (!D) return (0) ;
    nb = D->nb ; r = D->r ; s = D->s ; rr = D->rr ;
    sprank = rr [3] ;
    for (ns = 0, k = 0 ; k < nb ; k++)
    {
        ns += ((r [k+1] == r [k]+1) && (s [k+1] == s [k]+1)) ;
    }
    printf ("blocks: %g singletons: %g structural rank: %g\n",
        (double) nb, (double) ns, (double) sprank) ;
    cs_dfree (D) ;

    /* Cholesky left-looking */
    if (sym == 1 || sym == -1)
        A = make_sym (A) ;

    for (order = 0 ; order <= 0 ; order++)      /* natural and amd(A+A') */
    {
        /* --------------------------------------------------------- */
        /* update of A --> A' (test for is_matrix_1 / is_matrix1_up) */
        /* --------------------------------------------------------- */

        if (!order && m > 1000) continue ;
        printf ("cholesky left-looking update ") ;

        double *x ;
        iss *S ;
        csn *N ;
        cs *C ;
        csi *Perm, *pinv ;
        if (!CS_CSC (A) || !b) return (0) ;     /* check inputs */
        n = A->n ;

        Perm = cs_amd (order, A) ;     /* P = amd(A+A'), or natural */
        pinv = cs_pinv (Perm, n) ;     /* find inverse permutation */
        cs_free (Perm) ;
        if (order && !pinv) return (0) ;

        C = is_symperm (A, pinv, 1) ;

        S = is_left_schol (order, C) ;     /* ordering and symbolic analysis */
        N = is_left_chol (C, S) ;          /* numeric Cholesky factorization */
        x = cs_malloc (n, sizeof (double)) ;    /* get workspace */
        printf ("L :\n") ; cs_print (N->L, 0) ;

        A->x [13] = 8 ; A->x [14] = 2 ; A->x [20] = 8 ; A->x [24] = 2 ;
        csi *A_col_modified, *L_col_modified ;
        A_col_modified = cs_malloc (n, sizeof (csi)) ;
        for (k = 0 ; k < n ; k++)
            A_col_modified [k] = 0 ;
        A_col_modified [4] = 1 ;

        /*printf ("A_col_modified = [") ;
        for (k = 0 ; k < n ; k++)
            printf ("%td ", A_col_modified [k]) ;
        printf ("]\n") ;*/

        L_col_modified = is_pre_update (A_col_modified, n , S) ;
        /*printf ("L_col_modified = [") ;
        for (k = 0 ; k < n ; k++)
            printf ("%td ", L_col_modified [k]) ;
        printf ("]\n") ;*/

        N = is_left_cholupdate (A, S, N, L_col_modified) ;

        printf ("L_updated :\n") ; cs_print (N->L, 0) ;

        cs_free (x) ;
        cs_free (pinv) ;
        is_sfree (S) ;
        cs_nfree (N) ;
        cs_spfree (C) ;
        cs_free(A_col_modified) ;
        cs_free(L_col_modified) ;

    }
    printf ("------------------------------------------------------------------------------------- \n") ;
    return (1) ;
} 

/* -------------------------------------------------------------------------- */
/* ---------------------------------- MAIN ---------------------------------- */
/* -------------------------------------------------------------------------- */

int main (void)
{
    problem *Prob = get_problem (stdin, 0) ;
    update_demo (Prob) ;
    free_problem (Prob) ;
    printf ("\n") ;
    return (0) ;
}