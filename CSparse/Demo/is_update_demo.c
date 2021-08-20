#include "cs_demo.h"
#include "cs.h"
#include <time.h>
#include <string.h>

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
    //cs_dropzeros (A) ;                  /* drop zero entries */
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
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

csi is_ind_in_set (csi i, csi *set, csi set_size)
{
    csi k = 0, pos = -1 ;
    csi inf = 0, sup = set_size - 1, mil ;
    /* valeur à rechercher : i */
    while ((inf<=sup) && (pos==-1))
    {
        mil = (sup+inf)/2 + (sup+inf)%2 ;
        if (i < set [mil])
            sup = mil - 1 ;
        else if (i > set [mil])
            inf = mil + 1 ;
        else
            pos = mil ;
    }

    /* terminer */

    if (pos == -1)
        return 0 ;
    else 
        return 1 ;
}

/*
for (k = 0 ; k < n ; k++)
    Lpk [k] = is_init_Lpk (Lp [k], Lp [k+1], Li, I1, I1_size) ;
*/

csi is_init_Lpk (csi Lp_j, csi Lp_j1, csi *Li, csi *I1, csi I1_size)
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

/* load a set I0 and changed values of cs matrix A from a file */
csi *is_load_update_matrix2 (FILE *f, cs *A, csi *I0, csi *I0_size)
{
    csi i, j, k, max_size, ok = 1, count = 0;
    double x ;
    csi nb_col = A->n;
    if (!f) return (NULL) ; /* check inputs */
    csi *I0_bool;
    I0_bool = cs_malloc(nb_col, sizeof(csi)) ;
    memset(I0_bool, 0, nb_col*sizeof(csi)) ;

    while (fscanf (f, "%td %td %lf\n", &i, &j, &x) == 3)
    {
        A->x [(csi) i] = x ;
        if (I0_bool[j] == 0)
        {
            I0_bool[j] = 1;
            count++;
        }
    }
    *I0_size = count;
    I0 = cs_malloc (2*(*I0_size), sizeof (csi)) ;
    count = 0;
    for (k = 0 ; k < nb_col ; k++)
    {
        if (I0_bool [k] == 1)
            I0 [count++] = k;
    }   
    cs_free(I0_bool) ;
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

csi *is_pre_update2 (csi *I0, csi I0_size, csi *I1, csi*I1_size, const iss *S, csi nb_col)
{
    csi k, i, j, count, I1_max_size, ok = 1 ;
    csi *parent ;
    parent = S->parent ;
    //I1_max_size = I0_size ;
    ok = 1 ;

    count = 0 ;

    csi *I1_bool ;
    I1_bool = cs_malloc(nb_col, sizeof(csi)) ;
    memset(I1_bool, 0, nb_col*sizeof(csi)) ;

    for (k = 0 ; k < I0_size ; k++)
    {
        i = I0[k] ;
        for (j = i ; j != -1 ; j = parent[j])
        {
            if (I1_bool[j] == 1) // no need to climb up further, already done.
                break ;
            I1_bool[j] = 1 ;
            count++ ;
        }
    }

    I1_max_size = count;
    *I1_size = count;
    I1 = cs_malloc (2*I1_max_size, sizeof (csi)) ;
    count = 0;
    for (k = 0 ; k < nb_col ; k++)
    {
        if (I1_bool [k] == 1)
            I1 [count++] = k;
    }   
    cs_free(I1_bool) ;
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
        Lpk [k] = is_init_Lpk (Lp [k], Lp [k+1], Li, I1, I1_size) ;

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
            Lpk [j]++ ;
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

/* x=A\b where A is symmetric positive definite; b overwritten with solution */
csi is_left_cholsol_update (csi order, cs *A, double *b, FILE * filePtr)
{
    double *x, t ;
    iss *S ;
    csn *N ;
    cs *C ;
    csi n, ok, *Perm, *pinv ;
    csi *I0, *I1 ;
    csi I0_size, I1_size ;
    if (!CS_CSC (A) || !b) return (0) ;     /* check inputs */
    n = A->n ;

    Perm = cs_amd (order, A) ;     /* P = amd(A+A'), or natural */
    pinv = cs_pinv (Perm, n) ;     /* find inverse permutation */
    cs_free (Perm) ;
    if (order && !pinv) return (0) ;

    C = is_symperm (A, pinv, 1) ;      /* permuting matrix */
    S = is_left_schol (order, C) ;     /* ordering and symbolic analysis */
    N = is_left_chol (C, S) ;          /* numeric Cholesky factorization */
    x = cs_malloc (n, sizeof (double)) ;    /* get workspace */
    print_fill_in (C, N) ;
    //printf ("L = \n") ; cs_print (N->L, 0) ;

    I0_size = I1_size = 0 ;
    I0 = cs_malloc (2, sizeof (csi)) ;
    I0 = is_load_update_matrix2 (filePtr, A, I0, &I0_size) ;

    C = is_symperm (A, pinv, 1) ;
    t = tic () ;
    // I1 = cs_malloc (I0_size, sizeof (csi)) ;
    I1 = is_pre_update2 (I0, I0_size, I1, &I1_size, S, n) ;

    N = is_left_cholupdate (C, S, N, I1, I1_size) ;

    cs_free (I0) ;
    cs_free (I1) ;

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
    cs_spfree (C) ;
    printf (" %8.6f s | ", toc (t)) ;
    return (ok) ;
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/* ------------------------------- UPDATE DEMO ------------------------------ */
/* -------------------------------------------------------------------------- */

/* solve a linear system using Cholesky and after update factorization */
csi update_demo (problem *Prob, FILE * filePtr)
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

    if (sym == 1 || sym == -1)
        A = make_sym (A) ;

    printf ("--------------------------------------------------------------------------------------------- \n") ;
    printf ("          algorithm           |    permutation   |    fill-in    |    time    |    residue    \n") ;
    printf ("------------------------------|------------------|---------------|------------|-------------- \n") ;

    /* ------------------------------------------------------------ */

    /* Updated Cholesky left-looking factorisation */
    for (order = 0 ; order <= 1 ; order++)      /* natural and amd(A+A') */
    {
        if (!order && m > 1000) continue ;
        printf ("partial cholesky left-looking |     ") ;
        print_order (order) ; printf ("  |") ;
        rhs (x, b, m) ;                         /* compute right-hand side */
        //printf ("A = \n") ; cs_print (A, 0) ;
        ok = is_left_cholsol_update (order, A, x, filePtr) ;    /* solve Ax=b with Cholesky */
        print_resid (ok, A, x, b, resid) ;      /* print residual */
        // printf("\n");
    }

    /* Entire Cholesky left-looking factorization */
    for (order = 0 ; order <= 1 ; order++)      /* natural and amd(A+A') */
    {
        if (!order && m > 1000) continue ;
        printf ("entire cholesky left-looking  |     ") ;
        print_order (order) ; printf ("  |") ;
        rhs (x, b, m) ;                         /* compute right-hand side */
        t = tic () ;
        //printf ("A = \n") ; cs_print (A, 0) ;
        ok = is_left_cholsol (order, A, x) ;    /* solve Ax=b with Cholesky */
        printf (" %8.6f s | ", toc (t)) ;
        print_resid (ok, A, x, b, resid) ;      /* print residual */
        // printf("\n");
    }
    printf ("--------------------------------------------------------------------------------------------- \n") ;
    
    return (1) ;
} 

/* -------------------------------------------------------------------------- */
/* ---------------------------------- MAIN ---------------------------------- */
/* -------------------------------------------------------------------------- */

int main (int argc, char *argv[])
{
    if (argc < 2 || argc > 3)
    {
        fprintf ( stderr , " usage : ./is_update_demo vec_upt < matrix\n") ;
        exit (1) ;
    }
    else
    {
        char * fileName = NULL ;
        FILE * filePtr = NULL ;

        fileName = argv [1] ;
        filePtr = fopen (fileName, "r") ;

        problem *Prob = get_problem (stdin, 0) ;
        update_demo (Prob, filePtr) ;

        fclose (filePtr) ;
        free_problem (Prob) ;
        printf ("\n") ;
        return (0) ;
    }
}