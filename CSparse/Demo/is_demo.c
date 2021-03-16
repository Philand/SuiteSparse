#include "cs.h"
/* fichier de demo IS pour la compréhension de CSparse */

int main (void)
{
    cs *T, *A ;
    csi k, order, mn, m, n ;
    double *b ;

    /* ---------------------------------------------------------------------- */
    
    T = cs_load (stdin) ;               /* load triplet matrix T from stdin */
    printf ("T:\n") ; cs_print (T, 0) ; /* print T */
    A = cs_compress (T) ;               /* A = compressed-column form of T */
    printf ("A:\n") ; cs_print (A, 0) ; /* print A */
    cs_spfree (T) ;                     /* clear T */

    /* ---------------------------------------------------------------------- */

    printf ("Execution of the cs_cholsol function on the matrix A : \n") ;

    m = A->m ;
    n = A-> n ;
    mn = CS_MAX (m,n) ;
    b = cs_malloc (mn, sizeof (double)) ;
    for (k = 0 ; k < n ; k++)
        b [k] = 1 + ((double) k) / m ;

    order = 1 ;
    cs_cholsol (order, A, b) ; /* solve Ax=b with Cholesky */

    printf ("Successful execution !")

    /* ---------------------------------------------------------------------- */



    return (0) ;
}
