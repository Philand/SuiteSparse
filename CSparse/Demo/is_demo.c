#include "cs.h"
/* fichier de demo IS pour la compréhension de CSparse */

int main (void)
{
    cs *T, *A, *Eye, *AT, *C, *D ;
    csi i, m ;
    
    T = cs_load (stdin) ;               /* load triplet matrix T from stdin */
    printf ("T:\n") ; cs_print (T, 0) ; /* print T */
    A = cs_compress (T) ;               /* A = compressed-column form of T */
    printf ("A:\n") ; cs_print (A, 0) ; /* print A */
    cs_spfree (T) ;                     /* clear T */

    return (0) ;
}
