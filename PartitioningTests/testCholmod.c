#include <cholmod_blas.h>
#include <cholmod.h>

/* ordering methods */
#define NESTED_DISSECTION 0
#define METIS 1
#define BISECT 2

/* ff is a global variable so that it can be closed by my_handler */
FILE *ff ;

/* halt if an error occurs */
static void my_handler (int status, const char *file, int line,
    const char *message)
{
    printf ("cholmod error: file: %s line: %d status: %d: %s\n",
	    file, line, status, message) ;
    if (status < 0)
    {
	if (ff != NULL) fclose (ff) ;
	exit (0) ;
    }
}

/**************************************************************************\
SuiteSparse_long cholmod_nested_dissection	// returns # of components //
(
    // ---- input ---- //
    cholmod_sparse *A,  // matrix to order
    int *fset,		    // subset of 0:(A->ncol)-1
    size_t fsize,	    // size of fset

    // ---- output --- //
    int *Perm,		// size A->nrow, output permutation 
    int *CParent,	// size A->nrow.  On output, CParent [c] is the parent
                    // of component c, or EMPTY if c is a root, and where
                    // c is in the range 0 to # of components minus 1
    int *Cmember,	// size A->nrow.  Cmember [j] = c if node j of A is
                    // in component c 
                        
    // --------------- //
    cholmod_common *Common
) ;

int cholmod_metis
(
    // ---- input ---- //
    cholmod_sparse *A,	// matrix to order
    int *fset,		// subset of 0:(A->ncol)-1
    size_t fsize,	// size of fset
    int postorder,	// if TRUE, follow with etree or coletree postorder
    /* ---- output --- //
    int *Perm,		// size A->nrow, output permutation
    // --------------- //
    cholmod_common *Common
) ;
\**************************************************************************/

int main(int argc, char *argv[])
{
    // ---------------------------------------------------------------------- //
    // declaration of variables                                               //
    // ---------------------------------------------------------------------- //

    cholmod_sparse *A; // cholmod_sparse matrix
    FILE *f ;
    cholmod_common Common, *cm ;
    cholmod_factor *L ;

    int *fset, *perm, *iperm, *CParent, *Cmember, *partition ; 
    int nrow, orderingMethod, postorder, compress ;
    size_t fsize ;

    // ---------------------------------------------------------------------- //
    // terminal input processing                                              //
    // ---------------------------------------------------------------------- //

    ff = NULL ;
    if (argc == 4)
    {
        // get the file containing the input matrix
	    if ((f = fopen (argv [1], "r")) == NULL)
	    {
	        my_handler (CHOLMOD_INVALID, __FILE__, __LINE__,
		        "unable to open file") ;
	    }
	    ff = f ;
        fsize = atoi (argv [2]) ;
        orderingMethod = atoi (argv [3]) ;
    }   
    else
    {
        printf ("*************************************************************\n") ;
        printf ("Usage : testCholmod matrixFile nbrow orderingMethod\n"
                "you can find different matrix files in the repository /Matrix.\n"
                "nbrow : please indicate the size of the square matrix.\n"
                "orderingMethod possible values : \n"
                "---> 0 : nested dissection\n"
                "---> 1 : metis\n"
                "---> 2 : bisect\n"
                "\nexample : ./testCholmod MatrixTri/matrix15.tri 15 0\n") ;
        printf ("*************************************************************\n") ;
        exit (1) ;
    }

    printf ("*************************************************************\n") ;
    printf (" Testing program of Cholmod Partitionning\n\n") ;

    // ---------------------------------------------------------------------- //
    // start CHOLMOD and set parameters                                       //
    // ---------------------------------------------------------------------- //

    cm = &Common ;
    cholmod_start (cm) ;
    // CHOLMOD_FUNCTION_DEFAULTS ;
    cm->error_handler = my_handler ;

    // ---------------------------------------------------------------------- //
    // read in a matrix                                                       //
    // ---------------------------------------------------------------------- //

    A = cholmod_read_sparse (f, cm) ;    // read the matrix
    if (ff != NULL)
    {
        fclose (ff) ;
        ff = NULL ;
    }

    // ---------------------------------------------------------------------- //
    // allocate memory                                                        //
    // ---------------------------------------------------------------------- //

    nrow = A->nrow ;
    fset = malloc (fsize * sizeof(int)) ;
    perm = malloc (nrow * sizeof(int)) ;
    iperm = malloc (nrow * sizeof(int)) ;
    CParent = malloc (nrow * sizeof(int)) ;
    Cmember = malloc (nrow * sizeof(int)) ;
    partition = malloc (nrow * sizeof(int)) ;

    for (int k = 0; k < fsize; k++)
    {
        fset [k] = k ;
    }

    // ---------------------------------------------------------------------- //
    // execution of the main functions                                        //
    // ---------------------------------------------------------------------- //

    // cm->method [3] ;// NESDIS
    L = cholmod_analyze (A, cm) ;        // analyse
    cholmod_factorize (A, L, cm) ;       // factorize

    if (orderingMethod == NESTED_DISSECTION)
    {
        cholmod_nested_dissection (A, fset, fsize, perm, CParent, Cmember, cm) ;
    }
    else if (orderingMethod == METIS)
    {
        postorder = 1 ;
        cholmod_metis (A, fset, fsize, postorder, perm, cm) ;
    }
    else if (orderingMethod == BISECT)
    {
        compress = 0 ; // if 1, compress the graph --> only use it for big graph
        cholmod_bisect (A, fset, fsize, compress, partition, cm) ;
    }
    
    if (orderingMethod == NESTED_DISSECTION || orderingMethod == METIS)
    {    
        for (int k = 0; k < fsize; k++)
        {
            iperm [perm [k]] = k;
        }
    }

    // ---------------------------------------------------------------------- //
    // display the results                                                    //
    // ---------------------------------------------------------------------- //

    printf ("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n") ;

    cholmod_print_sparse (A, "A", cm) ;  // print the matrix
    cholmod_print_factor (L, "L", cm) ;  // print the factorization

    printf ("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n") ;

    // print the permutation
    if (orderingMethod == NESTED_DISSECTION || orderingMethod == METIS)
    {
        printf (" perm = [ ") ;
        for (int k = 0; k < fsize; k++)
        {
            printf ("%d ", perm [k]) ;
        }
        printf ("]\n");
    }

    // print the inverse of the permutation (graph of PtAP)
    if (orderingMethod == NESTED_DISSECTION || orderingMethod == METIS)
    {
        printf ("iperm = [ ") ;
        for (int k = 0; k < fsize; k++)
        {
            printf ("%d ", iperm [k]) ;
        }
        printf ("]\n");
    }

    // print the Cmember vector
    if (orderingMethod == NESTED_DISSECTION)
    {
        printf ("Cmember = [ ") ;
        for (int k = 0; k < fsize; k++)
        {
            printf ("%d ", Cmember [k]) ;
        }
        printf ("]\n");
    }

    // print the elimination tree
    if (orderingMethod == NESTED_DISSECTION)
    {
        printf ("CParent = [ ") ;
        for (int k = 0; k < fsize; k++)
        {
            printf ("%d ", CParent [k]) ;
        }
        printf ("]\n");
    }

    // print the partitions
    if (orderingMethod == BISECT)
    {
        printf ("partition = [ ") ;
        for (int k = 0; k < fsize; k++)
        {
            printf ("%d ", partition [k]) ;
        }
        printf ("]\n");
    }

    printf ("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n") ;

    // ---------------------------------------------------------------------- //
    // free memory and ending the function                                    //
    // ---------------------------------------------------------------------- //

    cholmod_free_factor (&L, cm) ; // free the factorization
    cholmod_free_sparse (&A, cm) ; // free the matrix
    cholmod_finish (cm) ;
    free (fset) ; free (perm) ; free (iperm) ; 
    free (CParent) ; free (Cmember) ; free (partition) ;

    printf ("\nTesting program complete, return\n");
    printf ("*************************************************************\n") ;
    return (0) ;  // Successfull exit 
}