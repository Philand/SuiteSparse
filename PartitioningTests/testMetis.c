#include <stdio.h>
#include <stdlib.h>
#include <metis.h>

int main(int argc, char *argv[])
{
    // ---------------------------------------------------------------------- //
    // declaration of variables                                               //
    // ---------------------------------------------------------------------- //

	idx_t nvtxs;		// number of vertices
	idx_t *xadj ;		// accumulation of the number of vertex neighbors 
	idx_t *adjncy ;		// vertex neighbors in ascending order  
	idx_t *vwgt;		// weight of the edges
	idx_t npes;			// number of partitions
	idx_t options[METIS_NOPTIONS];
	idx_t *perm;		// permutation of the graph (where goes the vertices)
	idx_t *invperm;		// inverse permutation (mapping of the permuted graph)
	idx_t *sizes;		// Sizes of the partitions & separators
	int result;			// result of the metis function (1 = METIS_OK)

    // ---------------------------------------------------------------------- //
    // different examples	                                                  //
    // ---------------------------------------------------------------------- //

	idx_t ex0_xadj[16] = {0, 2, 5, 8, 11, 13, 16, 20, 24, 28, 31, 33, 36, 39, 42, 44} ;
	idx_t ex0_adjncy[44] = { 1, 5, 0, 2, 6, 1, 3, 7, 2, 4, 8, 3, 9, 0, 6, 10, 1, 5, 7, 11, 2, 6, 8, 12, 3, 7, 9, 13, 4, 8, 14, 5, 11, 6, 10, 12, 7, 11, 13, 8, 12, 14, 9, 13} ;

	idx_t ex1_xadj[36] = {0, 2, 5, 8, 11, 13, 16, 20, 24, 28, 31, 34, 38, 42, 46, 49, 52, 56, 60, 64, 67, 70, 73, 77, 81, 85, 88, 92, 96, 100, 103, 105, 108, 111, 114, 116} ;
	idx_t ex1_adjncy[116] = {1, 5, 0, 2, 6, 1, 3, 7, 2, 4, 8, 3, 9, 0, 6, 10, 1, 5, 7, 11, 2, 6, 8, 12, 3, 7, 9, 13, 4, 8, 14, 5, 11, 15, 6, 10, 12, 16, 7, 11, 13, 17, 8, 12, 14, 18, 9, 13, 19, 10, 16, 20, 11, 15, 17, 21, 12, 16, 18, 22, 13, 17, 19, 23, 14, 18, 24, 15, 21, 25, 16, 20, 22, 26, 17, 21, 23, 27, 18, 22, 24, 28, 19, 23, 29, 20, 26, 30, 21, 25, 27, 31, 22, 26, 28, 32, 23, 27, 29, 33, 24, 28, 34, 25, 31, 26, 30, 32, 27, 31, 33, 28, 32, 34, 29, 33} ;

    // ---------------------------------------------------------------------- //
    // terminal input processing                                              //
    // ---------------------------------------------------------------------- //

	if (argc == 3)
	{
		npes = atoi(argv[1]); // number of partitions
		nvtxs = 35;
		METIS_SetDefaultOptions(options);
	}
	else
	{
		printf("Usage : testPartitionning nb_parts num_graph\n");
		printf ("nb_parts : choose a power of 2 (1, 2, 4, 8, 16, ...") ;
		printf ("graph possible values : \n") ;
	    printf ("---> 0 : 15 vertices & 22 edges\n") ;
        printf ("---> 1 : 35 vertices & 58 edges\n") ;
		exit(1);
	}

    printf ("*************************************************************\n") ;
    printf (" Testing METIS API functions\n\n") ;

    // ---------------------------------------------------------------------- //
    // allocate memory                                                        //
    // ---------------------------------------------------------------------- //
	
	perm = malloc(nvtxs*sizeof(long int));
	invperm = malloc(nvtxs*sizeof(long int));
	sizes = malloc(2*npes*sizeof(long int));

    // ---------------------------------------------------------------------- //
    // execution of the main functions                                        //
    // ---------------------------------------------------------------------- //
	
	if (npes > 0)
	{
		result = METIS_NodeNDP(nvtxs, &xadj[0], &adjncy[0], NULL, npes, NULL, &perm[0], &invperm[0], &sizes[0]);
	}
	else
	{
		result = METIS_NodeND(&nvtxs, &xadj[0], &adjncy[0], NULL, NULL, &perm[0], &invperm[0]);
	}
	
    // ---------------------------------------------------------------------- //
    // display the results                                                    //
    // ---------------------------------------------------------------------- //

	printf ("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n") ;

	printf("--> Return : %i\n", result);

	printf ("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n") ;


	// print the permutation
	printf("-->  perm = [ ");
	for (int k = 0; k < nvtxs; k++)
	{
		printf("%ld ", perm[k]);
	}
	printf("]\n");

	// print the inverse of the permutation (graph of PtAP)
	printf("--> iperm = [ ");
	for (int k = 0; k < nvtxs; k++)
	{
		printf("%ld ", invperm[k]);
	}
	printf("]\n");

	// print the sizes of the partitions & separators
	if (npes > 0)
	{	
		printf("--> psizes = [ ");
		for (int k = 0; k < npes*2-1; k++)
		{
			printf("%ld ", sizes[k]);
		}
		printf("]\n");
	}

	printf ("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n") ;

    // ---------------------------------------------------------------------- //
    // free memory and ending the function                                    //
    // ---------------------------------------------------------------------- //
	
	free(perm) ; free(invperm) ; free(sizes) ;
	free (xadj) ; free (adjncy) ;

    printf ("\nTesting program complete, return\n");
    printf ("*************************************************************\n") ;
	return 0;
}