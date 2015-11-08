#include <stdio.h>
#include <R.h>
#include <Rdefines.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>


#include "SERNLIB/FastSERN.h"
#include "SERNLIB/edgeprobfuncs.h"

#include <time.h>
#include <stdarg.h>
#include "SERNLIB/options.h"

#define NUM_ELEMS(X)           (sizeof(X)/sizeof(*(X)))

#define PROBABILITY_FUN_rhs 0
#define S_rhs               1
#define Q_rhs               2
#define N_rhs               3

#define DISTANCE_FUN_rhs    4
#define REGION_SHAPE_rhs    5
#define REGION_GEOMETRY_rhs 6


// we should have defaults for all these except seed
#define FIRST_DEFAULT       7
#define CONNECTED_rhs       7
#define M_rhs               8
#define THREADS_rhs         9
#define ALGORITHM_rhs       10
#define BUFFER_SIZE_rhs     11
#define SEED_rhs            12


#define CONNECTED_default   0
#define M_default           1
#define THREADS_default     1
#define ALGORITHM__default  0
#define BUFFER_SIZE_default 100000


void errIdAndTxt(const char * filename,
                 uint32_t line,
                 const char *fmt, ...)
{
    va_list argptr;
    Rprintf("\n%s line %d : ", filename, line);
    va_start(argptr,fmt);
    Rprintf(fmt, argptr);
    va_end(argptr);

}


SEXP generate(SEXP args)
{

    int  i;
    double *x, *y, *weights;

    NodeList nodes;
    EdgeList edges;
    Options options;

    GeometryStruct *geometry;
    PolygonStruct *polygon;

    VectorStruct unit_square[2] = {{0.0, 0.0}, {1.0, 1.0}};
    VectorStruct *vector;

    SEXP vec_id, vec_x, vec_y, vec_from, vec_to, vec_weights;
    SEXP list_edge_names, list_node_names;
    SEXP list_edge, list_node;
    SEXP cls_edge, cls_node;
    SEXP edge_row_names, node_row_names;
    SEXP ret;

    char *edge_names[3] = {"from", "to", "weight"};
    char *node_names[3] = {"id", "x", "y"};

    uint64_t N;
    int *from, *to, *id;

    // before we start clear all our structures
    memset(&nodes, 0, sizeof(NodeList));
    memset(&edges, 0, sizeof(EdgeList));
    memset(&options, 0, sizeof(Options));


    /* first tell conSERN what memory allocators */
    /* to use for data returned to us            */
    options.realloc = realloc;
    options.calloc = calloc;

    /* and how to display error and exit */
    options.errIdAndTxt = errIdAndTxt;




    /* check what optional outputs are required */
    //    /* only allocate space for the "distances" if required */
    //    if (nlhs >= 6) options.weights_enabled = 1;
    options.weights_enabled = 1;


    //    /* only allocate space for node component identifiers if required */
    //    if(nlhs == 7) options.components_enabled = 1;
    options.components_enabled = 0;



    //    /* which function pointer to use for calculating probabilities */
    //    if ((uint32_t)mxGetScalar(prhs[PROBABILITY_FUN_rhs]) >=
    //        NUM_ELEMS(probabilityFunctions))
    //    {
    //        mexErrMsgTxt("unimplemented probability function selected");
    //    }
    //
    //    options.probGivenDistance =
    //    probabilityFunctions[(uint32_t)mxGetScalar(prhs[PROBABILITY_FUN_rhs])];
    //
    options.probGivenDistance =  probabilityFunctions[0 /*waxman*/];

    //
    //    // there is a possibility that we have more than one shape parameter */
    //    data = mxDuplicateArray(prhs[S_rhs]);
    //
    //    // if we only have one use it
    //    if (mxGetN(data) == 1)
    //        options.s1 = mxGetScalar(prhs[S_rhs]);
    //    else // we have more than one
    //    {
    //        // two is ok
    //        if (mxGetN(data) == 2)
    //        {
    //            s = (double *) mxGetData(data);
    //            options.s1 = s[1];
    //            options.s2 = s[2];
    //        }
    //        else
    //        {
    //            // but any other number is no good
    //            mexErrMsgTxt("Only two values for s allowed");
    //        }
    //
    //    }
    options.s1 = REAL(args)[0];


    //    /* the thinning parameter */
    //    options.q = mxGetScalar(prhs[Q_rhs]);
    options.q = REAL(args)[1];



    //    /* number of nodes in the  graph */
    //    options.N = (uint32_t) mxGetScalar(prhs[N_rhs]);
    options.N = round(REAL(args)[2]);
    N = options.N;



    //    /* which function pointer to use for calculating distances */
    //    if ((uint32_t) mxGetScalar(prhs[DISTANCE_FUN_rhs]) >=
    //        NUM_ELEMS(distanceFunctions))
    //    {
    //        mexErrMsgTxt("unimplemented metric selected");
    //    }
    //
    //    /* distance given a diatace in x and one in y */
    //    options.distance =
    //    distanceFunctions[(uint32_t) mxGetScalar(prhs[DISTANCE_FUN_rhs])];
    options.distance = distanceFunctions[0 /* euclidean */];



    //
    //    /* default is to leave the graph disconnected */
    //    options.connected = (nrhs > CONNECTED_rhs) ?
    //    (uint32_t) mxGetScalar(prhs[CONNECTED_rhs]) : CONNECTED_default ;
    options.connected = 0;


    //    /* dimension of array of buckets */
    //    options.M = (nrhs > M_rhs) ?
    //    (uint32_t) mxGetScalar(prhs[M_rhs]) : M_default;
    options.M = 10;


    //    /* number of threads to use */
    //    options.ThreadCount = (nrhs > THREADS_rhs) ?
    //    (uint32_t) mxGetScalar(prhs[THREADS_rhs]) : THREADS_default;
    options.ThreadCount = 1;


    //    /* fast algorithm = 0 N^2 algorithm = 1 so default is fast */
    //    options.algorithm = (nrhs > ALGORITHM_rhs) ?
    //    (uint32_t) mxGetScalar(prhs[ALGORITHM_rhs]) : ALGORITHM__default;
    options.algorithm = 0;


    //    /* number of edges in buffer for each thread */
    //    options.BufferSize = (nrhs > BUFFER_SIZE_rhs)  ?
    //    (uint32_t) mxGetScalar(prhs[BUFFER_SIZE_rhs]) : BUFFER_SIZE_default;
    options.BufferSize =  BUFFER_SIZE_default;


    //    /* initialize random number generator */
    //    options.seedval = (nrhs > SEED_rhs) ?
    //    (uint32_t) mxGetScalar(prhs[SEED_rhs]) : (uint32_t) time(0);
    // read in the random number seed

    GetRNGstate();
    options.seedval = (uint32_t)(unif_rand() * (2^32 - 1));


    //
    //
    //    /* we can only handle the shape and the geometry once */
    //    /*  we have the options particularly M                */
    //    data = mxDuplicateArray(prhs[REGION_GEOMETRY_rhs]);
    //    vector  = (VectorStruct *) mxGetData(data);
    //    assert(mxGetM(data) == 2); // we already checked this
    //
    //    polygon = polygonNew();
    //
    //    for (i = 0;i <  mxGetN(data); i++)
    //        polygonAppend(&options, polygon, vector + i);
    //
    //    geometry = geometryGenerate( &options,
    //                                (GeometryType)
    //                                mxGetScalar(prhs[REGION_SHAPE_rhs]),
    //                                polygon);


    // force the unit square
    polygon = polygonNew();
    vector = unit_square;
    for (i = 0; i <  2; i++)
        polygonAppend(&options, polygon, vector + i);
    geometry = geometryGenerate( &options, rectangle, polygon);

    /* create the graph */
    // TODO function to set initial size close to that needed
    edges.growth = 4 * options.N;
    GenSERN(&nodes, &edges, &options, geometry);
    polygonFree(polygon);

    // creating an integer vector:
    PROTECT(vec_x = NEW_NUMERIC(N));
    x = NUMERIC_POINTER(vec_x);
    for (i = 0; i < N; i++)
        x[i] = nodes.x[i];
    free(nodes.x);

    // ... and a vector of real numbers:
    PROTECT(vec_y = NEW_NUMERIC(N));
    y = NUMERIC_POINTER(vec_y);
    for (i = 0; i < N; i++)
        y[i] = nodes.y[i];
    free(nodes.y);

    PROTECT(vec_id = NEW_INTEGER(N));
    id = INTEGER_POINTER(vec_id);
    for (i = 0; i < N; i++)
      id[i] = i+1;

    PROTECT(vec_from = NEW_INTEGER(edges.count));
    from = INTEGER_POINTER(vec_from);
    memcpy(from, edges.from, edges.count * sizeof(int));
    free(edges.from);

    PROTECT(vec_to = NEW_INTEGER(edges.count));
    to = INTEGER_POINTER(vec_to);
    memcpy(to, edges.to, edges.count * sizeof(int));
    free(edges.to);

    PROTECT(vec_weights = NEW_NUMERIC(edges.count));
    weights = NUMERIC_POINTER(vec_weights);

    memcpy(weights, edges.weight, edges.count * sizeof(double));
    free(edges.weight);

    // write out the random number seed
    PutRNGstate();

    // Creating a character string vector
    // of the "names" attribute of the
    // objects in out list:
    PROTECT(list_edge_names = allocVector(STRSXP,3));
    for(i = 0; i < 3; i++)
        SET_STRING_ELT(list_edge_names,i,mkChar(edge_names[i]));
    // Creating a list with 6 vector elements:
    PROTECT(list_edge = allocVector(VECSXP, 3));

    SET_VECTOR_ELT(list_edge, 0, vec_from);
    SET_VECTOR_ELT(list_edge, 1, vec_to);
    SET_VECTOR_ELT(list_edge, 2, vec_weights);


    PROTECT(list_node_names = allocVector(STRSXP,3));
    for(i = 0; i < 3; i++)
      SET_STRING_ELT(list_node_names,i,mkChar(node_names[i]));
    // Creating a list with 6 vector elements:
    PROTECT(list_node = allocVector(VECSXP, 3));

    SET_VECTOR_ELT(list_node, 0, vec_id);
    SET_VECTOR_ELT(list_node, 1, vec_x);
    SET_VECTOR_ELT(list_node, 2, vec_y);

    // and attaching the vector names:
    setAttrib(list_edge, R_NamesSymbol, list_edge_names);
    setAttrib(list_node, R_NamesSymbol, list_node_names);

    // make edge list into a dataframe
    PROTECT(cls_edge = allocVector(STRSXP, 1)); // class attribute
    SET_STRING_ELT(cls_edge, 0, mkChar("data.frame"));
    classgets(list_edge, cls_edge);

    // short hand way of adding automatic row names
    PROTECT(edge_row_names = allocVector(INTSXP, 2));
    INTEGER(edge_row_names)[0] = NA_INTEGER;
    INTEGER(edge_row_names)[1] = edges.count;
    setAttrib(list_edge, R_RowNamesSymbol, edge_row_names);


    // make node list into a dataframe
    PROTECT(cls_node = allocVector(STRSXP, 1)); // class attribute
    SET_STRING_ELT(cls_node, 0, mkChar("data.frame"));
    classgets(list_node, cls_node);

    // short hand way of adding automatic row names
    PROTECT(node_row_names = allocVector(INTSXP, 2));
    INTEGER(node_row_names)[0] = NA_INTEGER;
    INTEGER(node_row_names)[1] = N;
    setAttrib(list_node, R_RowNamesSymbol, node_row_names);

    // we now have two dataframes list_edge and list_node
    // we want to return them both in a list. list_edge in the
    // first element of the list and list_nodes in the second
    // element
    PROTECT(ret = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ret, 0, list_edge);
    SET_VECTOR_ELT(ret, 1, list_node);
    UNPROTECT(15);
    return ret;
}

