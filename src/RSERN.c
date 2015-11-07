#include <stdio.h>
#include <R.h> 
#include <Rdefines.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
//#include "rand64.h"
//#include "eric_rand.h"
//#include "geom_rand.h" using R random functions atm


typedef struct 
{
    double *x;
    double *y;
    uint64_t count;
    uint64_t start;
    uint64_t allocated;
    uint64_t growth;
    uint32_t i;
    uint32_t j;
} BucketStruct;


typedef struct 
{
    int *from; // R has no unsigned integer types
    int *to;
    double *weight;
    uint32_t weights_enabled;
    uint64_t count;
    uint64_t allocated;
    uint64_t growth;
} EdgeList;



inline uint32_t geom_rand2(double lambda)
{
    /* should probably check that 0<=p<1 */
	return	(fabs(lambda) < DBL_EPSILON) ? 	UINT32_MAX: 
											1 + floor(log(unif_rand()) / lambda);
}

inline uint32_t AddNode(BucketStruct *l, double x, double y)
{
    if (l->count >= l->allocated) 
    {
        l->allocated += l->growth;
        /* realloc on a null pointer acts like malloc */
        l->x = realloc(l->x, sizeof(double) * l->allocated);
        l->y = realloc(l->y, sizeof(double) * l->allocated);
        /* after the initial allocation add 10% each time */
        l->growth = 1 + l->allocated / 10; 
    }
    
    l->x[l->count] = x;
    l->y[l->count] = y;
    (l->count)++;
    
    return l->allocated;  
}

uint32_t AddEdge(EdgeList *l, int from, int to, double weight)
{
    if (l->count >= l->allocated) 
    {
        l->allocated += l->growth;
        /* realloc on a null pointer acts like malloc */
        l->from = realloc(l->from, sizeof(int) * l->allocated);
        l->to = realloc(l->to, sizeof(int) * l->allocated);
        
        /* only re-allocate space for the "distances" if required */
        if (l->weights_enabled) 
        {
            l->weight = realloc(l->weight,sizeof(double) * l->allocated);
        }
        
        l->growth =  l->allocated/10;
        
    }
    
    l->from[l->count] = from;
    l->to[l->count] = to;
    if (l->weights_enabled) l->weight[l->count] = weight;
    (l->count)++;
    
    return l->allocated;
}


void GenerateNodes(uint32_t M, uint32_t N, BucketStruct *buckets)
{ 
    
    double x_coord, y_coord;
    uint32_t bucket_x = 0, bucket_y = 0, node_number;
    
    if (M>1) 
    {
        /* choose node locations, and round them off */
        for (node_number = 0; node_number < N; node_number++) 
        {
            x_coord = unif_rand();
            y_coord = unif_rand();
            
            bucket_x = (uint32_t) ((double)M * x_coord);
            bucket_y = (uint32_t) ((double)M * y_coord);
            
            AddNode(&(buckets[bucket_y + M * bucket_x]), x_coord, y_coord);
        }
    }
    else
    {
        for (node_number = 0; node_number < N; node_number++) 
        {
            x_coord = unif_rand();
            y_coord = unif_rand();
            
            AddNode(&(buckets[bucket_y + M * bucket_x]), x_coord, y_coord);
        }
    }
}


BucketStruct *GenerateBuckets(uint32_t M, uint32_t N, double *x, double *y)
{
    uint32_t i;
    
    /* calloc used to zero all fields */
    BucketStruct * buckets = calloc(M * M, sizeof(BucketStruct));
    
    for(i = 0; i < M * M; i++)
    {
        /* the first bucket data is stored in the provided memory */
        /* this means we dont need to copy node coordinates into  */
        /* the provided memory if there is only one bucket        */ 
        if (i == 0)
        {
            buckets[0].start = 1;
            buckets[0].x = x;
            buckets[0].y = y;
            buckets[0].allocated = N;
            continue;        
        }
        /* for other buckets create our own memory */    
        /* set initial growth to be twice expected */
        //buckets[i].growth = ceil(2 * N / (M * M));
        buckets[i].growth = ceil((1 + 4 * M /sqrt(N)) * N/(M * M));
        buckets[i].start = 1;
        buckets[i].i = i % M;
        buckets[i].j = i / M;
    }
    
    GenerateNodes(M, N, buckets);
    
    /* calculate the starting node number for each bucket */
    for(i = 1; i < M * M; i++) 
        buckets[i].start = buckets[i - 1].start + buckets[i - 1].count;
    
    return buckets;
}



/* relies on the surface being a unit square */
double *CreateQ(uint32_t M, double s)
{
    double *Q;
    uint64_t i, j;
    double distance;
    
    /* calloc initialises with zeros */
    double *t = calloc(M, sizeof(double));
    
    for (i = 1; i < M; i++) t[i] = i - 1;    
    
    Q = malloc((size_t) sizeof(double) * M * M);
    
    for(i = 0; i < M; i++)
    {
        for(j = 0; j < M; j++)
        {
            distance = sqrt(t[i] * t[i] + t[j] * t[j]) / M ;
            Q[i * M + j] =  exp(-s*distance);
        }
    }
    
    free(t);
    return Q;
}


void StoreNodes(BucketStruct *buckets, uint32_t M, double *x, double *y)
{
    uint32_t i;
    
    for  (i = 1; i < M * M; i++)
    { 
        memcpy(x + buckets[i].start - 1,  
               buckets[i].x, sizeof(double )* buckets[i].count);
        free(buckets[i].x);
        buckets[i].x = x + buckets[i].start - 1;
    }
    
    for  (i = 1; i < M * M; i++)
    {
        memcpy(y + buckets[i].start - 1, 
               buckets[i].y, sizeof(double)* buckets[i].count);
        free(buckets[i].y);
        buckets[i].y = y + buckets[i].start - 1;
    }
}


/* the actual code to generate the Waxman graph */
int WaxmanGen(double s, double beta, int64_t N, int M,
              double* x, double* y,
              EdgeList* edges)
{
    
    int64_t i, j, k, S;
    uint32_t bucket_a, bucket_b;
    double p, lambda, distance, x_diff, y_diff, *Q;
    BucketStruct bucket_A, bucket_B, *buckets;
    
    
    Q = CreateQ(M, s);
    
    buckets = GenerateBuckets(M, N, x, y);
	// Rprintf("The number of buckets is %ld\n",M);
    
    
    StoreNodes(buckets, M, x, y);
    
    
    // Generate intra bucket edges in bucket_a 
    lambda = log(1 - beta);
    
  
    
    for (bucket_a = 0; bucket_a < M * M; bucket_a++)
    {
        k = -1;
        bucket_A = buckets[bucket_a];
        
        S = (bucket_A.count * (bucket_A.count - 1))/2;
        
        // find the next link that might exist  
        while ((k += geom_rand2(lambda)) < S) 
        {
            //Rprintf("\n k = %ld", k);
            // can we do an integer square root here ?
            j = 1 + ((((int64_t)sqrt(8 * k + 1)) - 1) / 2); 
            i = k - j * (j - 1) / 2;  
            
            x_diff = bucket_A.x[i] - bucket_A.x[j];
            y_diff = bucket_A.y[i] - bucket_A.y[j];
            distance = sqrt(x_diff * x_diff + y_diff * y_diff);
            
            if (unif_rand() < exp(-s * distance)) 
                AddEdge(edges, bucket_A.start + i, 
                        bucket_A.start + j, distance);
        } 
    } 
    
    
    
    // Generate inter bucket edges between bucket_a and bucket_b 
    for (bucket_a = 0; bucket_a < M * M; bucket_a++) 
    {
        bucket_A = buckets[bucket_a];
        // Rprintf("bucket_a is %d\n",bucket_a);
        for (bucket_b = bucket_a + 1; bucket_b < M * M; bucket_b++) 
        {
	        // Rprintf("bucket_b is %d\n",bucket_b);
            k = -1;
            bucket_B = buckets[bucket_b];
            
            p =  Q[abs(bucket_A.i - bucket_B.i) * M + 
                   abs(bucket_A.j - bucket_B.j)]; 
			
            lambda = log( 1 - p * beta); 
            
            S = bucket_A.count * bucket_B.count;
            
            // find the next link that might exist
            while ((k += geom_rand2(lambda)) < S) 
            {    
                
                i = k % bucket_A.count;
                j = k / bucket_A.count;
                
                x_diff = bucket_A.x[i] - bucket_B.x[j];
                y_diff = bucket_A.y[i] - bucket_B.y[j];
                distance = sqrt(x_diff * x_diff + y_diff * y_diff);
                
                if (unif_rand() * p < exp(-s * distance)) 
                    AddEdge(edges, bucket_A.start + i, 
                            bucket_B.start + j, distance);
            } 
        }
    }
     
    
    /* cleanup */
    free(buckets);
    free(Q);
    
    return(edges->allocated);
    
} /* WaxmanGen */


SEXP old_generate(SEXP args) 
    {
      int  i;
      double *x, *y, *weights;
      double beta, s;
      EdgeList edges = {NULL, NULL, NULL, 0, 0, 0, 0};
      SEXP vec_x, vec_y, vec_weights, vec_from, vec_to, list, list_names;
      char *names[5] = {"X", "Y", "from", "to", "weight"};
      uint64_t N;
      int *from, *to;
      
      // get the arguments
      s  = REAL(args)[0];
      beta = REAL(args)[1];
      N = round(REAL(args)[2]);
      
      // read in the random number seed
      GetRNGstate();
      
      
      // creating an integer vector: 
      PROTECT(vec_x = NEW_NUMERIC(N)); 
      x = NUMERIC_POINTER(vec_x);
      
      
      // ... and a vector of real numbers: 
      PROTECT(vec_y = NEW_NUMERIC(N)); 
      y = NUMERIC_POINTER(vec_y);
      
      edges.growth = 1.5 * N; // TODO don't ask this is a kludge atm
      edges.weights_enabled = 1;
      
      // create the graph TODO calculate M rather than hardwire to 10
      WaxmanGen(s, beta, N, 10, x, y, &edges);
      
      
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
      PROTECT(list_names = allocVector(STRSXP,5));
      for(i = 0; i < 5; i++)
        SET_STRING_ELT(list_names,i,mkChar(names[i]));
      // Creating a list with 5 vector elements: 
      PROTECT(list = allocVector(VECSXP, 5));
      
      SET_VECTOR_ELT(list, 0, vec_x);
      SET_VECTOR_ELT(list, 1, vec_y);
      SET_VECTOR_ELT(list, 2, vec_from);
      SET_VECTOR_ELT(list, 3, vec_to);
      SET_VECTOR_ELT(list, 4, vec_weights);
      
      // and attaching the vector names: 
      setAttrib(list, R_NamesSymbol, list_names); 
      UNPROTECT(7);
      return list;
    }

SEXP generate_1(SEXP args) 
{
  int  i;
  double *x, *y, *weights;
  double beta, s;
  EdgeList edges = {NULL, NULL, NULL, 0, 0, 0, 0};
  SEXP vec_id, vec_x, vec_y, vec_weights, vec_from, vec_to, list, list_names;
  char *names[6] = {"from", "to", "weight", "id", """x", "y"};
  uint64_t N;
  int *from, *to, *id;
  
  // get the arguments
  s  = REAL(args)[0];
  beta = REAL(args)[1];
  N = round(REAL(args)[2]);
  
  // read in the random number seed
  GetRNGstate();
  
  
  // creating an integer vector: 
  PROTECT(vec_x = NEW_NUMERIC(N)); 
  x = NUMERIC_POINTER(vec_x);
  
  
  // ... and a vector of real numbers: 
  PROTECT(vec_y = NEW_NUMERIC(N)); 
  y = NUMERIC_POINTER(vec_y);
  
  
  PROTECT(vec_id = NEW_INTEGER(N)); 
  id = INTEGER_POINTER(vec_id);
  for (i = 0; i < N; i++)
    id[i] = i+1;
  
  edges.growth = 1.5 * N; // TODO don't ask this is a kludge atm
  edges.weights_enabled = 1;
  
  // create the graph TODO calculate M rather than hardwire to 10
  WaxmanGen(s, beta, N, 10, x, y, &edges);
  
  
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
  PROTECT(list_names = allocVector(STRSXP,6));
  for(i = 0; i < 6; i++)
    SET_STRING_ELT(list_names,i,mkChar(names[i]));
  // Creating a list with 6 vector elements: 
  PROTECT(list = allocVector(VECSXP, 6));
  
  SET_VECTOR_ELT(list, 0, vec_from);
  SET_VECTOR_ELT(list, 1, vec_to);
  SET_VECTOR_ELT(list, 2, vec_weights);
  
  SET_VECTOR_ELT(list, 3, vec_id);
  SET_VECTOR_ELT(list, 4, vec_x);
  SET_VECTOR_ELT(list, 5, vec_y);
  
  // and attaching the vector names: 
  setAttrib(list, R_NamesSymbol, list_names); 
  UNPROTECT(8);
  return list;
}
SEXP generate(SEXP args) 
{
    int  i;
    double *x, *y, *weights;
    double beta, s;
    EdgeList edges = {NULL, NULL, NULL, 0, 0, 0, 0};
    
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
    
    // get the arguments
    s  = REAL(args)[0];
    beta = REAL(args)[1];
    N = round(REAL(args)[2]);
    
    // read in the random number seed
    GetRNGstate();
    
    
    // creating an integer vector: 
    PROTECT(vec_x = NEW_NUMERIC(N)); 
    x = NUMERIC_POINTER(vec_x);
    
    
    // ... and a vector of real numbers: 
    PROTECT(vec_y = NEW_NUMERIC(N)); 
    y = NUMERIC_POINTER(vec_y);
    
    
    PROTECT(vec_id = NEW_INTEGER(N)); 
    id = INTEGER_POINTER(vec_id);
    for (i = 0; i < N; i++)
      id[i] = i+1;
  
    edges.growth = 1.5 * N; // TODO don't ask this is a kludge atm
    edges.weights_enabled = 1;

    // create the graph TODO calculate M rather than hardwire to 10
    WaxmanGen(s, beta, N, 10, x, y, &edges);
    
    
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