#ifndef _BIGLOBAL_UTIL_H_
#define _BIGLOBAL_UTIL_H_

// -----------------------------------------------------------------------------
// Helpers for allocating vectors and matrices
// -----------------------------------------------------------------------------
void *vec_alloc(int, int, size_t);
void free_vec(void *, int, size_t);
void **mat2_alloc(int,int,int,size_t);
void free_mat2(void **,int,size_t);

double *dvec_alloc(int, int);
double *dvec_view_alloc(int ,int ,int );

double **dmat_alloc(int, int,int);
double **dmat_view_alloc(int , int , int ,int ,int );

double ***dmat3d_alloc(int, int, int, int);
double **create_dmat_from_dvec(double *,int,int,int);

int    *ivec_alloc(int, int);
int    **imat_alloc(int, int, int);
int    ***imat3d_alloc(int, int, int, int);

void free_dvec(double *,int);
void free_dmat(double **, int);
void free_dmat_view(double **,int ,int ,int );
void free_dmat3d(double ***,int);

void free_ivec(int *,  int);
void free_imat(int **,  int);
void free_imat3d(int ***, int);

void ***mat3_view_alloc( int , int , int ,int , int , int ,int , size_t);

#endif // _BIGLOBAL_UTIL_H_
