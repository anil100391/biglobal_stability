#ifndef _BIGLOBAL_DATA_H_
#define _BIGLOBAL_DATA_H_
/*
 *  This file contains definition of data structures and function prototypes
 *  for biglobal stability of viscous flows.
 */

#include <petscksp.h>

/*
 *  Description:
 *
 *  grid2d is a structure representing the 2-d grid.
 *  nx is number of intervals in x direction. In other words x direction has
 *  nx+1 points.
 *  ny is number of intervals in y direction. In other words y direction has
 *  ny+1 points.
 *  qx is degree of polynomial to be used for x direction in case fd-q method
 *  is used for grid generation.
 *  qy is degree of polynomial to be used for y direction in case fd-q method
 *  is used for grid generation.
 *  x and y are petsc vectors containing coordinates of x and y grid.
 *  
 */

/*  Description:
 *  
 *  base_flow_ is the structure which stores base flow properties i.e. density,
 *  velocity and temperature.
 */

/*
 *  Description:
 *
 *  ndr_data is the biggest datastructure in terms of memory usage for biglobal 
 *  problem. It contains differentiation matrices and final dispersion relation
 *  matrices. Dispersion relation for biglobal stability problem is given as
 *                            
 *                               A * q = omega * B * q
 * 
 *  All the matrices are of Petsc Matrix type i.e. Mat.
 *  
 */


typedef struct co_trnfm_s {
    double *dxdX;
    double *dydX;
    double *dxdY;
    double *dydY;
    double *d2xdX2;
    double *d2ydX2;
    double *d2xdY2;
    double *d2ydY2;
    double *d2xdXY;
    double *d2ydXY;
} co_trnfm_t;

typedef struct ndr_data_s {
    PetscScalar re;
    PetscScalar pr;
    PetscScalar gamma;
    int nx;
    int ny;
    double *x;
    double *y;
    double *rho;
    double *drhodx;
    double *drhody;
    double *d2rhodx;
    double *d2rhody;
    double *u;
    double *dudx;
    double *dudy;
    double *v;
    double *dvdx;
    double *dvdy;
    double *w;
    double *dwdx;
    double *dwdy;
    co_trnfm_t trnfm;
    Mat A;
    Mat B;
    Vec q;
} ndr_data_t;

/*
 *  Description:
 *
 *  grid_gen is a function pointer meant to hold fd_grid or cheb_grid functions.
 *  It takes data of type struct stability_par by reference and generates the x
 *  and y vectors based on grid parameters nx, ny, qx, qy. fd_grid generates x
 *  and y vectors based on fd-q method discussed in hermann and hernandez. 
 *  cheb_grid generates x and y vectors which are gauss-lobato points. qx and qy
 *  have no meaning and are redundant parameters for cheb_grid. If everything is
 *  correct function returns 0 otherwise a non zero value.
 *  
 */
int cheb_grid(ndr_data_t *);

/*
 *  Description:
 *
 *  diffr_mat is a function pointer meant to hold fd_diffr_mat or cheb_diffr_mat
 *  functions. It takes data of type struct stability_par and struct ndr_data by
 *  reference and generates the first and second order differentiation matrices
 *  with respect to x and y. Cross derivatives are not needed hence not generated.
 *  prefix fd and cheb in __diffr_mat are functions named based on method for which
 *  they generate diffferentiation matrices. If everything is correct function returns
 *  0 otherwise a non zero value.
 *
 */



int cheb_diffr_mat(ndr_data_t *);

void get_trnfm_der(int, int, co_trnfm_t *); 
void trnfm_geom(ndr_data_t *);
void cheb_d1dx(int,int,int,PetscScalar *); 
void cheb_d1dy(int,int,int,PetscScalar *); 
void cheb_d2dx(int,int,int,PetscScalar *); 
void cheb_d2dy(int,int,int,PetscScalar *); 
void cheb_d2dxy(int,int,int,PetscScalar *); 

void d1dx(int,int,int,co_trnfm_t *,PetscScalar *); 
void d1dy(int,int,int,co_trnfm_t *,PetscScalar *); 
void d2dx(int,int,int,co_trnfm_t *,PetscScalar *); 
void d2dy(int,int,int,co_trnfm_t *,PetscScalar *); 
void d2dxy(int,int,int,co_trnfm_t *,PetscScalar *); 
/*
 *  Note that separate functions are need for fd method and chebyshev method
 *  for grid generation, differentiation matrix creation and boundary condition
 *  implementation. Numerical dispersion relation however can be same as every
 *  calculation is done by Petsc matrix vector operations.
 */

int helmholtz(PetscScalar *, PetscScalar *, ndr_data_t *);
int lns_biglobal(int, PetscScalar *, PetscScalar *, ndr_data_t *);
int set_bc(ndr_data_t *);
int eigen_solver(ndr_data_t *);

void *vec_alloc(int, int, size_t);
void free_vec(void *, int, size_t);
void **mat2_alloc(int,int,int,size_t);
void free_mat2(void **,int,size_t);

#endif // _BIGLOBAL_DATA_H_
