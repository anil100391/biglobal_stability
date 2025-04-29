#ifndef _BIGLOBAL_DATA_H_
#define _BIGLOBAL_DATA_H_

// -----------------------------------------------------------------------------
// This file contains definition of data structures and function prototypes
// for biglobal stability of viscous flows.
// -----------------------------------------------------------------------------
// nx is number of intervals in x direction. In other words x direction has
// nx + 1 points.
// ny is number of intervals in x direction. In other words x direction has
// ny + 1 points.
// x and y are petsc vectors to store coordinates of x-y grid
// ndr_data is the main datastructure for biglobal problem. It contains
// differentiation matrices and final dispersion relation matrices. Dispersion
// relation for biglobal stability problem is given as
//
//                              A * q = omega * B * q
//
// All the matrices are of Petsc Matrix type i.e. Mat.
// -----------------------------------------------------------------------------

#include <petscksp.h>

// -----------------------------------------------------------------------------
// coordinate transformation structure
// -----------------------------------------------------------------------------
typedef struct co_trnfm_s
{
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

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
typedef struct ndr_data_s
{
    // Reynolds number
    PetscScalar re;

    // Prandtl number
    PetscScalar pr;

    PetscScalar gamma;

    // number of intervals in x direction i.e. nx + 1 points in x direction
    int nx;

    // number of intervals in y direction i.e. ny + 1 points in y direction
    int ny;

    // x-y grid
    double *x;
    double *y;

    double *rho;
    double *drhodx;
    double *drhody;
    double *d2rhodx;
    double *d2rhody;

    // velocity and it's partial derivatives
    double *u;
    double *dudx;
    double *dudy;
    double *v;
    double *dvdx;
    double *dvdy;
    double *w;
    double *dwdx;
    double *dwdy;

    // coordinate transformation
    co_trnfm_t trnfm;

    // dispersion relation matrices
    Mat A;
    Mat B;
    Vec q;
} ndr_data_t;

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// cheb_grid takes data of type struct ndr_data_t by reference and generates
// the x and y vectors based on grid parameters nx, ny. cheb_grid generates x
// and y vectors which are gauss-lobato points. Returns 0 if success
int cheb_grid(ndr_data_t *);

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// cheb_diffr_mat takes data of type struct ndr_data_t by reference and
// generates the first and second order differentiation matrices with respect
// to x and y. Returns 0 if success
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

int helmholtz(PetscScalar *, PetscScalar *, ndr_data_t *);
int lns_biglobal(int, PetscScalar *, PetscScalar *, ndr_data_t *);
int set_bc(ndr_data_t *);
int eigen_solver(ndr_data_t *);

void *vec_alloc(int, int, size_t);
void free_vec(void *, int, size_t);
void **mat2_alloc(int,int,int,size_t);
void free_mat2(void **,int,size_t);

#endif // _BIGLOBAL_DATA_H_
