/* -------------------------------------------------------------------------
 *   This file contains function definitions of following functions
 *
 *   1. cheb_grid - takes ndr_data_t as input by reference
 *                - creates the x and y grid based on gauss lobatto points
 *                - writes x and y grid to ndr_data_t.x and ndr_data_t.y
 *
 *   2. cheb_diffr_mat - takes ndr_data_t as input by reference
 *                     - creates differentiation matrices Dx, Dy, D2x, D2y
 *                       and Dxy and writes them to corresponding elements
 *                       of ndr_data_t structure.
 * -------------------------------------------------------------------------
 */

#include <data.h>
#include <math.h>

/*****************************************************************************/
int cheb_grid(ndr_data_t *arg)
{
    int i;
    int nx=arg->nx, ny=arg->ny;

    //------- Create Vectors x and y -------//
    arg->x = (double *)vec_alloc(nx+1,0,sizeof(double));
    arg->y = (double *)vec_alloc(ny+1,0,sizeof(double));

    //------- Assign Values to Vectors -------//   
    for (i=0; i<ny+1; i++) {
        arg->y[i] = cos((double)i*M_PI/ny);
    }

    for (i=0; i<nx+1; i++) {
        arg->x[i] = cos((double)i*M_PI/nx);
        //printf("%lf\n",arg->
    }
    
    return 0;
}

/*****************************************************************************/


void cheb_d1dy(int row,int nx,int ny,PetscScalar *diff) 
{
    int i,j,locrow,loccol;
    int ci, cj;
    
    locrow = row%(ny+1);
    loccol = row/(ny+1);
    
    for (i=0; i<nx+1; i++) {
        for (j=0; j<ny+1; j++) {
            diff[i*(ny+1)+j] = 0.0;
        }
    }
  
    for (j=0; j<ny+1; j++) {
        if (j==locrow) {
            if (locrow==0) diff[loccol*(ny+1)+j] = (2*ny*ny+1)/6.0;
            else if (locrow==ny) diff[loccol*(ny+1)+j] = -(2*ny*ny+1)/6.0;
            else diff[loccol*(ny+1)+j] = -cos(j*M_PI/ny) / (2*pow(sin(j*M_PI/ny),2));
        }
        else {
            ci = (locrow==0 || locrow==ny) ? 2:1;
            cj = (j==0 || j==ny) ? 2:1;
            i = locrow;
            diff[loccol*(ny+1)+j] = -ci*pow(-1,i+j) / (2*cj*sin((i-j)*M_PI_2/ny)*sin((i+j)*M_PI_2/ny));
        }
    }

}

void cheb_d1dx(int row,int nx,int ny,PetscScalar *diff) 
{
    int i,j,locrow,loccol;
    int ci, cj;
    double *dx;
    
    loccol = row%(ny+1);
    locrow = row/(ny+1);
    
    dx   = (double *)vec_alloc(nx+1,0,sizeof(double));
    
    for (i=0; i<nx+1; i++) {
        for (j=0; j<ny+1; j++) {
            diff[i*(ny+1)+j] = 0.0;
        }
    }
    
    for (j=0; j<nx+1; j++) {
        if (locrow==j) {
            if (locrow==0) dx[j] = (2*nx*nx+1)/6.0;
            else if (locrow==nx) dx[j] = -(2*nx*nx+1)/6.0;
            else dx[j] = -cos(j*M_PI/nx) / (2*pow(sin(j*M_PI/nx),2));
        }
        else {
            ci = (locrow==0 || locrow==nx) ? 2:1;
            cj = (j==0 || j==nx) ? 2:1;
            i = locrow;
            dx[j] = -ci*pow(-1,i+j) / (2*cj*sin((i-j)*M_PI_2/nx)*sin((i+j)*M_PI_2/nx));
        }
    }
    i = 0;
    for (j=loccol; j<(nx+1)*(ny+1); j=j+(ny+1)) {
        diff[j]=dx[i++];
    }
    
    free_vec((void *)dx,0,sizeof(double));

}

void cheb_d2dy(int row,int nx,int ny,PetscScalar *diff) 
{
    int i,j,locrow,loccol;
    int cj;
    double xi,xj;
    
    locrow = row%(ny+1);
    loccol = row/(ny+1);
    
    for (i=0; i<nx+1; i++) {
        for (j=0; j<ny+1; j++) {
            diff[i*(ny+1)+j] = 0.0;
        }
    }

    xi = cos(locrow*M_PI/ny);
    
    for (j=0; j<ny+1; j++) {
        xj = cos(j*M_PI/ny);
        if (locrow<ny && locrow>0) {
            if (j!=locrow) {
                cj = (j==0 || j==ny) ? 2:1;
                diff[loccol*(ny+1)+j] = pow(-1.0,locrow+j)*(xi*xi+xi*xj-2)/(cj*(1-xi*xi)*(xi-xj)*(xi-xj));
            }
            else {
                diff[loccol*(ny+1)+j] = -((ny*ny-1)*(1-xi*xi)+3.0)/(3.0*(1-xi*xi)*(1-xi*xi));
            }
        }
        else if (locrow==0) {
            if (j>0 && j<=ny) {
                cj = (j==0 || j==ny) ? 2:1;
                diff[loccol*(ny+1)+j] = 2.0*pow(-1.0,j)*((2*ny*ny+1)*(1-xj)-6)/(3.0*cj*(1-xj)*(1-xj));
            }
            else if (j==0) {
                diff[loccol*(ny+1)+j] = (pow(ny,4) - 1)/15.0;
            }
        }
        else if (locrow==ny) {
            if (j>=0 && j<ny) {
                cj = (j==0 || j==ny) ? 2:1;
                diff[loccol*(ny+1)+j] = 2.0*pow(-1.0,j+ny)*((2*ny*ny+1)*(1+xj)-6)/(3.0*cj*(1+xj)*(1+xj));
            }
            else if (j==ny) {
                diff[loccol*(ny+1)+j] = (pow(ny,4) - 1)/15.0;
            }
        }
    }
            

}

void cheb_d2dx(int row,int nx,int ny,PetscScalar *diff) 
{
    int i,j,locrow,loccol;
    int cj;
    double *d2x;
    double xi,xj;
    
    loccol = row%(ny+1);
    locrow = row/(ny+1);
    
    d2x   = (double *)vec_alloc(nx+1,0,sizeof(double));
    
    for (i=0; i<nx+1; i++) {
        for (j=0; j<ny+1; j++) {
            diff[i*(ny+1)+j] = 0.0;
        }
    }
    
    xi = cos(locrow*M_PI/nx);
    
    for (j=0; j<nx+1; j++) {
        xj = cos(j*M_PI/nx);
        if (locrow<nx && locrow>0) {
            if (j!=locrow) {
                cj = (j==0 || j==nx) ? 2:1;
                d2x[j] = pow(-1.0,locrow+j)*(xi*xi+xi*xj-2)/(cj*(1-xi*xi)*(xi-xj)*(xi-xj));
            }
            else {
                d2x[j] = -((nx*nx-1)*(1-xi*xi)+3)/(3.0*(1-xi*xi)*(1-xi*xi));
            }
        }
        else if (locrow==0) {
            if (j>0 && j<=nx) {
                cj = (j==0 || j==nx) ? 2:1;
                d2x[j] = 2.0*pow(-1.0,j)*((2*nx*nx+1)*(1-xj)-6)/(3.0*cj*(1-xj)*(1-xj));
            }
            else if (j==0) {
                d2x[j] = (pow(nx,4) - 1)/15.0;
            }
        }
        else if (locrow==nx) {
            if (j>=0 && j<nx) {
                cj = (j==0 || j==nx) ? 2:1;
                d2x[j] = 2.0*pow(-1.0,j+nx)*((2*nx*nx+1)*(1+xj)-6)/(3.0*cj*(1+xj)*(1+xj));
            }
            else if (j==nx) {
                d2x[j] = (pow(nx,4) - 1)/15.0;
            }
        }
    }
    i=0;
    for (j=loccol; j<(nx+1)*(ny+1); j=j+(ny+1)) {
        diff[j]=d2x[i++];
    }
    
    free_vec((void *)d2x,0,sizeof(double));

}

void cheb_d2dxy(int row,int nx,int ny,PetscScalar *diff) 
{
    int i,j,locrow;
    int ci, cj;
    double *dx,**dy;
    
    locrow = row/(ny+1);
    
    dx   = (double *)vec_alloc(nx+1,0,sizeof(double));
    dy   = (double **)mat2_alloc(ny+1,ny+1,0,sizeof(double));
    
    for (i=0; i<nx+1; i++) {
        for (j=0; j<ny+1; j++) {
            diff[i*(ny+1)+j] = 0.0;
        }
    }
    
    ///////////////// Create dx ////////////////
    for (j=0; j<nx+1; j++) {
        if (locrow==j) {
            if (locrow==0) dx[j] = (2*nx*nx+1)/6.0;
            else if (locrow==nx) dx[j] = -(2*nx*nx+1)/6.0;
            else dx[j] = -cos(j*M_PI/nx) / (2*pow(sin(j*M_PI/nx),2));
        }
        else {
            ci = (locrow==0 || locrow==nx) ? 2:1;
            cj = (j==0 || j==nx) ? 2:1;
            i = locrow;
            dx[j] = -ci*pow(-1,i+j) / (2*cj*sin((i-j)*M_PI_2/nx)*sin((i+j)*M_PI_2/nx));
        }
    }

    /////////////// Create dy /////////////////
    for (i=0; i<ny+1; i++) {    
        for (j=0; j<ny+1; j++) {
            if (j==i) {
                if (i==0) dy[i][j] = (2*ny*ny+1)/6.0;
                else if (i==ny) dy[i][j] = -(2*ny*ny+1)/6.0;
                else dy[i][j] = -cos(j*M_PI/ny) / (2*pow(sin(j*M_PI/ny),2));
            }
            else {
                ci = (i==0 || i==ny) ? 2:1;
                cj = (j==0 || j==ny) ? 2:1;
                dy[i][j] = -ci*pow(-1,i+j) / (2*cj*sin((i-j)*M_PI_2/ny)*sin((i+j)*M_PI_2/ny));
            }
        }
    }

    // reassign localrow
    locrow = row%(ny+1);
   
    for (i=0; i<(nx+1); i++) {
        for (j=0; j<(ny+1); j++) {
            diff[i*(ny+1)+j] = dx[i]*dy[locrow][j];
        }
    }
    
    free_vec((void *)dx,0,sizeof(double));
    free_mat2((void **)dy,0,sizeof(double));

}


void get_trnfm_der(int nx, int ny, co_trnfm_t *arg) 
{
    int i;
    arg->dxdX = (double *)vec_alloc((nx+1)*(ny+1),0,sizeof(double));
    arg->dydX = (double *)vec_alloc((nx+1)*(ny+1),0,sizeof(double));
    arg->dxdY = (double *)vec_alloc((nx+1)*(ny+1),0,sizeof(double));
    arg->dydY = (double *)vec_alloc((nx+1)*(ny+1),0,sizeof(double));
    arg->d2xdX2 = (double *)vec_alloc((nx+1)*(ny+1),0,sizeof(double));
    arg->d2ydX2 = (double *)vec_alloc((nx+1)*(ny+1),0,sizeof(double));
    arg->d2xdY2 = (double *)vec_alloc((nx+1)*(ny+1),0,sizeof(double));
    arg->d2ydY2 = (double *)vec_alloc((nx+1)*(ny+1),0,sizeof(double));
    arg->d2xdXY = (double *)vec_alloc((nx+1)*(ny+1),0,sizeof(double));
    arg->d2ydXY = (double *)vec_alloc((nx+1)*(ny+1),0,sizeof(double));
    
    for (i=0; i<(nx+1)*(ny+1); i++) {
        arg->dxdX[i] = 0.025;
        arg->dydX[i] = 0.0;
        arg->dxdY[i] = 0.0;
        arg->dydY[i] = 0.2;
        arg->d2xdX2[i] = 0.0;
        arg->d2ydX2[i] = 0.0;
        arg->d2xdY2[i] = 0.0;
        arg->d2ydY2[i] = 0.0;
        arg->d2xdXY[i] = 0.0;
        arg->d2ydXY[i] = 0.0;
    }

}

void trnfm_geom(ndr_data_t *arg)
{
    int nx = arg->nx, ny=arg->ny;
    int i,j;

    for (i=0; i<nx+1; i++) {
        arg->x[i] = 40.0*(arg->x[i] + 1);
    }
    for (j=0; j<ny+1; j++) {
        arg->y[j] = 5.0*(arg->y[j] + 1);
    }
}

void d1dx(int row,int nx,int ny,co_trnfm_t *geom, PetscScalar *diff)
{
    int i,locrow;
    PetscScalar *derx, *dery;
    double xX, yX;

    locrow = row%(ny+1);

    derx = (PetscScalar *)vec_alloc((nx+1)*(ny+1),0,sizeof(PetscScalar));
    dery = (PetscScalar *)vec_alloc((nx+1)*(ny+1),0,sizeof(PetscScalar));
    
    cheb_d1dx(row,nx,ny,derx);
    cheb_d1dy(row,nx,ny,dery);

    xX = (*geom).dxdX[locrow];
    yX = (*geom).dydX[locrow];
    for (i=0; i<(nx+1)*(ny+1); i++) {
        diff[i] = xX*derx[i] + yX*dery[i];
    }

    free_vec((void *)derx,0,sizeof(double));
    free_vec((void *)dery,0,sizeof(double));    
}

 
void d1dy(int row,int nx,int ny,co_trnfm_t *geom, PetscScalar *diff)
{
    int i,locrow;
    PetscScalar *derx, *dery;
    double xY, yY;

    locrow = row%(ny+1);

    derx = (PetscScalar *)vec_alloc((nx+1)*(ny+1),0,sizeof(PetscScalar));
    dery = (PetscScalar *)vec_alloc((nx+1)*(ny+1),0,sizeof(PetscScalar));
    cheb_d1dx(row,nx,ny,derx);
    cheb_d1dy(row,nx,ny,dery);

    xY = (*geom).dxdY[locrow];
    yY = (*geom).dydY[locrow];
    for (i=0; i<(nx+1)*(ny+1); i++) {
        diff[i] = xY*derx[i] + yY*dery[i];
    }

    free_vec((void *)derx,0,sizeof(double));
    free_vec((void *)dery,0,sizeof(double));    
}

 
void d2dx(int row,int nx,int ny,co_trnfm_t *geom, PetscScalar *diff)
{
    int i,locrow;
    PetscScalar *derx1, *derx2, *dery1, *dery2, *derxy;
    double xX, yX, xX2, yX2;

    locrow = row%(ny+1);

    derx1 = (PetscScalar *)vec_alloc((nx+1)*(ny+1),0,sizeof(PetscScalar));
    dery1 = (PetscScalar *)vec_alloc((nx+1)*(ny+1),0,sizeof(PetscScalar));
    derx2 = (PetscScalar *)vec_alloc((nx+1)*(ny+1),0,sizeof(PetscScalar));
    dery2 = (PetscScalar *)vec_alloc((nx+1)*(ny+1),0,sizeof(PetscScalar));
    derxy = (PetscScalar *)vec_alloc((nx+1)*(ny+1),0,sizeof(PetscScalar));

    cheb_d1dx(row,nx,ny,derx1);
    cheb_d1dy(row,nx,ny,dery1);
    cheb_d2dx(row,nx,ny,derx2);
    cheb_d2dy(row,nx,ny,dery2);
    cheb_d2dxy(row,nx,ny,derxy);

    xX = (*geom).dxdX[locrow];
    yX = (*geom).dydX[locrow];
    xX2 = (*geom).d2xdX2[locrow];
    yX2 = (*geom).d2ydX2[locrow];

    for (i=0; i<(nx+1)*(ny+1); i++) {
        diff[i] = xX*xX*derx2[i] + 2*xX*yX*derxy[i] + yX*yX*dery2[i] + xX2*derx1[i] + yX2*dery1[i];
    }

    free_vec((void *)derx1,0,sizeof(double));
    free_vec((void *)dery1,0,sizeof(double));    
    free_vec((void *)derx2,0,sizeof(double));
    free_vec((void *)dery2,0,sizeof(double));    
    free_vec((void *)derxy,0,sizeof(double));

} 


void d2dy(int row,int nx,int ny,co_trnfm_t *geom, PetscScalar *diff)
{
    int i,locrow;
    PetscScalar *derx1, *derx2, *dery1, *dery2, *derxy;
    double xY, yY, xY2, yY2;

    locrow = row%(ny+1);

    derx1 = (PetscScalar *)vec_alloc((nx+1)*(ny+1),0,sizeof(PetscScalar));
    dery1 = (PetscScalar *)vec_alloc((nx+1)*(ny+1),0,sizeof(PetscScalar));
    derx2 = (PetscScalar *)vec_alloc((nx+1)*(ny+1),0,sizeof(PetscScalar));
    dery2 = (PetscScalar *)vec_alloc((nx+1)*(ny+1),0,sizeof(PetscScalar));
    derxy = (PetscScalar *)vec_alloc((nx+1)*(ny+1),0,sizeof(PetscScalar));

    cheb_d1dx(row,nx,ny,derx1);
    cheb_d1dy(row,nx,ny,dery1);
    cheb_d2dx(row,nx,ny,derx2);
    cheb_d2dy(row,nx,ny,dery2);
    cheb_d2dxy(row,nx,ny,derxy);

    xY = (*geom).dxdY[locrow];
    yY = (*geom).dydY[locrow];
    xY2 = (*geom).d2xdY2[locrow];
    yY2 = (*geom).d2ydY2[locrow];

    for (i=0; i<(nx+1)*(ny+1); i++) {
        diff[i] = xY*xY*derx2[i] + 2*xY*yY*derxy[i] + yY*yY*dery2[i] + xY2*derx1[i] + yY2*dery1[i];
    }

    free_vec((void *)derx1,0,sizeof(double));
    free_vec((void *)dery1,0,sizeof(double));    
    free_vec((void *)derx2,0,sizeof(double));
    free_vec((void *)dery2,0,sizeof(double));    
    free_vec((void *)derxy,0,sizeof(double));
} 


void d2dxy(int row,int nx,int ny,co_trnfm_t *geom, PetscScalar *diff)
{
    int i,locrow;
    PetscScalar *derx1, *derx2, *dery1, *dery2, *derxy;
    double xX, xY, yX, yY, xXY, yXY;

    locrow = row%(ny+1);

    derx1 = (PetscScalar *)vec_alloc((nx+1)*(ny+1),0,sizeof(PetscScalar));
    dery1 = (PetscScalar *)vec_alloc((nx+1)*(ny+1),0,sizeof(PetscScalar));
    derx2 = (PetscScalar *)vec_alloc((nx+1)*(ny+1),0,sizeof(PetscScalar));
    dery2 = (PetscScalar *)vec_alloc((nx+1)*(ny+1),0,sizeof(PetscScalar));
    derxy = (PetscScalar *)vec_alloc((nx+1)*(ny+1),0,sizeof(PetscScalar));

    cheb_d1dx(row,nx,ny,derx1);
    cheb_d1dy(row,nx,ny,dery1);
    cheb_d2dx(row,nx,ny,derx2);
    cheb_d2dy(row,nx,ny,dery2);
    cheb_d2dxy(row,nx,ny,derxy);

    xX = (*geom).dxdX[locrow];
    xY = (*geom).dxdY[locrow];
    yX = (*geom).dydX[locrow];
    yY = (*geom).dydY[locrow];
    xXY = (*geom).d2xdXY[locrow];
    yXY = (*geom).d2ydXY[locrow];

    for (i=0; i<(nx+1)*(ny+1); i++) {
        diff[i] = xX*xY*derx2[i] + (xX*yY+xY*yX)*derxy[i] + yX*yY*dery2[i] + xXY*derx1[i] + yXY*dery1[i];
    }

    free_vec((void *)derx1,0,sizeof(double));
    free_vec((void *)dery1,0,sizeof(double));    
    free_vec((void *)derx2,0,sizeof(double));
    free_vec((void *)dery2,0,sizeof(double));    
    free_vec((void *)derxy,0,sizeof(double));
} 
