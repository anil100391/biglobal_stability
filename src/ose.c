#include <stdio.h>
#include <data.h>
#include <slepceps.h>
#include <math.h>

int orr_sommerfeld(PetscScalar *alpha, PetscScalar *omega, ndr_data_t *arg)
{
    int i,j,Ii,Istart,Iend;
    int ny = arg->ny;
    int M = ny+1;
    double *y;
    double **T, **dT, **d2T, **d3T, **d4T;
    int *col;
    PetscScalar *val,L1,L2,asqplbsq,beta=0.0;
    PetscErrorCode ierr;

    y = (double *)vec_alloc(ny+1,0,sizeof(double));
    for (j=0; j<ny+1; j++)
        y[j] = cos(j*M_PI/ny);

    T = (double **)mat2_alloc(ny+1,ny+1,0,sizeof(double));
    dT = (double **)mat2_alloc(ny+1,ny+1,0,sizeof(double));
    d2T = (double **)mat2_alloc(ny+1,ny+1,0,sizeof(double));
    d3T = (double **)mat2_alloc(ny+1,ny+1,0,sizeof(double));
    d4T = (double **)mat2_alloc(ny+1,ny+1,0,sizeof(double));

    for (i=0; i<ny+1; i++) {
        for (j=0; j<ny+1; j++) {
            T[i][j] = 0;
            dT[i][j] = 0;
            d2T[i][j] = 0;
            d3T[i][j] = 0;
            d4T[i][j] = 0;
        }
    }

    for (i=2; i<=ny-2; i++) {
        T[i][0] = 1.0;
        T[i][1] = y[i];
    }
    for (i=2; i<=ny-2; i++)
        for (j=2; j<ny+1; j++)
            T[i][j] = 2.0*y[i]*T[i][j-1]-T[i][j-2];
    for (i=2; i<=ny-2; i++) {
        dT[i][0] = 0.0;
        dT[i][1] = 1.0;
        dT[i][2] = 4.0*y[i];
    }    
    for (i=2; i<=ny-2; i++)
        for (j=3; j<ny+1; j++)
            dT[i][j] = 2.0*j*T[i][j-1]+(1.0*j/(j-2))*dT[i][j-2];
        
    for (i=2; i<=ny-2; i++) {
        d2T[i][0] = 0.0;
        d2T[i][1] = 0.0;
        d2T[i][2] = 4.0;
    }    
    for (i=2; i<=ny-2; i++)
        for (j=3; j<ny+1; j++)
            d2T[i][j] = 2.0*j*dT[i][j-1]+(1.0*j/(j-2))*d2T[i][j-2];

    for (i=2; i<=ny-2; i++) {
        d3T[i][0] = 0.0;
        d3T[i][1] = 0.0;
        d3T[i][2] = 0.0;
    }    
    for (i=2; i<=ny-2; i++)
        for (j=3; j<ny+1; j++)
            d3T[i][j] = 2.0*j*d2T[i][j-1]+(1.0*j/(j-2))*d3T[i][j-2];

    for (i=2; i<=ny-2; i++) {
        d4T[i][0] = 0.0;
        d4T[i][1] = 0.0;
        d4T[i][2] = 0.0;
    }    
    for (i=2; i<=ny-2; i++)
        for (j=3; j<ny+1; j++)
            d4T[i][j] = 2.0*j*d3T[i][j-1]+(1.0*j/(j-2))*d4T[i][j-2];

    col = (int *)vec_alloc(ny+1,0,sizeof(int));
    val = (PetscScalar *)vec_alloc(ny+1,0,sizeof(PetscScalar));
    for (j=0; j<ny+1; j++) 
        col[j] = j;

    ierr = MatCreate(PETSC_COMM_WORLD,&(arg->B)); CHKERRQ(ierr);
    ierr = MatSetSizes(arg->B,PETSC_DECIDE,PETSC_DECIDE,M,M); CHKERRQ(ierr);
    ierr = MatSetFromOptions(arg->B); CHKERRQ(ierr);
    ierr = MatSetUp(arg->B); CHKERRQ(ierr);
    
    ierr = MatCreate(PETSC_COMM_WORLD,&(arg->A)); CHKERRQ(ierr);
    ierr = MatSetSizes(arg->A,PETSC_DECIDE,PETSC_DECIDE,M,M); CHKERRQ(ierr);
    ierr = MatSetFromOptions(arg->A); CHKERRQ(ierr);
    ierr = MatSetUp(arg->A); CHKERRQ(ierr);
    
    ierr = MatGetOwnershipRange(arg->B,&Istart,&Iend); CHKERRQ(ierr);
    
    for (Ii=Istart; Ii<Iend; Ii++) {
        // boundary conditions //
        if (Ii==0) {
             for (j=0; j<ny+1; j++) val[j] = 1.0;
             ierr = MatSetValues(arg->B,1,&Ii,ny+1,col,val,INSERT_VALUES); CHKERRQ(ierr);
        }
        else if (Ii==1) {
             for (j=0; j<ny+1; j++) val[j] = j*j;
             ierr = MatSetValues(arg->B,1,&Ii,ny+1,col,val,INSERT_VALUES); CHKERRQ(ierr);
        }
        else if (Ii==ny-1) {
             for (j=0; j<ny+1; j++) val[j] = j*j*pow(-1,j-1);
             ierr = MatSetValues(arg->B,1,&Ii,ny+1,col,val,INSERT_VALUES); CHKERRQ(ierr);
        }
        else if (Ii==ny) {
             for (j=0; j<ny+1; j++) val[j] = pow(-1,j);
             ierr = MatSetValues(arg->B,1,&Ii,ny+1,col,val,INSERT_VALUES); CHKERRQ(ierr);
        }
        // other rows filling //
        else {
             for (j=0; j<ny+1; j++) {
                 val[j] = d2T[Ii][j] - ((*alpha)*(*alpha) + beta*beta)*T[Ii][j];
             }
             ierr = MatSetValues(arg->B,1,&Ii,ny+1,col,val,INSERT_VALUES); CHKERRQ(ierr);
        }
      

    }


    ierr = MatGetOwnershipRange(arg->A,&Istart,&Iend); CHKERRQ(ierr);
    asqplbsq = ((*alpha)*(*alpha)+beta*beta);
    for (Ii=Istart; Ii<Iend; Ii++) {
        if (Ii>1 && Ii<ny-1) {
            L1 = -(1-y[Ii]*y[Ii])*asqplbsq + 2 + asqplbsq*asqplbsq*I/(*alpha*arg->re);
            L2 = (1-y[Ii]*y[Ii]) - 2.0*asqplbsq*I/(*alpha*arg->re);
            for (j=0; j<ny+1; j++) {
                val[j] = L1*T[Ii][j] + L2*d2T[Ii][j] + d4T[Ii][j]*I/(*alpha*arg->re);
            }
            ierr = MatSetValues(arg->A,1,&Ii,ny+1,col,val,INSERT_VALUES); CHKERRQ(ierr);
        }
    }


    ierr = MatAssemblyBegin(arg->A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(arg->A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyBegin(arg->B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(arg->B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    return 0;
}
