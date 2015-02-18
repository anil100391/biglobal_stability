// This program implements biglobal stability analysis.
static char help[] = "You are doomed";
#include <stdio.h>
#include <math.h>
#include <slepceps.h>
#include <data.h>

int read_base_flow(ndr_data_t *);
#undef __FUNCT__
#define __FUNCTI__ "main"
int main(int argc, char **args)
{
//    int i,j;
    int z;
//    PetscScalar b,o;
    ndr_data_t ndrdata;
    PetscErrorCode ierr;
    PetscViewer plotmat;
//    PetscViewer writemat;
    PetscScalar beta, omega;
    PetscMPIInt rank,size;
//    double *vec;
//    PetscScalar *der;
//    double sum = 0;
//    double x,y;
//    int nx, ny;
    ndrdata.nx = 10;
    ndrdata.ny = 10;
    ndrdata.re = 1000;
    ndrdata.pr = 1.0;
    ndrdata.gamma = 1.4;
//    nx = ndrdata.nx;
//    ny = ndrdata.ny;
    
    beta = 1.0*M_PI;
//    o = 0;
    
/*  how to print complex numbers
    printf("%lf + %lfi\n", creal(data.alpha), cimag(data.alpha));
*/    
    SlepcInitialize(&argc,&args,(char*)0,help);
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);
    printf("Size of petsc world is %d\n",size);

    cheb_grid(&ndrdata);
    read_base_flow(&ndrdata);
    lns_biglobal(&beta,&omega,&ndrdata);
//    helmholtz(&beta,&omega,&ndrdata);
    //eigen_solver(&ndrdata); 
/*    vec = (double *)vec_alloc((ndrdata.nx+1)*(ndrdata.ny+1),0,sizeof(double));
    der = (PetscScalar *)vec_alloc((ndrdata.nx+1)*(ndrdata.ny+1),0,sizeof(PetscScalar));
    for (i=0; i<(nx+1); i++) {
        for (j=0; j<(ny+1); j++) {
            x = ndrdata.x[i]; y = ndrdata.y[j];
            vec[i*(ny+1)+j] = exp(x*x+y*y);
        }
    }
    z = 0;
    for (z=0; z<(nx+1)*(ny+1); z++) {
    d2dxy(z,ndrdata.nx,ndrdata.ny,der);
   // for (i=0; i<(nx+1)*(ny+1); i++) printf("%lf\n",der[i]);
    sum = 0;
    for (i=0; i<(nx+1)*(ny+1); i++) {
        sum = sum + der[i]*vec[i];
    }
    x = ndrdata.x[z/(ny+1)]; y = ndrdata.y[z%(ny+1)];
    printf("dq_theo - dq_num = %lf \n",4*x*y*exp(x*x+y*y)-sum);
   }
*/
     
    ierr = PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,"A",0,0,800,800,&plotmat);
    MatView(ndrdata.A,plotmat);
//    ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"A",FILE_MODE_WRITE,&writemat);
//    MatView(data.A,writemat);
//    ierr = PetscViewerDestroy(&writemat); CHKERRQ(ierr);
    eigen_solver(&ndrdata);
//    scanf("%d",&z); 
    printf("All done\n");
    ierr = SlepcFinalize(); CHKERRQ(ierr);
    return 0;
}




