// This program implements biglobal stability analysis.
static char help[] = "Persist! results will come right";
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
    PetscViewer writemat;
    PetscScalar beta, omega;
    PetscMPIInt rank,size;
    PetscScalar *test;
//    double *vec;
//    PetscScalar *der;
//    double sum = 0;
//    double x,y;
//    int nx, ny;
    ndrdata.nx = 40;
    ndrdata.ny = 40;
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
    //printf("Size of petsc world is %d\n",size);

    cheb_grid(&ndrdata);
    trnfm_geom(&ndrdata);
    get_trnfm_der(ndrdata.nx,ndrdata.ny,&(ndrdata.trnfm));
    read_base_flow(&ndrdata);



    //test = (PetscScalar *)vec_alloc((ndrdata.nx+1)*(ndrdata.ny+1),0,sizeof(PetscScalar));
    //printf("%p\n",&(ndrdata.trnfm));
    //d1dx(1,ndrdata.nx,ndrdata.ny,&(ndrdata.trnfm),test);


    lns_biglobal(1,&beta,&omega,&ndrdata);
    set_bc(&ndrdata);
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
     
//    ierr = PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,"A",0,0,800,800,&plotmat);
    //MatView(ndrdata.A,plotmat);
/*    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"A",FILE_MODE_WRITE,&writemat);
    MatView(ndrdata.A,writemat);
    ierr = PetscViewerDestroy(&writemat); CHKERRQ(ierr);
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"B",FILE_MODE_WRITE,&writemat);
    MatView(ndrdata.B,writemat);
    ierr = PetscViewerDestroy(&writemat); CHKERRQ(ierr);
*/
    eigen_solver(&ndrdata);
//    scanf("%d",&z); 
    printf("All done\n");
    ierr = SlepcFinalize(); CHKERRQ(ierr);
    return 0;
}

int read_base_flow(ndr_data_t *arg)
{

    int nx=arg->nx, ny=arg->ny;
    int i;
    FILE *f;
    
    /*********************** Read Base Flow from files **********************/
    arg->w = (double *)vec_alloc((nx+1)*(ny+1),0,sizeof(double));
    f = fopen("w.txt","r");
    if (f==NULL) {
        printf("Error while opening the file.\n");
    }
    i=0;
    while (fscanf(f,"%lf \n",&(arg->w[i])) !=EOF) {
        i = i+1;
    }

    fclose(f);

    arg->dwdx = (double *)vec_alloc((nx+1)*(ny+1),0,sizeof(double));
    f = fopen("dwdx.txt","r");
    if (f==NULL) {
        printf("Error while opening the file.\n");
    }
    i=0;
    while (fscanf(f,"%lf \n",&(arg->dwdx[i])) !=EOF) {
        i = i+1;
    }
   
    fclose(f);


    arg->dwdy = (double *)vec_alloc((nx+1)*(ny+1),0,sizeof(double));
    f = fopen("dwdy.txt","r");
    if (f==NULL) {
        printf("Error while opening the file.\n");
    }
    i=0;
    while (fscanf(f,"%lf \n",&(arg->dwdy[i])) !=EOF) {
        i = i+1;
    }
    
    fclose(f);

    // Set other known quantities for rectangular channel flow
     

    /*************************** Baseflow matrices set **************************/
    
    arg->rho = (double *)vec_alloc((nx+1)*(ny+1),0,sizeof(double));
    arg->drhodx = (double *)vec_alloc((nx+1)*(ny+1),0,sizeof(double));
    arg->drhody = (double *)vec_alloc((nx+1)*(ny+1),0,sizeof(double));
    arg->d2rhodx = (double *)vec_alloc((nx+1)*(ny+1),0,sizeof(double));
    arg->d2rhody = (double *)vec_alloc((nx+1)*(ny+1),0,sizeof(double));
    arg->u = (double *)vec_alloc((nx+1)*(ny+1),0,sizeof(double));
    arg->dudx = (double *)vec_alloc((nx+1)*(ny+1),0,sizeof(double));
    arg->dudy = (double *)vec_alloc((nx+1)*(ny+1),0,sizeof(double));
    arg->v = (double *)vec_alloc((nx+1)*(ny+1),0,sizeof(double));
    arg->dvdx = (double *)vec_alloc((nx+1)*(ny+1),0,sizeof(double));
    arg->dvdy = (double *)vec_alloc((nx+1)*(ny+1),0,sizeof(double));

    for (i=0; i<(nx+1)*(ny+1); i++) {
        arg->u[i] = 0; arg->dudx[i] = 0; arg->dudy[i] = 0;
        arg->v[i] = 0; arg->dvdx[i] = 0; arg->dvdy[i] = 0;
        arg->rho[i] = 1.0; arg->drhodx[i] = 0; arg->drhody[i] = 0;
        arg->d2rhodx[i] = 0; arg->d2rhody[i] = 0;
    }
    return 0;
}



