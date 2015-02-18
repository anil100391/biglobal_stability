// This program implements biglobal stability analysis for coflow wakes//
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
    int z;
    ndr_data_t ndrdata;
    PetscErrorCode ierr;
    PetscViewer plotmat;
    PetscViewer writemat;
    PetscScalar beta, omega;
    PetscMPIInt rank,size;
    
  EPS            eps;             /* eigenproblem solver context */
  ST             st;
  EPSType        type;
  PetscReal      tol;
  Vec            xr,xi,*Iv,*Cv;
  PetscInt       nev,maxit,i,its,lits,nconv,nini=0,ncon=0;
  char           filename[PETSC_MAX_PATH_LEN];
  PetscViewer    viewer;
  PetscBool      flg,evecs,ishermitian;
  PetscMPIInt    Istart,Iend;
    
    
    ndrdata.nx = 70;
    ndrdata.ny = 30;
    ndrdata.re = 100;
    ndrdata.pr = 1.0;
    ndrdata.gamma = 1.4;
    beta = 0;
    
    SlepcInitialize(&argc,&args,(char*)0,help);
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);
    printf("Size of petsc world is %d\n",size);

    cheb_grid(&ndrdata);
    trnfm_geom(&ndrdata);
    get_trnfm_der(ndrdata.nx,ndrdata.ny,&(ndrdata.trnfm));
    read_base_flow(&ndrdata);
    lns_biglobal(1,&beta,&omega,&ndrdata);
    printf("setting boundary conditions\n");
    set_bc(&ndrdata);
    ierr = PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,"A",0,0,800,800,&plotmat);
    //MatView(ndrdata.A,plotmat);
/*    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"A",FILE_MODE_WRITE,&writemat);
    MatView(ndrdata.A,writemat);
    ierr = PetscViewerDestroy(&writemat); CHKERRQ(ierr);
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"B",FILE_MODE_WRITE,&writemat);
    MatView(ndrdata.B,writemat);
    ierr = PetscViewerDestroy(&writemat); CHKERRQ(ierr);
*/
    


  

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Load the matrices that define the eigensystem, Ax=kBx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nGeneralized eigenproblem.\n\n");CHKERRQ(ierr);
  ierr = MatGetVecs(ndrdata.A,NULL,&xr);CHKERRQ(ierr);
  ierr = MatGetVecs(ndrdata.A,NULL,&xi);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);
  ierr = MatGetOwnershipRange(ndrdata.A,&Istart,&Iend);
  printf("process %d is responsible for rows %d to %d \n",rank,Istart,Iend);
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create eigensolver context
  */
  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
  //ierr = EPSGetST(eps, &st);
  //ierr = STSetType(st,STSINVERT);
  //ierr = STSetShift(st,4.0+4.0*I);
  /*
     Set operators. In this case, it is a generalized eigenvalue problem
  */
  ierr = EPSSetOperators(eps,ndrdata.A,ndrdata.B);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps, EPS_GNHEP); CHKERRQ(ierr);
  //ierr = EPSSetTolerances(EPS eps,PetscReal tol,PetscInt maxits); CHKERRQ(ierr);

  /*
     Set solver parameters at runtime
  */
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSSolve(eps);CHKERRQ(ierr);

  /*
     Optional: Get some information from the solver and display it
  */
  ierr = EPSGetIterationNumber(eps,&its);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %D\n",its);CHKERRQ(ierr);
  ierr = EPSGetST(eps,&st);CHKERRQ(ierr);
  ierr = STGetOperationCounters(st,NULL,&lits);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of linear iterations of the method: %D\n",lits);CHKERRQ(ierr);
  ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
  ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);CHKERRQ(ierr);
  ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%D\n",(double)tol,maxit);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSPrintSolution(eps,NULL);CHKERRQ(ierr);
  /*
     Save eigenvectors, if requested
  */
  ierr = PetscOptionsGetString(NULL,"-evecs",filename,PETSC_MAX_PATH_LEN,&evecs);CHKERRQ(ierr);
  ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
  if (nconv>0 && evecs) {
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer);
    //ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
    ierr = EPSIsHermitian(eps,&ishermitian);CHKERRQ(ierr);
    for (i=0;i<nconv;i++) {
      ierr = EPSGetEigenvector(eps,i,xr,xi);CHKERRQ(ierr);
      ierr = VecView(xr,viewer);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
      if (!ishermitian) { ierr = VecView(xi,viewer);CHKERRQ(ierr); }
#endif
    }
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }

  /*
     Free work space
  */
  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = VecDestroy(&xr);CHKERRQ(ierr);
  ierr = VecDestroy(&xi);CHKERRQ(ierr);
  if (nini > 0) {
    ierr = VecDestroyVecs(nini,&Iv);CHKERRQ(ierr);
  }
  if (ncon > 0) {
    ierr = VecDestroyVecs(ncon,&Cv);CHKERRQ(ierr);
  }






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
    arg->u = (double *)vec_alloc((nx+1)*(ny+1),0,sizeof(double));
    f = fopen("u.txt","r");
    if (f==NULL) {
        printf("Error while opening the file.\n");
    }
    i=0;
    while (fscanf(f,"%lf \n",&(arg->u[i])) !=EOF) {
        i = i+1;
    }

    fclose(f);

    arg->dudx = (double *)vec_alloc((nx+1)*(ny+1),0,sizeof(double));
    f = fopen("dudx.txt","r");
    if (f==NULL) {
        printf("Error while opening the file.\n");
    }
    i=0;
    while (fscanf(f,"%lf \n",&(arg->dudx[i])) !=EOF) {
        i = i+1;
    }
   
    fclose(f);


    arg->dudy = (double *)vec_alloc((nx+1)*(ny+1),0,sizeof(double));
    f = fopen("dudy.txt","r");
    if (f==NULL) {
        printf("Error while opening the file.\n");
    }
    i=0;
    while (fscanf(f,"%lf \n",&(arg->dudy[i])) !=EOF) {
        i = i+1;
    }
    
    fclose(f);

    arg->v = (double *)vec_alloc((nx+1)*(ny+1),0,sizeof(double));
    f = fopen("v.txt","r");
    if (f==NULL) {
        printf("Error while opening the file.\n");
    }
    i=0;
    while (fscanf(f,"%lf \n",&(arg->v[i])) !=EOF) {
        i = i+1;
    }

    fclose(f);

    arg->dvdx = (double *)vec_alloc((nx+1)*(ny+1),0,sizeof(double));
    f = fopen("dvdx.txt","r");
    if (f==NULL) {
        printf("Error while opening the file.\n");
    }
    i=0;
    while (fscanf(f,"%lf \n",&(arg->dvdx[i])) !=EOF) {
        i = i+1;
    }
   
    fclose(f);


    arg->dvdy = (double *)vec_alloc((nx+1)*(ny+1),0,sizeof(double));
    f = fopen("dvdy.txt","r");
    if (f==NULL) {
        printf("Error while opening the file.\n");
    }
    i=0;
    while (fscanf(f,"%lf \n",&(arg->dvdy[i])) !=EOF) {
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
    arg->w = (double *)vec_alloc((nx+1)*(ny+1),0,sizeof(double));
    arg->dwdx = (double *)vec_alloc((nx+1)*(ny+1),0,sizeof(double));
    arg->dwdy = (double *)vec_alloc((nx+1)*(ny+1),0,sizeof(double));
    
    for (i=0; i<(nx+1)*(ny+1); i++) {
        arg->w[i] = 0; arg->dwdx[i] = 0; arg->dwdy[i] = 0;
        arg->rho[i] = 1.0; arg->drhodx[i] = 0; arg->drhody[i] = 0;
        arg->d2rhodx[i] = 0; arg->d2rhody[i] = 0;
    }

    printf("Base flow reading done from files\n");
    return 0;
}



