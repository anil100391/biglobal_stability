/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2014, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.

   SLEPc is free software: you can redistribute it and/or modify it under  the
   terms of version 3 of the GNU Lesser General Public License as published by
   the Free Software Foundation.

   SLEPc  is  distributed in the hope that it will be useful, but WITHOUT  ANY
   WARRANTY;  without even the implied warranty of MERCHANTABILITY or  FITNESS
   FOR  A  PARTICULAR PURPOSE. See the GNU Lesser General Public  License  for
   more details.

   You  should have received a copy of the GNU Lesser General  Public  License
   along with SLEPc. If not, see <http://www.gnu.org/licenses/>.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Solves a generalized eigensystem Ax=kBx with matrices passed as arguments.\n"
  "The command line options are:\n"
  "  -evecs <filename>, output file to save computed eigenvectors.\n"
  "  -ninitial <nini>, number of user-provided initial guesses.\n"
  "  -finitial <filename>, binary file containing <nini> vectors.\n"
  "  -nconstr <ncon>, number of user-provided constraints.\n"
  "  -fconstr <filename>, binary file containing <ncon> vectors.\n\n";

#include <data.h>
#include <slepceps.h>

int eigen_solver(ndr_data_t *arg)
{
  EPS            eps;             /* eigenproblem solver context */
  ST             st;
  EPSType        type;
  PetscReal      tol;
  Vec            xr,xi,*Iv,*Cv;
  PetscInt       nev,maxit,i,its,lits,nconv,nini=0,ncon=0;
  char           filename[PETSC_MAX_PATH_LEN];
  PetscViewer    viewer;
  PetscBool      flg,evecs,ishermitian;
  PetscErrorCode ierr;
  PetscMPIInt    rank,size,Istart,Iend;
  PetscScalar kr,ki;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Load the matrices that define the eigensystem, Ax=kBx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nGeneralized eigenproblem.\n\n");CHKERRQ(ierr);

  ierr = MatGetVecs(arg->A,NULL,&xr);CHKERRQ(ierr);
  ierr = MatGetVecs(arg->A,NULL,&xi);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);
  ierr = MatGetOwnershipRange(arg->A,&Istart,&Iend);
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
  ierr = EPSSetOperators(eps,arg->A,arg->B);CHKERRQ(ierr);
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

//  ierr = EPSPrintSolution(eps,NULL);CHKERRQ(ierr);
  /*
     Save eigenvectors, if requested
  */
  ierr = PetscOptionsGetString(NULL,"-evecs",filename,PETSC_MAX_PATH_LEN,&evecs);CHKERRQ(ierr);
  ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);

  for (i=0; i<nconv; i++) {
    ierr = EPSGetEigenpair(eps,i,&kr,&ki,xr,xi); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"%f\t %f\n",PetscRealPart(kr),PetscImaginaryPart(kr));
  }
    


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
  return 0;
}

