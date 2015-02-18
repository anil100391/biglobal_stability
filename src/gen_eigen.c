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

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Load the matrices that define the eigensystem, Ax=kBx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nGeneralized eigenproblem.\n\n");CHKERRQ(ierr);

  ierr = MatGetVecs(arg->A,NULL,&xr);CHKERRQ(ierr);
  ierr = MatGetVecs(arg->A,NULL,&xi);CHKERRQ(ierr);

  /*
     Read user constraints if available
  */
  ierr = PetscOptionsGetInt(NULL,"-nconstr",&ncon,&flg);CHKERRQ(ierr);
  if (flg) {
    if (ncon<=0) SETERRQ(PETSC_COMM_WORLD,1,"The number of constraints must be >0");
    ierr = PetscOptionsGetString(NULL,"-fconstr",filename,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
    if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must specify the name of the file storing the constraints");
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
    ierr = VecDuplicateVecs(xr,ncon,&Cv);CHKERRQ(ierr);
    for (i=0;i<ncon;i++) {
      ierr = VecLoad(Cv[i],viewer);CHKERRQ(ierr);
    }
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }

  /*
     Read initial guesses if available
  */
  ierr = PetscOptionsGetInt(NULL,"-ninitial",&nini,&flg);CHKERRQ(ierr);
  if (flg) {
    if (nini<=0) SETERRQ(PETSC_COMM_WORLD,1,"The number of initial vectors must be >0");
    ierr = PetscOptionsGetString(NULL,"-finitial",filename,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
    if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must specify the name of the file containing the initial vectors");
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
    ierr = VecDuplicateVecs(xr,nini,&Iv);CHKERRQ(ierr);
    for (i=0;i<nini;i++) {
      ierr = VecLoad(Iv[i],viewer);CHKERRQ(ierr);
    }
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create eigensolver context
  */
  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
  ierr = EPSGetST(eps, &st);
  ierr = STSetType(st,STSINVERT);
  ierr = STSetShift(st,0.0001);
  /*
     Set operators. In this case, it is a generalized eigenvalue problem
  */
  ierr = EPSSetOperators(eps,arg->A,arg->B);CHKERRQ(ierr);

  /*
     If the user provided initial guesses or constraints, pass them here
  */
  ierr = EPSSetInitialSpace(eps,nini,Iv);CHKERRQ(ierr);
  ierr = EPSSetDeflationSpace(eps,ncon,Cv);CHKERRQ(ierr);

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
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
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

