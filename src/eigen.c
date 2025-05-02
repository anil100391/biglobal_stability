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

// Taken from slepsc tutorials ex7.c

#include <data.h>
#include <slepceps.h>

int eigen_solver(ndr_data_t *arg)
{
  EPS            eps;             /* eigenproblem solver context */
  ST             st;
  KSP            ksp;
  EPSType        type;
  PetscReal      tol;
  Vec            xr,xi,*Iv,*Cv;
  PetscInt       nev,maxit,i,its,lits,nconv,nini=0,ncon=0;
  char           filename[PETSC_MAX_PATH_LEN];
  PetscViewer    viewer;
  PetscBool      evecs,ishermitian;
  PetscMPIInt    rank,size,Istart,Iend;
  PetscScalar    kr,ki;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Load the matrices that define the eigensystem, Ax=kBx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nGeneralized eigenproblem.\n\n"));

  PetscCall(MatCreateVecs(arg->A,NULL,&xr));
  PetscCall(MatCreateVecs(arg->A,NULL,&xi));
  PetscCall(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));
  PetscCall(MPI_Comm_size(PETSC_COMM_WORLD,&size));
  PetscCall(MatGetOwnershipRange(arg->A,&Istart,&Iend));
  printf("process %d is responsible for rows %d to %d \n",rank,Istart,Iend);
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create eigensolver context
  */
  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));
  //PetscCall(EPSGetST(eps, &st);
  //PetscCall(STSetType(st,STSINVERT);
  //PetscCall(STSetShift(st,4.0+4.0*I);
  /*
     Set operators. In this case, it is a generalized eigenvalue problem
  */
  PetscCall(EPSSetOperators(eps,arg->A,arg->B));
  PetscCall(EPSSetProblemType(eps, EPS_GNHEP));
  //PetscCall(EPSSetTolerances(EPS eps,PetscReal tol,PetscInt maxits); CHKERRQ(ierr);

  /*
     Set solver parameters at runtime
  */
  PetscCall(EPSSetFromOptions(eps));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSSolve(eps));

  /*
     Optional: Get some information from the solver and display it
  */
  PetscCall(EPSGetIterationNumber(eps,&its));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %" PetscInt_FMT "\n",its));
  PetscCall(EPSGetST(eps,&st));
  PetscCall(STGetKSP(st,&ksp));
  PetscCall(KSPGetTotalIterations(ksp,&lits));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of linear iterations of the method: %" PetscInt_FMT "\n",lits));
  PetscCall(EPSGetType(eps,&type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type));
  PetscCall(EPSGetDimensions(eps,&nev,NULL,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %" PetscInt_FMT "\n",nev));
  PetscCall(EPSGetTolerances(eps,&tol,&maxit));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%" PetscInt_FMT "\n",(double)tol,maxit));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

//  PetscCall(EPSPrintSolution(eps,NULL));
  /*
     Save eigenvectors, if requested
  */
  PetscCall(PetscOptionsGetString(NULL, NULL,"-evecs",filename,PETSC_MAX_PATH_LEN,&evecs));
  PetscCall(EPSGetConverged(eps,&nconv));

  for (i=0; i<nconv; i++) {
    PetscCall(EPSGetEigenpair(eps,i,&kr,&ki,xr,xi));
    PetscPrintf(PETSC_COMM_WORLD,"%f\t %f\n",PetscRealPart(kr),PetscImaginaryPart(kr));
  }


  if (nconv>0 && evecs) {
    PetscCall(PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer));
    //PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&viewer));
    PetscCall(EPSIsHermitian(eps,&ishermitian));
    for (i=0;i<nconv;i++) {
      PetscCall(EPSGetEigenvector(eps,i,xr,xi));
      PetscCall(VecView(xr,viewer));
#if !defined(PETSC_USE_COMPLEX)
      if (!ishermitian) { PetscCall(VecView(xi,viewer)); }
#endif
    }
    PetscCall(PetscViewerDestroy(&viewer));
  }

  /*
     Free work space
  */
  PetscCall(EPSDestroy(&eps));
  PetscCall(VecDestroy(&xr));
  PetscCall(VecDestroy(&xi));
  if (nini > 0) {
    PetscCall(VecDestroyVecs(nini,&Iv));
  }
  if (ncon > 0) {
    PetscCall(VecDestroyVecs(ncon,&Cv));
  }
  return 0;
}

