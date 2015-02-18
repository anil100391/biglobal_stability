#include <data.h>

void get_shift(int,int,int *);

int cheb_bc_helmholtz(ndr_data_t *arg2) {

    PetscErrorCode ierr;
    int nx=arg2->grid.nx, ny=arg2->grid.ny;
    int i,j,k,pos;
    PetscInt Ii,Istart, Iend;
    PetscMPIInt rank;
    PetscScalar one=1, zero=0;

    // Set boundary conditions
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    ierr = MatGetOwnershipRange(arg2->A,&Istart,&Iend);
    printf("proc %d says I have rows %d to %d\n",rank,Istart,Iend);

    for (Ii=Istart; Ii<Iend; Ii++) {
        
    }

    for(i=0; i<nx+1; i=i+nx) {
        for (j=0; j<ny+1; j++) {

            for (k=i*(ny+1); k<i*(ny+1)+ny+1; k++) {
                pos = i*(ny+1) + j;
                if (k == pos)
                { ierr = MatSetValues(arg2->A,1,&pos,1,&k,&one,INSERT_VALUES); CHKERRQ(ierr);}
                else
                { ierr = MatSetValues(arg2->A,1,&pos,1,&k,&zero,INSERT_VALUES); CHKERRQ(ierr);}
            }

            if (i==0) {
                for (k=ny+1+j; k<(nx+1)*(ny+1); k=k+ny+1) {
                    pos = i*(ny+1) + j;
                    ierr = MatSetValues(arg2->A,1,&pos,1,&k,&zero,INSERT_VALUES); CHKERRQ(ierr);
                }
            }
            else if (i==nx) {
                for (k=0+j; k<nx*(ny+1); k=k+ny+1) {
                    pos = i*(ny+1) + j;
                    ierr = MatSetValues(arg2->A,1,&pos,1,&k,&zero,INSERT_VALUES); CHKERRQ(ierr);
                }
            }
        }
    }

    for (i=1; i<nx; i++) {
        for (j=0; j<ny+1; j=j+ny) {

            for (k=0+j; k<(nx+1)*(ny+1); k=k+ny+1) {
                pos = i*(ny+1) + j;
                ierr = MatSetValues(arg2->A,1,&pos,1,&k,&zero,INSERT_VALUES); CHKERRQ(ierr);
            }
            for (k=i*(ny+1); k<i*(ny+1)+ny+1; k++) {
                pos = i*(ny+1) + j;
                if (k == pos)
                { ierr = MatSetValues(arg2->A,1,&pos,1,&k,&one,INSERT_VALUES); CHKERRQ(ierr);}
                else
                { ierr = MatSetValues(arg2->A,1,&pos,1,&k,&zero,INSERT_VALUES); CHKERRQ(ierr);}
            }
        }
    }

//  Assemble matrix A
    ierr = MatAssemblyBegin(arg2->A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(arg2->A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    return 0;
}





/*int fd_bc_helmholtz(ndr_data_t *arg2) {

    PetscErrorCode ierr;
    int nx=arg2->grid.nx, ny=arg2->grid.ny;
    int qx=arg2->grid.qx, qy=arg2->grid.qy;
    int i,j,k,pos;
    int *siy, *six;
    PetscScalar one=1, zero=0;

    siy = (int *)vec_alloc(ny+1,0,sizeof(int));
    six = (int *)vec_alloc(nx+1,0,sizeof(int));
    get_shift(qy,ny,siy);
    get_shift(qx,nx,six);
    // Set boundary conditions

    for(i=0; i<nx+1; i=i+nx) {
        for (j=0; j<ny+1; j++) {

            for (k=i*(ny+1)+siy[j]; k<i*(ny+1)+siy[j]+qy+1; k++) {

                pos = i*(ny+1) + j;
                if (k == pos)
                { ierr = MatSetValues(arg2->A,1,&pos,1,&k,&one,INSERT_VALUES); CHKERRQ(ierr);}
                else
                { ierr = MatSetValues(arg2->A,1,&pos,1,&k,&zero,INSERT_VALUES); CHKERRQ(ierr);}
            }

            if (i==0) {
                for (k=ny+1+j+six[i]; k<(qx+1)*(ny+1); k=k+ny+1) {
                    pos = i*(ny+1) + j;
                    ierr = MatSetValues(arg2->A,1,&pos,1,&k,&zero,INSERT_VALUES); CHKERRQ(ierr);
                }
            }
            else if (i==nx) {
                for (k=0+j+six[i]*(ny+1); k<nx*(ny+1); k=k+ny+1) {
                    pos = i*(ny+1) + j;
                    ierr = MatSetValues(arg2->A,1,&pos,1,&k,&zero,INSERT_VALUES); CHKERRQ(ierr);
                }
            }
        }
    }

    for (i=1; i<nx; i++) {
        for (j=0; j<ny+1; j=j+ny) {

            for (k=0+j+six[i]*(ny+1); k<(qx+1)*(ny+1)+six[i]*(ny+1); k=k+ny+1) {
                pos = i*(ny+1) + j;
                ierr = MatSetValues(arg2->A,1,&pos,1,&k,&zero,INSERT_VALUES); CHKERRQ(ierr);
            }
            for (k=i*(ny+1)+siy[j]; k<i*(ny+1)+siy[j]+qy+1; k++) {
                pos = i*(ny+1) + j;
                if (k == pos)
                { ierr = MatSetValues(arg2->A,1,&pos,1,&k,&one,INSERT_VALUES); CHKERRQ(ierr);}
                else
                { ierr = MatSetValues(arg2->A,1,&pos,1,&k,&zero,INSERT_VALUES); CHKERRQ(ierr);}
            }
        }
    }
//  Assemble matrix A
    ierr = MatAssemblyBegin(arg2->A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(arg2->A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    return 0;
}
*/
