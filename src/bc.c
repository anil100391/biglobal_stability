#include <data.h>
#include <stdlib.h>


int set_bc(ndr_data_t *arg) {

    FILE *f;
    int nx = arg->nx, ny = arg->ny;
    PetscInt i,j,Ii,locIi,Istart, Iend;
    PetscErrorCode ierr;
    int eqn,bctype,nz;
    int bcinfo[3][4];
    char ch,var,name;
    PetscScalar *der, one=1.0, zero=0.0;
    int *col,k;

    f = fopen("bc.conf","r");
    while (ch = fgetc(f) != '\n') {
    
    }

    for (Ii=0; Ii<12; Ii++) {
        var = fgetc(f); fgetc(f); name=fgetc(f); fgetc(f);
        if (var=='u') i=0;
        else if (var=='v') i=1;
        else if (var=='w') i=2;
        if (name=='b') j=0;
        else if (name=='r') j=1;
        else if (name=='t') j=2;
        else if (name=='l') j=3;
        ch  = fgetc(f);
        bcinfo[i][j] = ch - '0'; fgetc(f);

        printf("%c\t%c\t%d\n",var,name,bcinfo[i][j]);
    }
    fclose(f);

    ierr = MatGetOwnershipRange(arg->A,&Istart,&Iend); CHKERRQ(ierr);

    for (Ii=Istart; Ii<Iend; Ii++) {
        eqn = Ii/((nx+1)*(ny+1));
        locIi = Ii%((nx+1)*(ny+1));

        if (eqn > 2) continue;

        i = locIi/(ny+1); j=locIi%(ny+1);
        if (i==0) {
            //This is the left boundary (id=3)
            bctype = bcinfo[eqn][3];
            if (bctype == 1) {
                nz = 1;
                ierr = MatSetValues(arg->A,1,&Ii,1,&Ii,&one,INSERT_VALUES); CHKERRQ(ierr);
            }
            else if (bctype == 2) {
                nz = (nx+1)*(ny+1);
                der = (PetscScalar *)vec_alloc((nx+1)*(ny+1),0,sizeof(PetscScalar));
                col = (int *)vec_alloc((nx+1)*(ny+1),0,sizeof(int));
                for (k=0; k<(nx+1)*(ny+1); k++) col[k] = eqn*(nx+1)*(ny+1) + k;
                d1dx((int)locIi,nx,ny,&(arg->trnfm),der);
                ierr = MatSetValues(arg->A,1,&Ii,nz,col,der,INSERT_VALUES); CHKERRQ(ierr);
                free_vec((void *)col,0,sizeof(int));
                free_vec((void *)der,0,sizeof(PetscScalar));
            }
        }
        else if (i==nx) {
            //This is the right boundary
            bctype = bcinfo[eqn][1];
            if (bctype == 1) {
                nz = 1;
                ierr = MatSetValues(arg->A,1,&Ii,1,&Ii,&one,INSERT_VALUES); CHKERRQ(ierr);
            }
            else if(bctype == 2) {
                nz = (nx+1)*(ny+1);
                der = (PetscScalar *)vec_alloc((nx+1)*(ny+1),0,sizeof(PetscScalar));
                col = (int *)vec_alloc((nx+1)*(ny+1),0,sizeof(int));
                for (k=0; k<(nx+1)*(ny+1); k++) col[k] = eqn*(nx+1)*(ny+1) + k;
                d1dx((int)locIi,nx,ny,&(arg->trnfm),der);
                ierr = MatSetValues(arg->A,1,&Ii,nz,col,der,INSERT_VALUES); CHKERRQ(ierr);
                free_vec((void *)col,0,sizeof(int));
                free_vec((void *)der,0,sizeof(PetscScalar));
            }

        }
        else if (j==0) {
            //This is the bottom boundary
            bctype = bcinfo[eqn][0];
            if (bctype == 1) {
                nz = 1;
                ierr = MatSetValues(arg->A,1,&Ii,1,&Ii,&one,INSERT_VALUES); CHKERRQ(ierr);
            }
            else if(bctype == 2) {
                nz = (nx+1)*(ny+1);
                der = (PetscScalar *)vec_alloc((nx+1)*(ny+1),0,sizeof(PetscScalar));
                col = (int *)vec_alloc((nx+1)*(ny+1),0,sizeof(int));
                for (k=0; k<(nx+1)*(ny+1); k++) col[k] = eqn*(nx+1)*(ny+1) + k;
                d1dy((int)locIi,nx,ny,&(arg->trnfm),der);
                ierr = MatSetValues(arg->A,1,&Ii,nz,col,der,INSERT_VALUES); CHKERRQ(ierr);
                free_vec((void *)col,0,sizeof(int));
                free_vec((void *)der,0,sizeof(PetscScalar));
            }

        }
        else if (j==ny) {
            //This is the top boundary
            bctype = bcinfo[eqn][2];
            if (bctype == 1) {
                nz = 1;
                ierr = MatSetValues(arg->A,1,&Ii,1,&Ii,&one,INSERT_VALUES); CHKERRQ(ierr);
            }
            else if(bctype == 2) {
                nz = (nx+1)*(ny+1);
                der = (PetscScalar *)vec_alloc((nx+1)*(ny+1),0,sizeof(PetscScalar));
                col = (int *)vec_alloc((nx+1)*(ny+1),0,sizeof(int));
                for (k=0; k<(nx+1)*(ny+1); k++) col[k] = eqn*(nx+1)*(ny+1) + k;
                d1dy((int)locIi,nx,ny,&(arg->trnfm),der);
                ierr = MatSetValues(arg->A,1,&Ii,nz,col,der,INSERT_VALUES); CHKERRQ(ierr);
                free_vec((void *)col,0,sizeof(int));
                free_vec((void *)der,0,sizeof(PetscScalar));
            }

        }    
                
    }

    ierr = MatAssemblyBegin(arg->A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(arg->A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyBegin(arg->B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(arg->B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    return 0;

}


