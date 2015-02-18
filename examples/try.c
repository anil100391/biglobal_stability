// This program implements orr-sommerfeld stability equations for pipe flow.
static char help[] = "help yourself";
#include <stdio.h>
#include <math.h>
#include <slepceps.h>
#include <data.h>

#undef __FUNCT__
#define __FUNCTI__ "main"
int main(int argc, char **args)
{
    ndr_data_t ndrdata;
    PetscErrorCode ierr;
    PetscViewer plotmat;
    PetscViewer writemat;
    PetscScalar alpha, omega;
    int M=3;
    PetscInt Ii,Istart,Iend;
    PetscScalar a[3][3]={{1,2,1},
                         {3,1,4},
                         {1,3,5}};
    PetscScalar b[3][3]= {{1,2,0},
                          {0,1,0},
                          {1,0,1}};
    int col[3] = {0,1,2};
    int row;
    
    SlepcInitialize(&argc,&args,(char*)0,help);

    ierr = MatCreate(PETSC_COMM_WORLD,&(ndrdata.B)); CHKERRQ(ierr);
    ierr = MatSetSizes(ndrdata.B,PETSC_DECIDE,PETSC_DECIDE,M,M); CHKERRQ(ierr);
    ierr = MatSetFromOptions(ndrdata.B); CHKERRQ(ierr);
    ierr = MatSetUp(ndrdata.B); CHKERRQ(ierr);

    ierr = MatCreate(PETSC_COMM_WORLD,&(ndrdata.A)); CHKERRQ(ierr);
    ierr = MatSetSizes(ndrdata.A,PETSC_DECIDE,PETSC_DECIDE,M,M); CHKERRQ(ierr);
    ierr = MatSetFromOptions(ndrdata.A); CHKERRQ(ierr);
    ierr = MatSetUp(ndrdata.A); CHKERRQ(ierr);

    ierr = MatGetOwnershipRange(ndrdata.B,&Istart,&Iend); CHKERRQ(ierr);

    for (Ii=0; Ii<3; Ii++) {
        ierr = MatSetValues(ndrdata.B,1,&Ii,3,col,&b[Ii][0],INSERT_VALUES); CHKERRQ(ierr);
    }

    ierr = MatGetOwnershipRange(ndrdata.A,&Istart,&Iend); CHKERRQ(ierr);

    for (Ii=0; Ii<3; Ii++) {
        ierr = MatSetValues(ndrdata.A,1,&Ii,3,col,&a[Ii][0],INSERT_VALUES); CHKERRQ(ierr);
    }
   
    ierr = MatAssemblyBegin(ndrdata.A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(ndrdata.A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyBegin(ndrdata.B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(ndrdata.B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

 
    ierr = PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,"A",0,0,800,800,&plotmat);
//    MatView(ndrdata.B,plotmat);
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




