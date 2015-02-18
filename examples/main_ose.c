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
    
    ndrdata.ny = 300;
    ndrdata.re = 10000;
    alpha = 1.0; 
    
    SlepcInitialize(&argc,&args,(char*)0,help);

    orr_sommerfeld(&alpha,&omega,&ndrdata);
     
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




