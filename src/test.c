// This program implements biglobal stability analysis.
static char help[] = "You are doomed";
#include <stdio.h>
#include <math.h>
#include <slepceps.h>
#include <data.h>


int main(int argc, char **args)
{
    int nx=20, ny=20, M=(nx+1)*(ny+1),nz;
    int *col,shift;
    int eqn;
    PetscErrorCode ierr;
    PetscViewer plotmat;
    Mat A;
    PetscMPIInt rank,size;
    PetscScalar *derx,*dery,*val;
    PetscInt i,j,Ii,locIi,Istart,Iend;

    SlepcInitialize(&argc,&args,(char*)0,help);
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);
    printf("Size of petsc world is %d\n",size);

    /*************************** start here ************************/
    
    ierr = MatCreate(PETSC_COMM_WORLD,&A); CHKERRQ(ierr);
    ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,M,M); CHKERRQ(ierr);
    ierr = MatSetFromOptions(A); CHKERRQ(ierr);
    nz = (ny+1);
    ierr = MatMPIAIJSetPreallocation(A,nz,NULL,nz,NULL); CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(A,nz,NULL); CHKERRQ(ierr);
    ierr = MatSetUp(A); CHKERRQ(ierr);

    ierr = MatGetOwnershipRange(A,&Istart,&Iend); CHKERRQ(ierr);

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    printf("proc %d is responsible for row %d to %d\n",rank,Istart,Iend);

    dery = (PetscScalar *)vec_alloc((nx+1)*(ny+1),0,sizeof(PetscScalar));
    col = (int *)vec_alloc((ny+1),0,sizeof(int));
    val = (PetscScalar *)vec_alloc((ny+1),0,sizeof(PetscScalar));

    for (Ii=Istart; Ii<Iend; Ii++) {
        d2dy((int)Ii,nx,ny,dery);
        j = 0;
        shift = Ii/(ny+1);
        for (i=0; i<(ny+1); i++) {
            col[j] = shift*(ny+1) + i;
            val[j] = dery[col[j]];
            j = j + 1;
        }
         ierr = MatSetValues(A,1,&Ii,(ny+1),col,val,INSERT_VALUES); CHKERRQ(ierr);
         
    }

    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    free_vec((void *)dery,0,sizeof(PetscScalar));
    free_vec((void *)col,0,sizeof(int));
    free_vec((void *)val,0,sizeof(PetscScalar));
    /*************************** end here ************************/    

    ierr = PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,"A",0,0,800,800,&plotmat);
    MatView(A,plotmat);
//    ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"A",FILE_MODE_WRITE,&writemat);
//    MatView(data.A,writemat);
//    ierr = PetscViewerDestroy(&writemat); CHKERRQ(ierr);
    printf("All done\n");
    scanf("%d",&eqn);
    ierr = SlepcFinalize(); CHKERRQ(ierr);
    
    return 0;
}




