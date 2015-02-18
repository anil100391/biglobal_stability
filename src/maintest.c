// This program implements biglobal stability analysis.
static char help[] = "You are doomed";
#include <stdio.h>
#include <math.h>
#include <slepceps.h>
#include <data.h>
#include <time.h>


int main(int argc, char **args)
{
    //PetscErrorCode ierr;
    //SlepcInitialize(&argc,&args,(char*)0,help);
    //ierr = SlepcFinalize(); CHKERRQ(ierr);
    
    int nx=50, ny=50, M=(nx+1)*(ny+1),nz;
    int *col,col1,shift;
    int k,eqn;
    PetscErrorCode ierr;
    PetscViewer plotmat;
    Mat A;
    Vec v1,v2;
    double x,y;
    PetscMPIInt rank,size;
    PetscScalar *derx,*dery,*val,val1;
    PetscInt i,j,Ii,locIi,Istart,Iend;
    clock_t begin,end;
    double time_spent; 

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



    ierr = VecCreate(PETSC_COMM_WORLD,&v1); CHKERRQ(ierr);
    ierr = VecSetSizes(v1,PETSC_DECIDE,(nx+1)*(ny+1)); CHKERRQ(ierr);
    ierr = VecSetFromOptions(v1); CHKERRQ(ierr);


    ierr = VecGetOwnershipRange(v1,&Istart,&Iend);  CHKERRQ(ierr);
    for (Ii=Istart; Ii<Iend; Ii++) {
        col1 =  Ii;
        i = Ii/(ny+1); j = Ii%(ny+1);
        x = cos(i*M_PI/nx); y = cos(j*M_PI/ny);
        val1 = exp(x+y);
        VecSetValues(v1,1,&col1,&val1,INSERT_VALUES); CHKERRQ(ierr);

    }
    VecAssemblyBegin(v1);
    VecAssemblyEnd(v1);

    
    begin = clock();
    for (k=0; k<50000; k++) {
        ierr = VecCreate(PETSC_COMM_WORLD,&v2); CHKERRQ(ierr);
        ierr = VecSetSizes(v2,PETSC_DECIDE,(nx+1)*(ny+1)); CHKERRQ(ierr);
        ierr = VecSetFromOptions(v2); CHKERRQ(ierr);
        ierr = MatMult(A,v1,v2);  CHKERRQ(ierr);

        VecAssemblyBegin(v2);
        VecAssemblyEnd(v2);
        VecDestroy(&v2);
    }
    end = clock();
    time_spent = (double)(end-begin) / CLOCKS_PER_SEC;
    printf("Time taken for 10000 matmul is %lf\n",time_spent);
    


    ierr = PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,"A",0,0,800,800,&plotmat);
//    VecView(v2,PETSC_VIEWER_STDOUT_WORLD);
//    MatView(A,plotmat);
//    ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"A",FILE_MODE_WRITE,&writemat);
//    MatView(data.A,writemat);
//    ierr = PetscViewerDestroy(&writemat); CHKERRQ(ierr);
    printf("All done\n");
//    scanf("%d",&eqn);
    ierr = SlepcFinalize(); CHKERRQ(ierr);
    
    return 0;
}




