#include <data.h>

//****************************************************************************//
/*int helmholtz(PetscScalar *beta, PetscScalar *omega, ndr_data_t *arg)
{
    PetscErrorCode ierr;
    int nx = arg->nx, ny = arg->ny;
    int M = (nx+1)*(ny+1);
    int *col,shift,nz;
    PetscScalar *derx,*dery,*val;
    PetscInt i,j,Ii,Istart,Iend;
    PetscMPIInt rank;
    
    ierr = MatCreate(PETSC_COMM_WORLD,&(arg->A)); CHKERRQ(ierr);
    ierr = MatSetSizes(arg->A,PETSC_DECIDE,PETSC_DECIDE,M,M); CHKERRQ(ierr);
    ierr = MatSetFromOptions(arg->A); CHKERRQ(ierr);
    nz = (ny+nx+1);
    ierr = MatMPIAIJSetPreallocation(arg->A,nz,NULL,nz,NULL); CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(arg->A,nz,NULL); CHKERRQ(ierr);
    ierr = MatSetUp(arg->A); CHKERRQ(ierr);

    ierr = MatGetOwnershipRange(arg->A,&Istart,&Iend); CHKERRQ(ierr);

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    printf("proc %d is responsible for row %d to %d\n",rank,Istart,Iend);

    derx = (PetscScalar *)vec_alloc(M,0,sizeof(PetscScalar));
    dery = (PetscScalar *)vec_alloc(M,0,sizeof(PetscScalar));
    col = (int *)vec_alloc((nx+ny+1),0,sizeof(int));
    val = (PetscScalar *)vec_alloc((nx+ny+1),0,sizeof(PetscScalar));

    for (Ii=Istart; Ii<Iend; Ii++) {
        d2dx((int)Ii,nx,ny,derx);
        d2dy((int)Ii,nx,ny,dery);
        for(i=0; i<M; i++) {
            derx[i] = derx[i] + dery[i];
        }
        j = 0;
        shift = Ii/(ny+1);
        for (i=0; i<shift; i++) {
            col[j] = Ii%(ny+1) + i*(ny+1);
            val[j] = derx[col[j]];
            j = j + 1;
        }
        for (i=0; i<(ny+1); i++) {
            col[j] = shift*(ny+1) + i;
            val[j] = derx[col[j]];
            j = j + 1;
        }
        for (i=shift+1; i<(nx+1); i++) {
            col[j] = Ii%(ny+1) + i*(ny+1);
            val[j] = derx[col[j]];
            j = j + 1;
        }
        ierr = MatSetValues(arg->A,1,&Ii,(nx+ny+1),col,val,INSERT_VALUES); CHKERRQ(ierr);
         
    }

    free_vec((void *)derx,0,sizeof(PetscScalar));
    free_vec((void *)dery,0,sizeof(PetscScalar));
    free_vec((void *)col,0,sizeof(int));
    free_vec((void *)val,0,sizeof(PetscScalar));
    
    /////////////////// Boundary Conditions //////////////////////

    col = (int *)vec_alloc((nx+ny+1),0,sizeof(int));
    val = (PetscScalar *)vec_alloc((nx+ny+1),0,sizeof(PetscScalar));
    for (Ii=Istart; Ii<Iend; Ii++) {
        i = Ii/(ny+1); j = Ii%(ny+1);
        if (i*j==0 || i==nx || j==ny) {
            j = 0;
            shift = Ii/(ny+1);
            for (i=0; i<shift; i++) {
                col[j] = Ii%(ny+1) + i*(ny+1);
                val[j] = 0.0;
                j = j + 1;
            }
            for (i=0; i<(ny+1); i++) {
                col[j] = shift*(ny+1) + i;
                if (col[j]==Ii)
                    val[j] = 1.0;
                else
                    val[j] = 0.0;
                j = j + 1;
            }
            for (i=shift+1; i<(nx+1); i++) {
                col[j] = Ii%(ny+1) + i*(ny+1);
                val[j] = 0.0;
                j = j + 1;
            }
            ierr = MatSetValues(arg->A,1,&Ii,(nx+ny+1),col,val,INSERT_VALUES); CHKERRQ(ierr);
        }
    }


    ierr = MatAssemblyBegin(arg->A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(arg->A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    free_vec((void *)col,0,sizeof(int));
    free_vec((void *)val,0,sizeof(PetscScalar));

  
    return 0;
}
*/


int lns_biglobal(int flag, PetscScalar *beta, PetscScalar *omega, ndr_data_t *arg) 
{
    PetscErrorCode ierr;
    int nx = arg->nx, ny = arg->ny;
    int M = 5*(nx+1)*(ny+1),nz;
    int *col,shift;
    int eqn;
    double rho,u,v,w;
    PetscScalar *derx,*dery,*val;
    PetscInt i,j,Ii,locIi,Istart,Iend;
    PetscMPIInt rank;
    PetscScalar one=1.0,zero=0.0,art=1e-12;
    
    ierr = MatCreate(PETSC_COMM_WORLD,&(arg->A)); CHKERRQ(ierr);
    ierr = MatSetSizes(arg->A,PETSC_DECIDE,PETSC_DECIDE,M,M); CHKERRQ(ierr);
    ierr = MatSetFromOptions(arg->A); CHKERRQ(ierr);
    if (ny>nx) nz = (ny+nx+1)+(nx+1)*(ny+1)+2*(ny+1)+1;
    else  nz = (ny+nx+1)+(nx+1)*(ny+1)+2*(nx+1)+1;
    ierr = MatMPIAIJSetPreallocation(arg->A,nz,NULL,nz,NULL); CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(arg->A,nz,NULL); CHKERRQ(ierr);
    ierr = MatSetUp(arg->A); CHKERRQ(ierr);

    ierr = MatGetOwnershipRange(arg->A,&Istart,&Iend); CHKERRQ(ierr);

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
//    printf("Constructing A\n");
//    printf("proc %d is responsible for row %d to %d\n",rank,Istart,Iend);

    derx = (PetscScalar *)vec_alloc(M,0,sizeof(PetscScalar));
    dery = (PetscScalar *)vec_alloc(M,0,sizeof(PetscScalar));


    for (Ii=Istart; Ii<Iend; Ii++) {
        eqn = Ii/((nx+1)*(ny+1));
        locIi = Ii%((nx+1)*(ny+1));
        rho = arg->rho[locIi];
        u = arg->u[locIi];
        v = arg->v[locIi];
        w = arg->w[locIi];
        switch (eqn) 
        {
        case 0:

            /////////////////// Boundary Conditions //////////////////////
            i = locIi/(ny+1); j=locIi%(ny+1);
            nz = 1;
            if (i*j==0 || i==nx || j==ny) {
                //ierr = MatSetValues(arg->A,1,&Ii,nz,&Ii,&one,INSERT_VALUES); CHKERRQ(ierr);
                continue;
            }
            else {
    
            nz = (nx+ny+1)+(nx+1)*(ny+1)+(nx+1)+1+(nx+1);
            col = (int *)vec_alloc(nz,0,sizeof(int));
            val = (PetscScalar *)vec_alloc(nz,0,sizeof(PetscScalar));
            ////////////// x-momentum block 1 ///////////////
            d1dx((int)locIi,nx,ny,&(arg->trnfm),derx);
            for (i=0; i<(nx+1)*(ny+1); i++) 
                derx[i] = rho*u*derx[i];
            d1dy((int)locIi,nx,ny,&(arg->trnfm),dery);
            for (i=0; i<(nx+1)*(ny+1); i++) 
                dery[i] = rho*v*dery[i];
            for(i=0; i<(nx+1)*(ny+1); i++)
                derx[i] = derx[i] + dery[i];
            derx[locIi] = derx[locIi] + (*beta)*rho*w*I + flag*(*beta)*(*beta)/arg->re + rho*arg->dudx[locIi];
            d2dx((int)locIi,nx,ny,&(arg->trnfm),dery);
            for(i=0; i<(nx+1)*(ny+1); i++)
                derx[i] = derx[i] - 4.0*flag*dery[i]/(3.0*arg->re);
            d2dy((int)locIi,nx,ny,&(arg->trnfm),dery);
            for(i=0; i<(nx+1)*(ny+1); i++)
                derx[i] = derx[i] - flag*dery[i]/arg->re;
            
            j = 0;
            shift = locIi/(ny+1);
            for (i=0; i<shift; i++) {
                col[j] = locIi%(ny+1) + i*(ny+1);
                val[j] = derx[col[j]];
                j = j + 1;
            }
            for (i=0; i<(ny+1); i++) {
                col[j] = shift*(ny+1) + i;
                val[j] = derx[col[j]];
                j = j + 1;
            }
            for (i=shift+1; i<(nx+1); i++) {
                col[j] = locIi%(ny+1) + i*(ny+1);
                val[j] = derx[col[j]];
                j = j + 1;
            }
            ////////////// x-momentum block 2 ///////////////
            d2dxy((int)locIi,nx,ny,&(arg->trnfm),derx);
            for (i=0; i<(nx+1)*(ny+1); i++)
                derx[i] = -(1.0*flag/(3.0*arg->re))*derx[i];
            derx[locIi] = derx[locIi] + rho*arg->dudy[locIi];
            
            for (i=0; i<(nx+1)*(ny+1); i++) {
                col[j] = (nx+1)*(ny+1) + i;
                val[j] = derx[i];
                j = j + 1;
            }
            ////////////// x-momentum block 3 ///////////////
            d1dx((int)locIi,nx,ny,&(arg->trnfm),derx);
            for (i=0; i<(nx+1)*(ny+1); i++)
                derx[i] = -((*beta)*flag*derx[i]/(3.0*arg->re))*I;
            for (i=0; i<(nx+1); i++) {
                col[j] = 2*(nx+1)*(ny+1) + locIi%(ny+1) + i*(ny+1);
                val[j] = derx[col[j]-2*(nx+1)*(ny+1)];
                j = j + 1;
            }
            ////////////// x-momentum block 4 ///////////////
            col[j] = 3*(nx+1)*(ny+1) + locIi;
            val[j] = u*arg->dudx[locIi] + v*arg->dudy[locIi];
            j = j + 1;
            ////////////// x-momentum block 5 ///////////////

            d1dx((int)locIi,nx,ny,&(arg->trnfm),derx);
            for (i=0; i<(nx+1)*(ny+1); i++)
                derx[i] = derx[i]/arg->gamma;
            for (i=0; i<(nx+1); i++) {
                col[j] = 4*(nx+1)*(ny+1) + locIi%(ny+1) + i*(ny+1);
                val[j] = derx[col[j]-4*(nx+1)*(ny+1)];
                j = j + 1;
            }
            nz = (ny+nx+1)+(ny+1)*(nx+1)+(nx+1)+1+(nx+1);
            ierr = MatSetValues(arg->A,1,&Ii,nz,col,val,INSERT_VALUES); CHKERRQ(ierr);
            free_vec((void *)col,0,sizeof(int));
            free_vec((void *)val,0,sizeof(PetscScalar));
            }
            break; 
        
        case 1:

            /////////////////// Boundary Conditions //////////////////////
            i = locIi/(ny+1); j=locIi%(ny+1);
            nz = 1;
            if (i*j==0 || i==nx || j==ny) {
                //ierr = MatSetValues(arg->A,1,&Ii,nz,&Ii,&one,INSERT_VALUES); CHKERRQ(ierr);
                continue;
            }
            else {

            nz = (nx+ny+1)+(nx+1)*(ny+1)+(ny+1)+1+(ny+1);
            col = (int *)vec_alloc(nz,0,sizeof(int));
            val = (PetscScalar *)vec_alloc(nz,0,sizeof(PetscScalar));
            ////////////// y-momentum block 1 ///////////////
            d2dxy((int)locIi,nx,ny,&(arg->trnfm),derx);
            for (i=0; i<(nx+1)*(ny+1); i++)
                derx[i] = (-1.0*flag/(3.0*arg->re))*derx[i];
            derx[locIi] = derx[locIi] + rho*arg->dvdx[locIi];
            j = 0;
            shift = locIi/(ny+1);
            for (i=0; i<(nx+1)*(ny+1); i++) {
                col[j] = i;
                val[j] = derx[col[j]];
                j = j + 1;
            }
            ////////////// y-momentum block 2 ///////////////
            d1dx((int)locIi,nx,ny,&(arg->trnfm),derx);
            for (i=0; i<(nx+1)*(ny+1); i++) 
                derx[i] = rho*u*derx[i];
            d1dy((int)locIi,nx,ny,&(arg->trnfm),dery);
            for (i=0; i<(nx+1)*(ny+1); i++) 
                dery[i] = rho*v*dery[i];
            for(i=0; i<(nx+1)*(ny+1); i++)
                derx[i] = derx[i] + dery[i];
            derx[locIi] = derx[locIi] + (*beta)*rho*w*I + flag*(*beta)*(*beta)/arg->re + rho*arg->dvdy[locIi];
            d2dx((int)locIi,nx,ny,&(arg->trnfm),dery);
            for(i=0; i<(nx+1)*(ny+1); i++)
                derx[i] = derx[i] - flag*dery[i]/arg->re;
            d2dy((int)locIi,nx,ny,&(arg->trnfm),dery);
            for(i=0; i<(nx+1)*(ny+1); i++)
                derx[i] = derx[i] -4.0*flag*dery[i]/(3.0*arg->re);

            shift = locIi/(ny+1);
            for (i=0; i<shift; i++) {
                col[j] = (nx+1)*(ny+1) + locIi%(ny+1) + i*(ny+1);
                val[j] = derx[col[j]-(nx+1)*(ny+1)];
                j = j + 1;
            }
            for (i=0; i<(ny+1); i++) {
                col[j] = (nx+1)*(ny+1)+shift*(ny+1) + i;
                val[j] = derx[col[j]-(nx+1)*(ny+1)];
                j = j + 1;
            }
            for (i=shift+1; i<(nx+1); i++) {
                col[j] = (nx+1)*(ny+1) + locIi%(ny+1) + i*(ny+1);
                val[j] = derx[col[j]-(nx+1)*(ny+1)];
                j = j + 1;
            }
            
            ////////////// y-momentum block 3 ///////////////
            d1dy((int)locIi,nx,ny,&(arg->trnfm),dery);
            for (i=0; i<(nx+1)*(ny+1); i++)
                dery[i] = -((*beta)*flag*dery[i]/(3.0*arg->re))*I;
            for (i=0; i<(ny+1); i++) {
                col[j] = 2*(nx+1)*(ny+1) + shift*(ny+1) + i;
                val[j] = dery[col[j]-2*(nx+1)*(ny+1)];
                j = j + 1;
            }
            ////////////// y-momentum block 4 ///////////////
            col[j] = 3*(nx+1)*(ny+1) + locIi;
            val[j] = u*arg->dvdx[locIi] + v*arg->dvdy[locIi];
            j = j + 1;
            
            ////////////// y-momentum block 5 ///////////////
            d1dy((int)locIi,nx,ny,&(arg->trnfm),dery);
            for (i=0; i<(nx+1)*(ny+1); i++)
                dery[i] = dery[i]/arg->gamma;
            for (i=0; i<(ny+1); i++) {
                col[j] = 4*(nx+1)*(ny+1) + shift*(ny+1) + i;
                val[j] = dery[col[j]-4*(nx+1)*(ny+1)];
                j = j + 1;
            }
            nz = (ny+1)*(nx+1)+(ny+nx+1)+(ny+1)+1+(ny+1);
            ierr = MatSetValues(arg->A,1,&Ii,nz,col,val,INSERT_VALUES); CHKERRQ(ierr);
            free_vec((void *)col,0,sizeof(int));
            free_vec((void *)val,0,sizeof(PetscScalar));
            }
            break;

        case 2:
            /////////////////// Boundary Conditions //////////////////////
            i = locIi/(ny+1); j=locIi%(ny+1);
            nz = 1;
            if (i*j==0 || i==nx || j==ny) {
                //ierr = MatSetValues(arg->A,1,&Ii,nz,&Ii,&one,INSERT_VALUES); CHKERRQ(ierr);
                continue;
            }
            else {

            nz = (nx+ny+1)+(nx+1)*(ny+1)+(ny+1)+1+1;
            col = (int *)vec_alloc(nz,0,sizeof(int));
            val = (PetscScalar *)vec_alloc(nz,0,sizeof(PetscScalar));
            ////////////// z-momentum block 1 ///////////////
            d1dx((int)locIi,nx,ny,&(arg->trnfm),derx);
            for (i=0; i<(nx+1)*(ny+1); i++)
                derx[i] = (-1.0*flag/(3.0*arg->re))*derx[i]*I;
            derx[locIi] = derx[locIi] + rho*arg->dwdx[locIi];
            j = 0;
            shift = locIi/(ny+1);
            for (i=0; i<(nx+1); i++) {
                col[j] = locIi%(ny+1) + i*(ny+1);
                val[j] = derx[col[j]];
                j = j + 1;
            }
            
            ////////////// z-momentum block 2 ///////////////
            d1dy((int)locIi,nx,ny,&(arg->trnfm),dery);
            for (i=0; i<(nx+1)*(ny+1); i++)
                dery[i] = (-1.0*(*beta)*flag/(3.0*arg->re))*dery[i]*I;
            dery[locIi] = dery[locIi] + rho*arg->dwdy[locIi];
            for (i=0; i<(ny+1); i++) {
                col[j] = (nx+1)*(ny+1) + shift*(ny+1) + i;
                val[j] = dery[col[j]-(nx+1)*(ny+1)];
                j = j + 1;
            }
            ////////////// z-momentum block 3 ///////////////
            d1dx((int)locIi,nx,ny,&(arg->trnfm),derx);
            for (i=0; i<(nx+1)*(ny+1); i++) 
                derx[i] = rho*u*derx[i];
            d1dy((int)locIi,nx,ny,&(arg->trnfm),dery);
            for (i=0; i<(nx+1)*(ny+1); i++) 
                dery[i] = rho*v*dery[i];
            for(i=0; i<(nx+1)*(ny+1); i++)
                derx[i] = derx[i] + dery[i];
            derx[locIi] = derx[locIi] + (*beta)*rho*w*I + flag*(*beta)*(*beta)/arg->re;
            d2dx((int)locIi,nx,ny,&(arg->trnfm),dery);
            for(i=0; i<(nx+1)*(ny+1); i++)
                derx[i] = derx[i] - flag*dery[i]/arg->re;
            d2dy((int)locIi,nx,ny,&(arg->trnfm),dery);
            for(i=0; i<(nx+1)*(ny+1); i++)
                derx[i] = derx[i] - flag*dery[i]/arg->re;
            derx[locIi] = derx[locIi] + flag*(*beta)*(*beta)/(3.0*arg->re);
            shift = locIi/(ny+1);
            for (i=0; i<shift; i++) {
                col[j] = 2*(nx+1)*(ny+1) + locIi%(ny+1) + i*(ny+1);
                val[j] = derx[col[j]-2*(nx+1)*(ny+1)];
                j = j + 1;
            }
            for (i=0; i<(ny+1); i++) {
                col[j] = 2*(nx+1)*(ny+1)+shift*(ny+1) + i;
                val[j] = derx[col[j]-2*(nx+1)*(ny+1)];
                j = j + 1;
            }
            for (i=shift+1; i<(nx+1); i++) {
                col[j] = 2*(nx+1)*(ny+1) + locIi%(ny+1) + i*(ny+1);
                val[j] = derx[col[j]-2*(nx+1)*(ny+1)];
                j = j + 1;
            }

            ////////////// z-momentum block 4 ///////////////
            col[j] = 3*(nx+1)*(ny+1) + locIi;
            val[j] = u*arg->dwdx[locIi] + v*arg->dwdy[locIi];
            j = j + 1;

            ////////////// z-momentum block 5 ///////////////
            col[j] = 4*(nx+1)*(ny+1) + locIi;
            val[j] = ((*beta)/arg->gamma)*I;
            j = j + 1;

            nz = (nx+1)+(ny+1)+(ny+nx+1)+1+1;
            ierr = MatSetValues(arg->A,1,&Ii,nz,col,val,INSERT_VALUES); CHKERRQ(ierr);
            free_vec((void *)col,0,sizeof(int));
            free_vec((void *)val,0,sizeof(PetscScalar));
            }
            break;

        case 3:
            nz = (nx+1) + (ny+1) + 1 + (ny+nx+1);
            col = (int *)vec_alloc(nz,0,sizeof(int));
            val = (PetscScalar *)vec_alloc(nz,0,sizeof(PetscScalar));
            ////////////// continuity eqn block 1 ///////////////
            d1dx((int)locIi,nx,ny,&(arg->trnfm),derx);
            for (i=0; i<(nx+1)*(ny+1); i++)
                derx[i] = rho*rho*derx[i];
            derx[locIi] = derx[locIi] + rho*arg->drhodx[locIi];
            j = 0;
            shift = locIi/(ny+1);
            for (i=0; i<(nx+1); i++) {
                col[j] = locIi%(ny+1) + i*(ny+1);
                val[j] = derx[col[j]];
                j = j + 1;
            }
            
            ////////////// continuity eqn block 2 ///////////////
            d1dy((int)locIi,nx,ny,&(arg->trnfm),dery);
            for (i=0; i<(nx+1)*(ny+1); i++)
                dery[i] = rho*rho*dery[i];
            dery[locIi] = dery[locIi] + rho*arg->drhody[locIi];
            for (i=0; i<(ny+1); i++) {
                col[j] = (nx+1)*(ny+1) + shift*(ny+1) + i;
                val[j] = dery[col[j]-(nx+1)*(ny+1)];
                j = j + 1;
            }
          
            ////////////// continuity eqn block 3 ///////////////
            col[j] = 2*(nx+1)*(ny+1) + locIi;
            val[j] = rho*rho*(*beta)*I;
            j = j + 1;
            ////////////// continuity eqn block 4 ///////////////
            d1dx((int)locIi,nx,ny,&(arg->trnfm),derx);
            for (i=0; i<(nx+1)*(ny+1); i++) 
                derx[i] = rho*u*derx[i];
            d1dy((int)locIi,nx,ny,&(arg->trnfm),dery);
            for (i=0; i<(nx+1)*(ny+1); i++) 
                dery[i] = rho*v*dery[i];
            for(i=0; i<(nx+1)*(ny+1); i++)
                derx[i] = derx[i] + dery[i];
            derx[locIi] = derx[locIi] + (*beta)*rho*w*I + rho*(arg->dudx[locIi] + arg->dvdy[locIi]);
            shift = locIi/(ny+1);
            for (i=0; i<shift; i++) {
                col[j] = 3*(nx+1)*(ny+1) + locIi%(ny+1) + i*(ny+1);
                val[j] = derx[col[j]-3*(nx+1)*(ny+1)];
                j = j + 1;
            }
            for (i=0; i<(ny+1); i++) {
                col[j] = 3*(nx+1)*(ny+1)+shift*(ny+1) + i;
                val[j] = derx[col[j]-3*(nx+1)*(ny+1)];
                j = j + 1;
            }
            for (i=shift+1; i<(nx+1); i++) {
                col[j] = 3*(nx+1)*(ny+1) + locIi%(ny+1) + i*(ny+1);
                val[j] = derx[col[j]-3*(nx+1)*(ny+1)];
                j = j + 1;
            }

            nz = (nx+1)+(ny+1)+1+(ny+nx+1);
            ierr = MatSetValues(arg->A,1,&Ii,nz,col,val,INSERT_VALUES); CHKERRQ(ierr);
            free_vec((void *)col,0,sizeof(int));
            free_vec((void *)val,0,sizeof(PetscScalar));
            break;
        
        case 4:
            nz = (nx+1) + (ny+1) + 1 + (ny+nx+1);
            col = (int *)vec_alloc(nz,0,sizeof(int));
            val = (PetscScalar *)vec_alloc(nz,0,sizeof(PetscScalar));
            ////////////// energy eqn block 1 ///////////////
            d1dx((int)locIi,nx,ny,&(arg->trnfm),derx);
            for (i=0; i<(nx+1)*(ny+1); i++)
                derx[i] = -arg->gamma*rho*rho*derx[i];
            j = 0;
            shift = locIi/(ny+1);
            for (i=0; i<(nx+1); i++) {
                col[j] = locIi%(ny+1) + i*(ny+1);
                val[j] = derx[col[j]];
                j = j + 1;
            }
            
            ////////////// energy eqn block 2 ///////////////
            d1dy((int)locIi,nx,ny,&(arg->trnfm),dery);
            for (i=0; i<(nx+1)*(ny+1); i++)
                dery[i] = -arg->gamma*rho*rho*dery[i];
            for (i=0; i<(ny+1); i++) {
                col[j] = (nx+1)*(ny+1) + shift*(ny+1) + i;
                val[j] = dery[col[j]-(nx+1)*(ny+1)];
                j = j + 1;
            }
          
            ////////////// energy eqn block 3 ///////////////
            col[j] = 2*(nx+1)*(ny+1) + locIi;
            val[j] = -arg->gamma*rho*rho*(*beta)*I;
            j = j + 1;
            ////////////// energy eqn block 4 ///////////////
            for (i=0; i<(nx+1)*(ny+1); i++) derx[i]=0;
            derx[locIi] = 2*arg->drhodx[locIi]*u + 2*arg->drhody[locIi]*v;
            derx[locIi] = derx[locIi] - (3*arg->gamma-2)*rho*(arg->dudx[locIi]+arg->dvdy[locIi]);
            derx[locIi] = derx[locIi] + (arg->gamma*flag/(arg->pr*arg->re))*((*beta)*(*beta) - arg->d2rhodx[locIi]/rho);
            derx[locIi] = derx[locIi] - (arg->gamma*flag/(arg->pr*arg->re))*(arg->d2rhody[locIi]/rho);
            d1dy((int)locIi,nx,ny,&(arg->trnfm),dery);
            for (i=0; i<(nx+1)*(ny+1); i++)
                derx[i] = derx[i] + (arg->gamma*flag/(arg->pr*arg->re))*(2*arg->drhody[locIi]/rho)*dery[i];
            d2dy((int)locIi,nx,ny,&(arg->trnfm),dery);
            for (i=0; i<(nx+1)*(ny+1); i++)
                derx[i] = derx[i] - (arg->gamma*flag/(arg->pr*arg->re))*dery[i];
            d1dx((int)locIi,nx,ny,&(arg->trnfm),dery);
            for (i=0; i<(nx+1)*(ny+1); i++)
                derx[i] = derx[i] + (arg->gamma*flag/(arg->pr*arg->re))*(2*arg->drhodx[locIi]/rho)*dery[i];
            d2dx((int)locIi,nx,ny,&(arg->trnfm),dery);
            for (i=0; i<(nx+1)*(ny+1); i++)
                derx[i] = derx[i] - (arg->gamma*flag/(arg->pr*arg->re))*dery[i];
            
            shift = locIi/(ny+1);
            for (i=0; i<shift; i++) {
                col[j] = 3*(nx+1)*(ny+1) + locIi%(ny+1) + i*(ny+1);
                val[j] = derx[col[j]-3*(nx+1)*(ny+1)];
                j = j + 1;
            }
            for (i=0; i<(ny+1); i++) {
                col[j] = 3*(nx+1)*(ny+1)+shift*(ny+1) + i;
                val[j] = derx[col[j]-3*(nx+1)*(ny+1)];
                j = j + 1;
            }
            for (i=shift+1; i<(nx+1); i++) {
                col[j] = 3*(nx+1)*(ny+1) + locIi%(ny+1) + i*(ny+1);
                val[j] = derx[col[j]-3*(nx+1)*(ny+1)];
                j = j + 1;
            }

            nz = (nx+1)+(ny+1)+1+(ny+nx+1);
            ierr = MatSetValues(arg->A,1,&Ii,nz,col,val,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValues(arg->A,1,&Ii,1,&Ii,&art,INSERT_VALUES); CHKERRQ(ierr);//artificial addition
            free_vec((void *)col,0,sizeof(int));
            free_vec((void *)val,0,sizeof(PetscScalar));
            break;    
        }

         
    }

    free_vec((void *)derx,0,sizeof(PetscScalar));
    free_vec((void *)dery,0,sizeof(PetscScalar));

    ///////////////////////// Construction of B //////////////////////////

    ierr = MatCreate(PETSC_COMM_WORLD,&(arg->B)); CHKERRQ(ierr);
    ierr = MatSetSizes(arg->B,PETSC_DECIDE,PETSC_DECIDE,M,M); CHKERRQ(ierr);
    ierr = MatSetFromOptions(arg->B); CHKERRQ(ierr);
    nz = 2;
    ierr = MatMPIAIJSetPreallocation(arg->B,nz,NULL,nz,NULL); CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(arg->B,nz,NULL); CHKERRQ(ierr);
    ierr = MatSetUp(arg->B); CHKERRQ(ierr);

    ierr = MatGetOwnershipRange(arg->B,&Istart,&Iend); CHKERRQ(ierr);

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
//    printf("Constructing B\n");
//    printf("proc %d is responsible for row %d to %d\n",rank,Istart,Iend);
    for (Ii=Istart; Ii<Iend; Ii++) {
        locIi = Ii%((nx+1)*(ny+1));
        eqn = Ii/((nx+1)*(ny+1));
        one = arg->rho[locIi]*I;
        switch (eqn) {
        case 3:
            ierr = MatSetValues(arg->B,1,&Ii,1,&Ii,&one,INSERT_VALUES); CHKERRQ(ierr);
            break;
        case 4:
            i = Ii;
            ierr = MatSetValues(arg->B,1,&Ii,1,&i,&art,INSERT_VALUES); CHKERRQ(ierr); //artifical addition
            break;
        default:
            i = locIi/(ny+1); j=locIi%(ny+1);
            if (i*j==0 || i==nx || j==ny) {
                ierr = MatSetValues(arg->B,1,&Ii,1,&Ii,&zero,INSERT_VALUES); CHKERRQ(ierr);
            }
            else {
                ierr = MatSetValues(arg->B,1,&Ii,1,&Ii,&one,INSERT_VALUES); CHKERRQ(ierr);
            }
        }
    }


    
    //ierr = MatAssemblyBegin(arg->A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    //ierr = MatAssemblyEnd(arg->A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    //ierr = MatAssemblyBegin(arg->B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    //ierr = MatAssemblyEnd(arg->B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    return 0;
}

