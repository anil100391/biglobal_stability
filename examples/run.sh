#mpirun -np 4 main_ose -eps_nev 40 -eps_ncv 84 -st_type sinvert -eps_target 2
#mpirun -np 4 coflow_wakes -eps_nev 5 -eps_ncv 12 -st_type sinvert -eps_target 3 -st_ksp_type preonly -st_pc_type lu -st_pc_factor_mat_solver_package mumps -mat_mumps_icntl_14 80
mpirun -np 4 coflow_wakes -eps_nev 5 -eps_ncv 12 -st_type sinvert -eps_target 3
