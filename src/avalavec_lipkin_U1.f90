PROGRAM avalavec_modelH_LIPKIN
  !
  ! Program to compute eigenvalues and eigenvectors of the Lipkin 
  ! Model Hamiltonian
  !
  ! Basis U2-U1
  !
  ! by Currix TM.
  !
  !
  USE nrtype
  !
  USE defparam_Lipkin
  !
  ! Lapack 95
  USE LA_PRECISION, ONLY: WP => DP
  USE F95_LAPACK, ONLY: LA_SYEVR
  !
  !
  IMPLICIT NONE
  !
  REAL(KIND = DP) :: xi_par, y_par, Delta ! Model Hamiltonian Parameters
  REAL(KIND = DP) :: Mx, Mz        ! Magnetization expected values
  !
  INTEGER(KIND = I4B) :: state_index, state_index_2, nt
  !
  !
  ! NAMELISTS
  NAMELIST/par_aux/ Iprint, Eigenvec_Log, Excitation_Log, Save_avec_Log
  NAMELIST/par_0/ N_val
  NAMELIST/par_1/ xi_par, y_par, Delta
  NAMELIST/par_2/ MXZ_Log
  !
  ! 
  ! Initialize time
  CALL CPU_TIME(time_check_ref)
  !
  !
  ! Read parameters
  !
  READ(UNIT = *, NML = par_aux)
  !
  IF (Iprint > 1) PRINT 10
  READ(UNIT = *, NML = par_0)
  !
  IF (Iprint > 1) PRINT 20
  READ(UNIT = *, NML = par_1)
  !
  IF (Iprint > 1) PRINT 30
  READ(UNIT = *, NML = par_2)
  !
  !
  IF (Iprint > 1) THEN
     WRITE(UNIT = *, FMT = 5) Iprint, Eigenvec_Log, Excitation_Log, Save_avec_Log
     WRITE(UNIT = *, FMT = 15) N_val
     WRITE(UNIT = *, FMT = 25) xi_par, y_par, delta
  ENDIF
  !
  ! HAMILTONIAN PARAMETERS
  !modham_param(1) => nt
  ModHam_parameter(1) = xi_par
  !modham_param(2) => Q0 x Q0
  ModHam_parameter(2) = y_par
  !modham_param(3) => Q0 (Perturbation to break degeneracy)
  ModHam_parameter(3) = delta
  !
  !
  ! BLOCK DIMENSION
  dim_block = N_val + 1
  !
  ! ALLOCATE BASIS
  ALLOCATE(U1_Basis(1:dim_block), STAT = IERR)    
  IF (IERR /= 0) THEN
     WRITE(UNIT = *, FMT = *) "U1_Basis allocation request denied."
     STOP
  ENDIF
  !
  CALL U1_BASIS_VIBRON(N_val, U1_Basis) ! Build U1 basis
  !
  ! Build Hamiltonian Matrix
  !
  ! ALLOCATE Hamiltonian Matrix 
  ALLOCATE(Ham_mat(1:dim_block,1:dim_block), STAT = IERR)    
  IF (IERR /= 0) THEN
     WRITE(UNIT = *, FMT = *) "Ham_mat allocation request denied."
     STOP
  ENDIF
  !
  !
  Ham_mat = 0.0_DP
  CALL Build_Mod_Ham_U1(N_val, dim_block, Ham_mat) 
  !
  !
  IF (Iprint > 2) THEN
     WRITE(*,*) ' '
     WRITE(*,*) ' Hamiltonian Matrix'
     DO state_index = 1, dim_block                         
        WRITE(*,*) (Ham_mat(state_index, state_index_2), state_index_2 = 1, dim_block)
     ENDDO
  ENDIF
  !
  ! Hamiltonian Diagonalization
  !
  !
  ! ALLOCATE EIGENVALUES VECTOR
  ALLOCATE(Eigenval_vector(1:dim_block), STAT = IERR)    
  IF (IERR /= 0) THEN
     WRITE(UNIT = *, FMT = *) "Eigenval_vector allocation request denied."
     STOP
  ENDIF
  !
  !      
  !
  ! Diagonalize Hamiltonian matrix (LAPACK95)
  IF (Eigenvec_Log .OR. Save_avec_Log .OR. MXZ_Log) THEN
     CALL LA_SYEVR(A=Ham_mat, W=Eigenval_vector, JOBZ='V', UPLO='U')
  ELSE
     CALL LA_SYEVR(A=Ham_mat, W=Eigenval_vector, JOBZ='N', UPLO='U')
  ENDIF
  !
  IF (Excitation_Log) THEN
     !
     GS_energy = Eigenval_vector(1)
     eigenval_vector = Eigenval_vector - GS_energy
     !
  ENDIF
  !
  !
  DO state_index = 1, dim_block
     !
     nt = U1_Basis(state_index)%n_U1_val
     !
     IF (MXZ_Log) THEN
        CALL magnetization_XZ_U1(N_val, dim_block, Ham_mat(:, state_index), MX, MZ)
        WRITE(UNIT = *, FMT = *) nt, Eigenval_vector(state_index)/REAL(N_val,DP), MX/N_val, MZ/N_val
     ELSE
        WRITE(UNIT = *, FMT = *) nt, Eigenval_vector(state_index)/REAL(N_val,DP)
     ENDIF
     !
     !
     IF (Eigenvec_Log .AND. Iprint > 0) THEN
        !
        ! Display eigenvectors
        DO state_index_2 = 1, dim_block
           !
           nt = U1_Basis(state_index_2)%n_U1_val
           !
           WRITE(UNIT = *, FMT = *) Ham_mat(state_index_2, state_index), "|",nt,">"
           !
        ENDDO
        !
     ENDIF
     !
  ENDDO
  !
  ! Save eigenvector components
  IF (Save_avec_Log) CALL SAVE_EIGENV_COMPONENTS(N_val, dim_block, "u1", Ham_mat)
  !
  !    
  ! DEALLOCATE EIGENVALUES VECTOR
  DEALLOCATE(Eigenval_vector, STAT = IERR)    
  IF (IERR /= 0) THEN
     WRITE(UNIT = *, FMT = *) "Eigenval_vector deallocation request denied."
     STOP
  ENDIF
  ! DEALLOCATE Hamiltonian Matrix
  DEALLOCATE(Ham_mat, STAT = IERR)    
  IF (IERR /= 0) THEN
     WRITE(UNIT = *, FMT = *) "Ham_mat allocation request denied."
     STOP
  ENDIF
  ! DEALLOCATE BASIS
  DEALLOCATE(U1_Basis, STAT = IERR)    
  IF (IERR /= 0) THEN
     WRITE(UNIT = *, FMT = *) "U1_Basis deallocation request denied."
     STOP
  ENDIF
  !
5 FORMAT(1X, " Iprint = ", I2, "; Eigenvec_LOG = ", L2, "; Excitation_Log = ", L2, "; Save_Avec_Log = ", L2)
10 FORMAT(1X, "Reading  N_val")
15 FORMAT(1X, "N_val = ", I6)
20 FORMAT(1X, "Reading xi, y")
25 FORMAT(1X, "xi = ", ES14.7, "; y = ", ES14.7, "; Delta = ", ES14.7)
30 FORMAT(1X, "Reading MXZ_Log")
  !
  !
END PROGRAM avalavec_modelH_LIPKIN
