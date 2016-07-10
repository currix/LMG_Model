MODULE defparam_Lipkin
  !
  USE nrtype
  !
  IMPLICIT NONE
  !
  ! Type definitions
  TYPE :: so2_bas
     INTEGER(KIND = I4B) :: N_U2_val
     INTEGER(KIND = I4B) :: m_SO2_val ! This label is really 2*m
  END TYPE so2_bas
  !
  TYPE :: u1_bas
     INTEGER(KIND = I4B) :: N_U2_val
     INTEGER(KIND = I4B) :: n_U1_val
  END TYPE u1_bas
  !
  TYPE :: exp_level
     REAL(KIND = DP) :: exp_term_value
     REAL(KIND = DP) :: exp_err_value
     INTEGER(KIND = I4B) :: exp_v
     CHARACTER(LEN = 3) :: reference
  END TYPE exp_level
  !
  TYPE term_values
     TYPE(exp_level) :: exper_level
     TYPE(term_values), POINTER :: next
  END TYPE term_values
  !
  !
  INTEGER(KIND = I4B) :: N_val ! U(2) [N] value
  !
  INTEGER(KIND = I4B) :: dim_block ! L Block dimension
  !
  ! Hamiltonian Matrix
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: Ham_mat
  ! Eigenvalue Matrix
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: Eigenval_vector ! Hamiltonian Eigenvalues
  ! Cylindrical Oscillator Basis
  TYPE(u1_bas), DIMENSION(:), ALLOCATABLE :: U1_Basis   !  U(2) > U(1) basis
  ! Displaced Oscillator Basis
  TYPE(so2_bas), DIMENSION(:), ALLOCATABLE :: SO2_Basis !  U(2) > SO(2) basis
  !
  ! Model Hamiltonian Parameters
  INTEGER(KIND = I4B), PARAMETER :: n_modham_param = 3
  REAL(KIND = DP), DIMENSION(1:n_modham_param) :: ModHam_parameter 
  !
  INTEGER(KIND = I4B) :: Ierr, Iprint
  !
  LOGICAL :: Eigenvec_Log     ! If .T. compute eigenvalues and eigenvectors
  LOGICAL :: Excitation_Log   ! If .T. compute excitation energy with respect to L = 0 g.s.
  LOGICAL :: Save_avec_Log   ! If .T. save eigenvector components.
  LOGICAL :: MXZ_Log   ! If .T. save magnetization
  !
  !
  REAL(KIND = DP) :: GS_energy = 0.0_DP ! L = 0 Ground state energy 
  !
  !
  REAL(KIND = DP), PARAMETER :: Zero_Parameter = 1.0E-20_DP ! Zero limit to discard a parameter evaluation
  !
  REAL(KIND = SP) :: time_check, time_check_ref ! Auxiliary parameters for benchmarking
  !
CONTAINS
  !
  SUBROUTINE U1_BASIS_VIBRON(N_val, U1_Basis) 
    !
    ! Subroutine to build the U(2) > U(1) basis in the Lipkin Model
    !
    !  |N 0>, |N 1>, |N 2>, ... , |N N> 
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val ! U(2) [N]
    !
    TYPE(u1_bas), DIMENSION(:), INTENT(OUT) :: U1_Basis
    !
    ! Local Variables
    INTEGER(KIND = I4B) :: index, nt = 0
    !
    index = 1_I4B
    !
    DO WHILE (nt <= N_val) 
       !
       U1_Basis(index)%n_U1_val = nt
       !
       index = index + 1_I4B
       nt = nt + 1_I4B
       !
    ENDDO
    !
    U1_Basis(1:Index-1)%N_U2_val = N_val
    !
  END SUBROUTINE U1_BASIS_VIBRON
  !
  !
  !
  SUBROUTINE SO2_BASIS_VIBRON(N_val, SO2_Basis) 
    !
    ! Subroutine to build the U(2) > SO(2) basis in the Lipkin model
    !
    !  |j -j>, |j -j+1>, |j -j+2>, ... , |j j> 
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val ! U(2) [N]
    !
    TYPE(so2_bas), DIMENSION(:), INTENT(OUT) :: SO2_Basis
    !
    ! Local Variables
    INTEGER(KIND = I4B) :: index, mval
    !
    index = 1_I4B
    mval = -N_val ! mval is 2*m
    !
    DO WHILE (mval <= N_val) 
       !
       SO2_Basis(index)%m_SO2_val = mval
       !
       index = index + 1_I4B
       mval = mval + 2_I4B
       !
    ENDDO
    !
    SO2_Basis(1:Index-1)%N_U2_val = N_val
    !
  END SUBROUTINE SO2_BASIS_VIBRON
  !
  !
  !
  SUBROUTINE Build_Mod_Ham_U1(N_val, dim_block, Ham_mat) 
    !
    !
    ! Subroutine to build the U(2) Lipkin Model Hamiltonian
    !
    ! U2-U1 basis
    ! Lipkin Model Hamiltonian:
    ! H_QQ = (1-xi) n_t  - (xi/N) Q^y Q^y + delta Q^y
    ! n_t = t^+t
    ! Q^y = t^+s + s^+t + y t^+t
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val ! U(2) [N]
    !
    INTEGER(KIND = I4B), INTENT(IN) :: dim_block ! Angular momentum L_val block dimension
    !
    !
    REAL(KIND = DP), DIMENSION(:,:), INTENT(OUT) :: Ham_mat ! Hamiltonian tridiagonal matrix
    !
    REAL(KIND = DP), DIMENSION(1:dim_block) :: diag, nt_values
    INTEGER(KIND = I4B) :: ntvalue, index
    REAL(KIND = DP) :: Nvalue
    !
    !
    Nvalue = REAL(N_val, DP)
    !
    !
    ! Build Hamiltonian
    !
    DO ntvalue = 0, N_val
       nt_values(ntvalue+1) = REAL(ntvalue, DP)
    ENDDO
    !
    ! Main Diagonal
    !
    diag = (1.0_DP-ModHam_parameter(1))*nt_values &
         - ModHam_parameter(1)*((2.0_DP*nt_values+1.0_DP) - (2.0_DP/Nvalue)*nt_values**2) &
         - (ModHam_parameter(1)/Nvalue)*(ModHam_parameter(2)**2)*nt_values**2 
    !
    DO index = 1, dim_block
       Ham_mat(index, index) = diag(index)
    ENDDO
    !
    !
    ! First Diagonal
    !
    diag = 0.0_DP
    diag(1:dim_block-1) = -(ModHam_parameter(1)/Nvalue)*ModHam_parameter(2)* &
         (2.0_DP*nt_values(1:dim_block-1) + 1.0_DP)*SQRT((Nvalue-nt_values(1:dim_block-1))*(nt_values(1:dim_block-1)+1.0_DP))+ &
         ModHam_parameter(3)*SQRT((Nvalue-nt_values(1:dim_block))*(nt_values(1:dim_block-1)+1.0_DP))
    !
    DO index = 1, dim_block - 1
       Ham_mat(index, index + 1) = diag(index)
    ENDDO
    !
    !
    ! Second Diagonal
    !
    diag = 0.0_DP
    diag(1:dim_block-2) = -(ModHam_parameter(1)/Nvalue)* &
         SQRT((Nvalue-nt_values(1:dim_block-2))*(Nvalue-nt_values(1:dim_block-2)-1.0_DP)* &
         (nt_values(1:dim_block-2)+1.0_DP)*(nt_values(1:dim_block-2)+2.0_DP))
    !
    DO index = 1, dim_block - 2
       Ham_mat(index, index + 2) = diag(index)
    ENDDO
    !
    RETURN
    !
    !
  END SUBROUTINE BUILD_MOD_HAM_U1
  !
  !
  SUBROUTINE Build_Mod_Ham_SO2(N_val, dim_block, SO2_Basis, Ham_mat) 
    !
    !
    ! Subroutine to build the U(2) Lipkin Model Hamiltonian
    !
    ! U2-SO2 basis
    ! Lipkin Model Hamiltonian:
    ! H_QQ = (1-xi) n_t  - (xi/N) Q^0 Q^0 + delta Q^0
    ! n_t = t^+t
    ! Q^0 = t^+s + s^+t 
    ! 
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val ! U(2) [N]
    !
    INTEGER(KIND = I4B), INTENT(IN) :: dim_block ! Angular momentum L_val block dimension
    !
    TYPE(so2_bas), DIMENSION(:), INTENT(IN) :: SO2_Basis
    !
    !
    REAL(KIND = DP), DIMENSION(:,:), INTENT(OUT) :: Ham_mat ! Hamiltonian tridiagonal matrix
    !
    REAL(KIND = DP), DIMENSION(1:dim_block) :: diag, m2_values
    INTEGER(KIND = I4B) :: index
    REAL(KIND = DP) :: Nvalue
    !
    !
    Nvalue = REAL(N_val, DP)
    !
    !
    ! Build Hamiltonian
    !
    ! Extract m values from the basis 2m values
    DO index = 1, dim_block
       m2_values(index) = REAL(SO2_Basis(index)%m_SO2_val, DP)
    ENDDO
    !
    ! Main Diagonal
    !
    diag = (1.0_DP-ModHam_parameter(1))*Nvalue/2.0_DP &
         - (ModHam_parameter(1)/Nvalue)*(m2_values**2) + &
         ModHam_parameter(3)*m2_values
    !
    DO index = 1, dim_block
       Ham_mat(index, index) = diag(index)
    ENDDO
    !
    !
    ! First Diagonal m' = m + 1
    !
    ! <j m + 1| H | j m >
    !
    diag = 0.0_DP
    diag(1:dim_block-1) = -((1.0_DP-ModHam_parameter(1))/(4.0_DP))* &
         SQRT((Nvalue - m2_values(1:dim_block - 1))*(Nvalue + m2_values(1:dim_block-1) + 2.0_DP))
    !
    DO index = 1, dim_block - 1
       Ham_mat(index, index + 1) = diag(index)
    ENDDO
    !
    RETURN
    !
    !
  END SUBROUTINE BUILD_MOD_HAM_SO2
  !
  !
  SUBROUTINE SAVE_EIGENV_COMPONENTS(N_val, dim_block, Basis, Ham_mat)
    !
    ! Subroutine to save the components of the eigenvectors 
    ! of the Lipkin Model Hamiltonian
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val ! U(2) [N]
    !
    INTEGER(KIND = I4B), INTENT(IN) :: dim_block ! Angular momentum L_val block dimension
    !
    CHARACTER(LEN=*), INTENT(IN) :: Basis
    !
    REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) :: Ham_mat ! Hamiltonian matrix
    !
    INTEGER(KIND = I4B) :: basis_index
    !
    CHARACTER(LEN=65) :: output_filename
    !
    !
    ! Build filename
    IF ( N_val < 10) THEN !to avoid spaces
       WRITE(output_filename, '("eigvec_",A,"_N",I1,".dat")') TRIM(Basis), N_val
    ELSE IF ( N_val < 100) THEN 
       WRITE(output_filename, '("eigvec_",A,"_N",I2,".dat")') TRIM(Basis), N_val
    ELSE IF ( N_val < 1000) THEN 
       WRITE(output_filename, '("eigvec_",A,"_N",I3,".dat")') TRIM(Basis), N_val
    ELSE IF ( N_val < 10000) THEN 
       WRITE(output_filename, '("eigvec_",A,"_N",I4,".dat")') TRIM(Basis), N_val
    ELSE IF ( N_val < 100000) THEN 
       WRITE(output_filename, '("eigvec_",A,"_N",I5,".dat")') TRIM(Basis), N_val
    ELSE
       WRITE(output_filename, '("eigvec_",A,"_N",I6,".dat")') TRIM(Basis), N_val
    ENDIF
    !
    OPEN(UNIT = 76, FILE = output_filename, STATUS = "UNKNOWN", ACTION = "WRITE")
    !
    WRITE(UNIT=76, FMT=*) "# N = ", N_val, " ; ", Basis,  " basis,  xi = ",  ModHam_parameter(1), "; y = ", ModHam_parameter(2)
    !
    DO basis_index = 1, dim_block
       WRITE(UNIT=76, FMT=*) Ham_mat(basis_index, 1:dim_block)
    ENDDO
    !
  END SUBROUTINE SAVE_EIGENV_COMPONENTS
  !
  SUBROUTINE magnetization_XZ_U1(N_val, dim, Eigenvector, Magx, Magz)
    !
    ! Subroutine to compute <Jx> and <Jz> for eigenvectors 
    ! of the Lipkin Model Hamiltonian expressed in the U(1) basis
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val ! U(2) [N]
    !
    INTEGER(KIND = I4B), INTENT(IN) :: dim ! Block dimension
    !
    REAL(KIND = DP), DIMENSION(:), INTENT(IN) :: Eigenvector ! Eigenvector components
    REAL(KIND = DP), INTENT(OUT) :: Magx, Magz ! Magnetization
    !
    ! We assume a basis |j mx> = |j -j>, |j -j + 1>, ..., |j j>
    INTEGER(KIND = I4B) :: basis_index, i, mx
    !
    ! Compute <psi|Mx|psi>
    Magx = 0.0_DP
    Magx = SUM( Eigenvector**2 *  (/ (i, i = N_val/2, -N_val/2, -1) /) )
    ! Compute <psi|Mz|psi>
    Magz = 0.0_DP
    DO basis_index = 2, dim ! mx = j-1 ... -j
       mx = N_val/2 - basis_index + 1
       Magz = Magz + Eigenvector(basis_index)*Eigenvector(basis_index-1)*JX(N_val/2, mx, mx+1)
    ENDDO
    DO basis_index = 1, dim-1 ! mx = j ... -j + 1
       mx = N_val/2 - basis_index + 1
       Magz = Magz + Eigenvector(basis_index)*Eigenvector(basis_index+1)*JX(N_val/2, mx, mx-1)
    ENDDO
  END SUBROUTINE magnetization_XZ_U1
  FUNCTION JX(J, M, MP)
    INTEGER(KIND = I4B), INTENT(IN) :: J, M, MP
    REAL(KIND = DP) :: JX
    !
    IF (MP == M + 1) THEN
       JX = 0.5_DP*SQRT(1.0_DP*(J + M + 1)*(J - M))
       ELSE IF (MP == M - 1) THEN
       JX = 0.5_DP*SQRT(1.0_DP*(J + M)*(J - M + 1))
    ELSE
       JX = 0.0_DP
    ENDIF
    !
  END FUNCTION JX
  !
  !
  SUBROUTINE magnetization_XZ_SO2(N_val, dim, Eigenvector, Magx, Magz)
    !
    ! Subroutine to compute <Jx> and <Jz> for eigenvectors 
    ! of the Lipkin Model Hamiltonian expressed in the SO(2) basis
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val ! U(2) [N]
    !
    INTEGER(KIND = I4B), INTENT(IN) :: dim ! Block dimension
    !
    REAL(KIND = DP), DIMENSION(:), INTENT(IN) :: Eigenvector ! Eigenvector components
    REAL(KIND = DP), INTENT(OUT) :: Magx, Magz ! Magnetization
    !
    ! We assume a basis |j mz> = |j -j>, |j -j + 1>, ..., 
    INTEGER(KIND = I4B) :: basis_index, mz
    !
    ! Compute <psi|Mz|psi>
    Magz = 0.0_DP
    Magz = SUM( Eigenvector**2 *  (/ (mz, mz = -N_val/2, N_val/2) /) )
    ! Compute <psi|Mx|psi>
    Magx = 0.0_DP
    DO basis_index = 2, dim ! mz = -j+1 ... j
       mz = -N_val/2 + basis_index - 1
       Magx = Magx + Eigenvector(basis_index)*Eigenvector(basis_index-1)*JX(N_val/2, mz, mz-1)
    ENDDO
    DO basis_index = 1, dim-1 ! mz = -j ... j - 1
       mz = -N_val/2 + basis_index - 1
       Magx = Magx + Eigenvector(basis_index)*Eigenvector(basis_index+1)*JX(N_val/2, mz, mz+1)
    ENDDO
  END SUBROUTINE magnetization_XZ_SO2
  !
END MODULE defparam_Lipkin
