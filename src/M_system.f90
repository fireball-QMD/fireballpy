module M_system

  integer :: max_scf_iterations = 200

  integer :: natoms
  real, dimension (:, :), allocatable :: ratom
  integer, dimension (:), allocatable :: degelec
  integer, dimension (:), allocatable :: imass
  character (len = 2), dimension (:), allocatable:: symbol
  integer :: mbeta_max  
  integer neigh_max

  integer :: icluster


  ! --- ETOT ----
  real :: etot
  real :: etotold, etotnew
  real :: etotper
  real :: atomic_energy, efermi
  real :: uxcdcc_hf
  real :: uiiuee
  real :: uxcdcc
  real :: uxcdcc_ols
  real :: etotxc_1c
  real :: ebs
  real :: eqmmm
  real :: dxcv
  real :: dc_v_intra_dip_1c
  integer :: Kscf  

  ! -- EWALD--
  real, dimension (:, :), allocatable :: ewald
  real, dimension (:, :, :), allocatable :: dewald
  real, dimension (:, :), allocatable :: fewald

  ! --- NEIGH ----
  integer, dimension (:, :), allocatable :: neigh_b  
  integer, dimension (:, :), allocatable :: neigh_j
  integer, dimension (:), allocatable :: neighn
  integer, dimension (:, :, :), allocatable :: neigh_comb
  integer, dimension (:, :, :), allocatable :: neigh_comj
  integer, dimension (:, :, :), allocatable :: neigh_com_ng
  integer, dimension (:, :), allocatable :: neigh_comm
  integer, dimension (:), allocatable :: neigh_comn 
  integer, dimension (:,:), allocatable :: neigh_back
  integer, dimension (:), allocatable :: neigh_self
  integer, dimension (:, :), allocatable :: nPP_b
  integer, dimension (:, :), allocatable :: nPP_j
  integer, dimension (:, :), allocatable :: nPP_map
  integer, dimension (:), allocatable :: nPPn
  integer, dimension (:), allocatable :: nPP_self
  integer, dimension (:, :), allocatable :: nPPx_b
  integer, dimension (:, :), allocatable :: nPPx_j
  integer, dimension (:, :), allocatable :: nPPx_map
  integer, dimension (:, :), allocatable :: nPPx_point
  integer, dimension (:), allocatable :: nPPxn
  integer, dimension (:), allocatable :: nPPx_self
  integer, dimension(:), allocatable :: neigh_pair_a1
  integer, dimension(:), allocatable :: neigh_pair_a2
  integer, dimension(:), allocatable :: neigh_pair_n1
  integer, dimension(:), allocatable :: neigh_pair_n2
  integer, dimension(:,:), allocatable    :: neighj_tot
  integer, dimension(:,:), allocatable    :: neighb_tot
  integer, dimension(:), allocatable      :: neighn_tot


  !CHARGES
  real, dimension (:, :), allocatable :: Qin
  real, dimension (:), allocatable :: Qinmixer
  real, dimension (:, :), allocatable :: Qout
  real, dimension (:), allocatable :: Qoutmixer
  real, dimension (:), allocatable :: dq

  real, dimension (:), allocatable :: QLowdin_TOT
  real, dimension (:), allocatable :: QMulliken_TOT
  real, allocatable, dimension(:,:) :: qaux
  !Dipoles for iqout = 7
  real    dip_x, dipQout_x, dipTot_x, dipProy_x, dipIntra_x, dip_res_x, dipQin_x, dipRes_x
  real    dip_y, dipQout_y, dipTot_y, dipProy_y, dipIntra_y, dip_res_y, dipQin_y, dipRes_y
  real    dip_z, dipQout_z, dipTot_z, dipProy_z, dipIntra_z, dip_res_z, dipQin_z, dipRes_z
  real    dip_tot, dip_proy, dipQin_tot, dipTot_tot, dipIntra_tot, dipQout_tot, dip_res_tot, dipRes_tot 
  real, dimension (:), allocatable :: dq_DP


  ! anderson iteration procedure
  real, allocatable, dimension(:,:) :: Fv   ! x_try-x_old 
  real, allocatable, dimension(:,:) :: Xv   ! x_old 
  real, allocatable, dimension(:,:) :: delF! F(m+1)-F(m) 
  real, allocatable, dimension(:,:) :: delX! X(m+1)-X(m) 
  real, allocatable, dimension(:)   :: r2_sav
  ! Broyden mixing
  real, allocatable, dimension(:,:) :: RJac
  ! Louie mixing
  real, allocatable, dimension(:,:) :: betaInvH
  real, allocatable, dimension(:,:) :: gamaH

  real, dimension (:,:), allocatable :: xl

  ! Lattice vectors
  real, dimension (3) :: a1vec
  real, dimension (3) :: a2vec
  real, dimension (3) :: a3vec

  !esto es otra peli
  ! allocate_f
  real, dimension (:, :), allocatable :: dusr


end module
