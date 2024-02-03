module M_system

  !=======================
  integer :: icluster = 1
  integer :: iforce = 0
  integer :: idipole = 1
  integer :: iqout = 7
  !======================

  integer :: max_scf_iterations = 200
  real, parameter ::  xc_overtol = 5.0d-5
  real, parameter :: smt_elect = 0.8d0 ! Ewald and electrostatic

  integer :: natoms
  real, dimension (:, :), allocatable :: ratom
  integer, dimension (:), allocatable :: degelec
  integer, dimension (:), allocatable :: imass
  character (len = 2), dimension (:), allocatable:: symbol
  integer :: mbeta_max  
  integer neigh_max

  !--- mandar a M_fdata ? AQUI
  integer, dimension (:), allocatable :: getmssh
  integer, dimension (:), allocatable :: getlssh
  integer, dimension (:), allocatable :: getissh
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
  integer :: itestrange !AQUI pensar
  integer :: testrange !AQUI pensar
  real :: Uexc_1c
  real ::  Umuxc_1c
  real, dimension (:, :, :, :), allocatable :: vxc_1c

  ! -- EWALD--
  real, dimension (:, :), allocatable :: ewald
  real, dimension (:, :, :), allocatable :: dewald
  real, dimension (:, :), allocatable :: fewald
  real, dimension (:, :, :, :), allocatable :: ewaldsr
  real, dimension (:, :, :, :), allocatable :: dip
  real, dimension (:, :, :, :, :), allocatable :: dipp
 
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
  integer :: numorb_max
  integer, dimension (:), allocatable :: neighPP_self
  integer, dimension (:), allocatable :: neighPPn
  integer, dimension (:, :), allocatable :: neighPP_b
  integer, dimension (:, :), allocatable :: neighPP_j
  
  !CHARGES
  real, dimension (:, :), allocatable :: Qin
  real, dimension (:), allocatable :: Qinmixer
  real, dimension (:, :), allocatable :: Qout
  real, dimension (:), allocatable :: Qoutmixer
  real, dimension (:), allocatable :: dq

  real, dimension (:), allocatable :: QLowdin_TOT
  real, dimension (:), allocatable :: QMulliken_TOT
  real, allocatable, dimension(:,:) :: qaux
  real    dip_x, dipQout_x, dipTot_x, dipProy_x, dipIntra_x, dip_res_x, dipQin_x, dipRes_x
  real    dip_y, dipQout_y, dipTot_y, dipProy_y, dipIntra_y, dip_res_y, dipQin_y, dipRes_y
  real    dip_z, dipQout_z, dipTot_z, dipProy_z, dipIntra_z, dip_res_z, dipQin_z, dipRes_z
  real    dip_tot, dip_proy, dipQin_tot, dipTot_tot, dipIntra_tot, dipQout_tot, dip_res_tot, dipRes_tot 
  real, dimension (:), allocatable :: dq_DP

  !interaccions
  real, dimension (:, :, :, :), allocatable :: vxc
  real, dimension (:, :, :, :), allocatable :: vxc_ca
  real, dimension (:, :, :, :), allocatable :: rho_off
  real, dimension (:, :, :, :), allocatable :: rhoij_off
  real, dimension (:, :, :, :), allocatable :: s_mat 
  real, dimension (:, :, :, :), allocatable :: sm_mat
  real, dimension (:, :, :, :, :), allocatable :: spm_mat
  real, dimension (:, :, :), allocatable :: rho_on
  real, dimension (:, :, :), allocatable :: arho_on
  real, dimension (:, :, :), allocatable :: rhoi_on
  real, dimension (:, :, :), allocatable :: arhoi_on
  real, dimension (:, :, :, :, :), allocatable :: arhop_on
  real, dimension (:, :, :, :, :), allocatable :: rhop_on
  real, dimension (:, :, :, :), allocatable :: arhoij_off
  real, dimension (:, :, :, :), allocatable :: arho_off
  real, dimension (:, :, :, :, :), allocatable :: arhopij_off 
  real, dimension (:, :, :, :, :), allocatable :: arhop_off
  real, dimension (:, :, :, :, :), allocatable :: rhop_off
  real, dimension (:, :, :, :, :), allocatable :: rhopij_off
  real, dimension (:, :, :, :), allocatable :: vca
  real, dimension (:,:,:,:,:,:), allocatable :: gvhxc
  real, dimension (:, :, :, :), allocatable :: ewaldlr
  real, dimension (:, :, :, :), allocatable :: h_mat
  real, dimension (:, :, :, :), allocatable :: t_mat
  real, dimension (:, :, :, :), allocatable :: vna
  real, dimension (:, :, :, :), allocatable :: ewaldqmmm
  real, dimension (:,:,:), allocatable :: Vdip_1c
  real, dimension (:, :, :, :, :), allocatable :: dipc

  real, allocatable, dimension(:,:) :: Fv   
  real, allocatable, dimension(:,:) :: Xv    
  real, allocatable, dimension(:,:) :: delF 
  real, allocatable, dimension(:,:) :: delX 
  real, allocatable, dimension(:)   :: r2_sav
 
  real, allocatable, dimension(:,:) :: RJac
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
