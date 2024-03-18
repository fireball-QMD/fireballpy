module M_system

  !=======================
  integer :: iforce   = 1
  integer :: idipole  = 1 !AQUI  creo que para 0 no funciona arreglar
  integer :: iqout    = 1 !1:Lowdin 2:Mulliken 3:NPA 4:M-dipole :7MD-pres..
  !ifixcharge funcione
  integer :: icluster = 1 !el 0 tiene que funcionar
  ! si icluster = 0 idipole = 0
  !k_points=gamma AQUI pensar
  ! pensar gamma ?
  logical :: iqmmm = .false.
  !======================

  integer :: ialgmix  = 1 !1:anderson 2:broyden 3:louie 4:pulay
  real, parameter :: xc_overtol = 5.0d-5
  real, parameter :: smt_elect = 0.8d0 ! Ewald and electrostatic
  integer, parameter :: ithetamax = 5
  integer, parameter :: idmix = 6
  integer, parameter :: max_scf_iterations = 200
  real :: tempfe = 100.0d0
  real :: bmix = 0.04d0
  real :: sigma = 0.0d0
  real :: sigmaold = 0.0d0
  real :: sigmatol = 1.0E-8
  logical ::  scf_achieved = .false.
 
  integer :: natoms
  real, dimension (:, :), allocatable :: ratom
  integer, dimension (:), allocatable :: degelec
  integer, dimension (:), allocatable :: imass
  character (len = 2), dimension (:), allocatable:: symbol
  integer :: mbeta_max  
  integer :: neigh_max
  integer :: ishiftO
  real, dimension (3) :: shifter

  !--- mandar a M_fdata ? AQUI
  integer, dimension (:), allocatable :: getmssh
  integer, dimension (:), allocatable :: getlssh
  integer, dimension (:), allocatable :: getissh

  ! --- ETOT ----
  real :: etot
  real :: etotold, etotnew
  real :: etotper
  real :: atomic_energy, efermi
  real :: uiiuee
  real :: uxcdcc
  real :: uxcdcc_ols
  real :: etotxc_1c
  real :: ebs
  real :: eqmmm
  real :: dc_v_intra_dip_1c
  integer :: Kscf  
  real :: Uexc_1c
  real ::  Umuxc_1c
  real, dimension (:, :, :, :), allocatable :: vxc_1c

  !Charges
  real, dimension(:), allocatable  :: Q0_TOT
  integer, dimension (:), allocatable :: nelectron

  !--diag--
  real, dimension (:, :, :), allocatable :: blowre
  real, dimension (:, :, :), allocatable :: bbnkre
  real, dimension (:, :, :), allocatable :: blowim
  real, dimension (:, :, :), allocatable :: bbnkim
  real, dimension (:, :), allocatable :: eigen_k
  real, dimension (:, :), allocatable :: special_k
  integer :: norbitals
  integer :: norbitals_new
  integer :: nkpoints  
  integer, dimension (:), allocatable :: getiatom
  integer, dimension (:,:),allocatable :: ioccupy_k !AQUI pensar allo en denmat
  integer, dimension (:), allocatable :: ioccupy   !AQUI
  real, dimension (:,:), allocatable :: foccupy !AQUI
  
  !scf
  real, dimension (:, :, :, :), allocatable :: cape
  real, dimension (:, :, :), allocatable :: blowre_o
  real, dimension (:, :, :), allocatable :: bbnkre_o  
  real, dimension (:, :, :, :), allocatable :: rhoPP
  real :: ztot
  real, dimension (:), allocatable :: weight_k
  integer ::  nssh_tot
  real, dimension (:), allocatable :: mwe
  real, dimension (:), allocatable :: drwe

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
  integer :: tot_pairs
  real, dimension (:, :, :, :), allocatable :: sVNL
  real, dimension (:, :, :, :, :), allocatable :: spVNL
  real, dimension (:, :, :, :, :), allocatable :: sp_mat
  real, dimension (:, :, :, :, :), allocatable :: tp_mat
  real, dimension (:, :, :), allocatable :: dipcm
  real, dimension (:, :, :, :), allocatable :: dippcm
  real, dimension (:, :, :, :, :, :), allocatable :: dippc
  real, dimension(3,3) :: eps2  
  real, dimension (:, :, :, :), allocatable :: vnl
  integer, dimension (:), allocatable :: neighPP_comn
  integer, dimension (:, :), allocatable :: neighPP_comm
  integer, dimension (:, :, :), allocatable :: neighPP_comj
  integer, dimension (:, :, :), allocatable :: neighPP_comb 
  integer,allocatable   :: neighj_aux(:,:)
  integer :: neighPP_max 
  integer :: num_neig_maxtot
 
  !CHARGES
  real, dimension (:, :), allocatable :: Qin
  real, dimension (:), allocatable :: Qinmixer
  real, dimension (:, :), allocatable :: Qout
  real, dimension (:), allocatable :: Qoutmixer
  real, dimension (:), allocatable :: dq
  real, dimension(:), allocatable  :: Q_partial
  real, dimension (:, :), allocatable :: Qin_es
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
  real, dimension (:, :, :, :), allocatable :: rho
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
  real, dimension (:, :, :, :), allocatable :: ewaldlr
  real, dimension (:, :, :, :), allocatable :: h_mat
  real, dimension (:, :, :, :), allocatable :: t_mat
  real, dimension (:, :, :, :), allocatable :: vna
  real, dimension (:, :, :, :), allocatable :: ewaldqmmm
  real, dimension (:, :, :, :, :), allocatable :: dipc
  integer, dimension(:,:), allocatable :: muR
  integer, dimension(:,:), allocatable :: nuR
  integer, dimension(:,:), allocatable :: alphaR
  integer, dimension(:,:), allocatable :: betaR
  real, dimension(:,:,:,:), allocatable   :: hr_box

  
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

  ! allocate_f
  real, dimension (:, :, :), allocatable :: fotnl
  real, dimension (:, :, :), allocatable :: fanl
  real, dimension (:, :, :), allocatable :: fotna
  real, dimension (:, :, :), allocatable :: fana
  real, dimension (:, :, :), allocatable :: faxc 
  real, dimension (:, :, :), allocatable :: faxc_ca 
  real, dimension (:, :, :), allocatable :: dxcdcc
  real, dimension (:, :), allocatable :: ft
  real, dimension (:, :), allocatable :: dusr
  real, dimension (:, :), allocatable :: special_k_orig
  real, dimension (:), allocatable :: weight_k_orig
  real, dimension (:, :), allocatable :: scale_k
  real, dimension (:, :, :), allocatable :: fotxc 
  real, dimension (:, :, :), allocatable :: fotxc_ca
  real, dimension (:, :, :), allocatable :: faca 
  real, dimension (:, :, :), allocatable :: fotca
  real, dimension (:, :), allocatable :: f3naa, f3nab, f3nac
  real, dimension (:, :), allocatable :: f3nla, f3nlb, f3nlc
  real, dimension (:, :), allocatable :: f3caa, f3cab, f3cac
  real, dimension (:, :), allocatable :: flrew
  real, dimension (:, :), allocatable :: f3xca_ca, f3xcb_ca, f3xcc_ca
  real, dimension (:, :), allocatable :: f3xca, f3xcb, f3xcc
  real, dimension (:, :), allocatable :: flrew_qmmm
  real, dimension (:, :), allocatable :: fro
  real, dimension (:, :), allocatable :: ftot
  real, dimension (:, :), allocatable :: dxcv

end module
