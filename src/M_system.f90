module M_system

  !========================
  integer :: iforce     = 1
  integer :: iqout      = 2 ! 1:Lowdin 2:Mulliken 3:NPA 4:M-dipole :7MD-pres..
  integer :: icluster   = 1 ! 
  integer :: idipole    = 1 
  integer :: igamma     = 1
  integer :: iqmmm      = 0
  integer :: ifixcharge = 0
  !========================
  !idipole icluster igamma 
  !1       1        1        !non periodic + only gamma kpts
  !0       1        1       
  !0       0        0        !periodic + kpoints
  !0       0        1        !periodic + only gamma kpts

  integer :: ialgmix  = 1 !1:anderson 2:broyden 3:louie 4:pulay
  real*8, parameter :: xc_overtol = 5.0d-5
  real*8, parameter :: smt_elect = 0.8d0 ! Ewald and electrostatic
  integer, parameter :: ithetamax = 5
  integer, parameter :: idmix = 6
  integer, parameter :: max_scf_iterations = 200
  real*8 :: tempfe = 100.0d0
  real*8 :: bmix = 0.04d0
  real*8 :: sigma = 0.0d0
  real*8 :: sigmaold = 0.0d0
  real*8 :: sigmatol = 1.0E-8
  logical ::  scf_achieved = .false.
 
  integer :: natoms
  real*8, dimension (:, :), allocatable :: ratom
  integer, dimension (:), allocatable :: degelec
  integer, dimension (:), allocatable :: imass
  character (len = 2), dimension (:), allocatable:: symbol
  integer :: mbeta_max  
  integer :: neigh_max
  integer :: ishiftO
  real*8, dimension (3) :: shifter

  !--- mandar a M_fdata ? AQUI
  integer, dimension (:), allocatable :: getmssh
  integer, dimension (:), allocatable :: getlssh
  integer, dimension (:), allocatable :: getissh

  ! --- ETOT ----
  real*8 :: etot
  real*8 :: etotold, etotnew
  real*8 :: etotper
  real*8 :: atomic_energy, efermi
  real*8 :: uiiuee
  real*8 :: uxcdcc
  real*8 :: uxcdcc_ols
  real*8 :: etotxc_1c
  real*8 :: ebs
  real*8 :: eqmmm
  real*8 :: dc_v_intra_dip_1c
  integer :: Kscf  
  real*8 :: Uexc_1c
  real*8 ::  Umuxc_1c
  real*8, dimension (:, :, :, :), allocatable :: vxc_1c

  !Charges
  real*8, dimension(:), allocatable  :: Q0_TOT
  integer, dimension (:), allocatable :: nelectron

  !--diag--
  real*8, dimension (:, :, :), allocatable :: blowre
  real*8, dimension (:, :, :), allocatable :: bbnkre
  real*8, dimension (:, :, :), allocatable :: blowim
  real*8, dimension (:, :, :), allocatable :: bbnkim
  real*8, dimension (:, :), allocatable :: eigen_k
  real*8, dimension (:, :), allocatable :: special_k
  integer :: norbitals
  integer :: norbitals_new
  integer :: nkpoints  
  integer, dimension (:), allocatable :: getiatom
  integer, dimension (:,:),allocatable :: ioccupy_k !AQUI pensar allo en denmat
  integer, dimension (:), allocatable :: ioccupy   !AQUI
  real*8, dimension (:,:), allocatable :: foccupy !AQUI
  
  !scf
  real*8, dimension (:, :, :, :), allocatable :: cape
  real*8, dimension (:, :, :), allocatable :: blowre_o
  real*8, dimension (:, :, :), allocatable :: bbnkre_o  
  real*8, dimension (:, :, :, :), allocatable :: rhoPP
  real*8 :: ztot
  real*8, dimension (:), allocatable :: weight_k
  integer ::  nssh_tot
  real*8, dimension (:), allocatable :: mwe
  real*8, dimension (:), allocatable :: drwe

  ! -- EWALD--
  real*8, dimension (:, :), allocatable :: ewald
  real*8, dimension (:, :, :), allocatable :: dewald
  real*8, dimension (:, :), allocatable :: fewald
  real*8, dimension (:, :, :, :), allocatable :: ewaldsr
  real*8, dimension (:, :, :, :), allocatable :: dip
  real*8, dimension (:, :, :, :, :), allocatable :: dipp
 
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
  real*8, dimension (:, :, :, :), allocatable :: sVNL
  real*8, dimension (:, :, :, :, :), allocatable :: spVNL
  real*8, dimension (:, :, :, :, :), allocatable :: sp_mat
  real*8, dimension (:, :, :, :, :), allocatable :: tp_mat
  real*8, dimension (:, :, :), allocatable :: dipcm
  real*8, dimension (:, :, :, :), allocatable :: dippcm
  real*8, dimension (:, :, :, :, :, :), allocatable :: dippc
  real*8, dimension(3,3) :: eps2  
  real*8, dimension (:, :, :, :), allocatable :: vnl
  integer, dimension (:), allocatable :: neighPP_comn
  integer, dimension (:, :), allocatable :: neighPP_comm
  integer, dimension (:, :, :), allocatable :: neighPP_comj
  integer, dimension (:, :, :), allocatable :: neighPP_comb 
  integer,allocatable   :: neighj_aux(:,:)
  integer :: neighPP_max 
  integer :: num_neig_maxtot
 
  !CHARGES
  real*8, dimension (:, :), allocatable :: Qin
  real*8, dimension (:), allocatable :: Qinmixer
  real*8, dimension (:, :), allocatable :: Qout
  real*8, dimension (:), allocatable :: Qoutmixer
  real*8, dimension (:), allocatable :: dq
  real*8, dimension(:), allocatable  :: Q_partial
  real*8, dimension (:, :), allocatable :: Qin_es
  real*8, dimension (:), allocatable :: QLowdin_TOT
  real*8, dimension (:), allocatable :: QMulliken_TOT
  real*8, allocatable, dimension(:,:) :: qaux
  real*8    dip_x, dipQout_x, dipTot_x, dipProy_x, dipIntra_x, dip_res_x, dipQin_x, dipRes_x
  real*8    dip_y, dipQout_y, dipTot_y, dipProy_y, dipIntra_y, dip_res_y, dipQin_y, dipRes_y
  real*8    dip_z, dipQout_z, dipTot_z, dipProy_z, dipIntra_z, dip_res_z, dipQin_z, dipRes_z
  real*8    dip_tot, dip_proy, dipQin_tot, dipTot_tot, dipIntra_tot, dipQout_tot, dip_res_tot, dipRes_tot 
  real*8, dimension (:), allocatable :: dq_DP

  !interaccions
  real*8, dimension (:, :, :, :), allocatable :: vxc
  real*8, dimension (:, :, :, :), allocatable :: vxc_ca
  real*8, dimension (:, :, :, :), allocatable :: rho
  real*8, dimension (:, :, :, :), allocatable :: rho_off
  real*8, dimension (:, :, :, :), allocatable :: rhoij_off
  real*8, dimension (:, :, :, :), allocatable :: s_mat 
  real*8, dimension (:, :, :, :), allocatable :: sm_mat
  real*8, dimension (:, :, :, :, :), allocatable :: spm_mat
  real*8, dimension (:, :, :), allocatable :: rho_on
  real*8, dimension (:, :, :), allocatable :: arho_on
  real*8, dimension (:, :, :), allocatable :: rhoi_on
  real*8, dimension (:, :, :), allocatable :: arhoi_on
  real*8, dimension (:, :, :, :, :), allocatable :: arhop_on
  real*8, dimension (:, :, :, :, :), allocatable :: rhop_on
  real*8, dimension (:, :, :, :), allocatable :: arhoij_off
  real*8, dimension (:, :, :, :), allocatable :: arho_off
  real*8, dimension (:, :, :, :, :), allocatable :: arhopij_off 
  real*8, dimension (:, :, :, :, :), allocatable :: arhop_off
  real*8, dimension (:, :, :, :, :), allocatable :: rhop_off
  real*8, dimension (:, :, :, :, :), allocatable :: rhopij_off
  real*8, dimension (:, :, :, :), allocatable :: vca
  real*8, dimension (:, :, :, :), allocatable :: ewaldlr
  real*8, dimension (:, :, :, :), allocatable :: h_mat
  real*8, dimension (:, :, :, :), allocatable :: t_mat
  real*8, dimension (:, :, :, :), allocatable :: vna
  real*8, dimension (:, :, :, :), allocatable :: ewaldqmmm
  real*8, dimension (:, :, :, :, :), allocatable :: dipc
  integer, dimension(:,:), allocatable :: muR
  integer, dimension(:,:), allocatable :: nuR
  integer, dimension(:,:), allocatable :: alphaR
  integer, dimension(:,:), allocatable :: betaR
  real*8, dimension(:,:,:,:), allocatable   :: hr_box

  
  real*8, allocatable, dimension(:,:) :: Fv   
  real*8, allocatable, dimension(:,:) :: Xv    
  real*8, allocatable, dimension(:,:) :: delF 
  real*8, allocatable, dimension(:,:) :: delX 
  real*8, allocatable, dimension(:)   :: r2_sav
 
  real*8, allocatable, dimension(:,:) :: RJac
  real*8, allocatable, dimension(:,:) :: betaInvH
  real*8, allocatable, dimension(:,:) :: gamaH
  real*8, dimension (:,:), allocatable :: xl

  ! Lattice vectors
  real*8, dimension (3) :: a1vec
  real*8, dimension (3) :: a2vec
  real*8, dimension (3) :: a3vec

  ! allocate_f
  real*8, dimension (:, :, :), allocatable :: fotnl
  real*8, dimension (:, :, :), allocatable :: fanl
  real*8, dimension (:, :, :), allocatable :: fotna
  real*8, dimension (:, :, :), allocatable :: fana
  real*8, dimension (:, :, :), allocatable :: faxc 
  real*8, dimension (:, :, :), allocatable :: faxc_ca 
  real*8, dimension (:, :, :), allocatable :: dxcdcc
  real*8, dimension (:, :), allocatable :: ft
  real*8, dimension (:, :), allocatable :: dusr
  real*8, dimension (:, :), allocatable :: special_k_orig
  real*8, dimension (:), allocatable :: weight_k_orig
  real*8, dimension (:, :), allocatable :: scale_k
  real*8, dimension (:, :, :), allocatable :: fotxc 
  real*8, dimension (:, :, :), allocatable :: fotxc_ca
  real*8, dimension (:, :, :), allocatable :: faca 
  real*8, dimension (:, :, :), allocatable :: fotca
  real*8, dimension (:, :), allocatable :: f3naa, f3nab, f3nac
  real*8, dimension (:, :), allocatable :: f3nla, f3nlb, f3nlc
  real*8, dimension (:, :), allocatable :: f3caa, f3cab, f3cac
  real*8, dimension (:, :), allocatable :: flrew
  real*8, dimension (:, :), allocatable :: f3xca_ca, f3xcb_ca, f3xcc_ca
  real*8, dimension (:, :), allocatable :: f3xca, f3xcb, f3xcc
  real*8, dimension (:, :), allocatable :: flrew_qmmm
  real*8, dimension (:, :), allocatable :: fro
  real*8, dimension (:, :), allocatable :: ftot
  real*8, dimension (:, :), allocatable :: dxcv

end module
