module M_system
  use M_constants, only: wp

  !========================
  integer :: iforce     = 1
  integer :: iqout      != 2 ! 1:Lowdin 2:Mulliken 3:NPA 4:M-dipole :7MD-pres..
  integer :: icluster   != 1 ! 
  integer :: idipole    != 1 
  integer :: igamma     != 1
  integer :: iqmmm      = 0
  integer :: ifixcharge = 0
  !========================
  !               idipole icluster igamma 
  !molecule       1       1        1        !R non periodic + only gamma kpts
  !molecule_test  0       1        1        !R solo para comparar, hasta que idipole este en sistemas periodicos
  !periodic       0       0        0        !I periodic + kpoints
  !periodic_gamma 0       0        1        !R periodic + only gamma kpts

  integer :: ialgmix  = 1 !1:anderson 2:johnson 3:custom
  real(wp), parameter :: xc_overtol = 5.0d-5
  real(wp), parameter :: smt_elect = 0.8d0 ! Ewald and electrostatic
  integer, parameter :: ithetamax = 5
  integer :: idmix = 6
  integer :: max_scf_iterations = 200
  real(wp) :: tempfe = 100.0d0
  real(wp) :: bmix = 0.1
  real(wp) :: sigma = 0.0d0
  real(wp) :: sigmaold = 0.0d0
  real(wp) :: sigmatol = 1.0E-8
  real(wp) :: sigmabest
  real(wp) :: w02 = 0.0
  logical :: scf_achieved = .false.
 
  integer :: natoms
  real(wp), dimension (:, :), allocatable :: ratom
  integer, dimension (:), allocatable :: degelec
  integer, dimension (:), allocatable :: imass
  character (len = 2), dimension (:), allocatable:: symbol
  integer :: mbeta_max  
  integer :: neigh_max
  integer :: ishiftO
  real(wp), dimension (3) :: shifter

  !--- mandar a M_fdata ? AQUI
  integer, dimension (:), allocatable :: getmssh
  integer, dimension (:), allocatable :: getlssh
  integer, dimension (:), allocatable :: getissh

  ! --- ETOT ----
  real(wp) :: etot
  real(wp) :: etotold, etotnew
  real(wp) :: etotper
  real(wp) :: atomic_energy, efermi
  real(wp) :: uiiuee
  real(wp) :: uxcdcc
  real(wp) :: uxcdcc_ols
  real(wp) :: etotxc_1c
  real(wp) :: ebs
  real(wp) :: eqmmm
  real(wp) :: dc_v_intra_dip_1c
  integer :: Kscf  
  real(wp) :: Uexc_1c
  real(wp) ::  Umuxc_1c
  real(wp), dimension (:, :, :, :), allocatable :: vxc_1c

  !Charges
  real(wp), dimension(:), allocatable  :: Q0_TOT
  integer, dimension (:), allocatable :: nelectron

  !--diag--
  real(wp), dimension (:, :, :), allocatable :: blowre
  real(wp), dimension (:, :, :), allocatable :: bbnkre
  real(wp), dimension (:, :, :), allocatable :: blowim
  real(wp), dimension (:, :, :), allocatable :: bbnkim
  real(wp), dimension (:, :), allocatable :: sm12_real
  complex(wp), dimension (:, :, :), allocatable :: sm12_complex
  

  real(wp), dimension (:, :), allocatable :: eigen_k
  real(wp), dimension (:, :), allocatable :: special_k
  integer :: norbitals
  integer :: norbitals_new
  integer :: nkpoints  
  integer, dimension (:), allocatable :: getiatom
  integer, dimension (:,:),allocatable :: ioccupy_k !AQUI pensar allo en denmat
  integer, dimension (:), allocatable :: ioccupy   !AQUI
  real(wp), dimension (:,:), allocatable :: foccupy !AQUI
  
  !scf
  real(wp), dimension (:, :, :, :), allocatable :: cape
  real(wp), dimension (:, :, :, :), allocatable :: rhoPP
  real(wp) :: ztot
  real(wp), dimension (:), allocatable :: weight_k
  integer ::  nssh_tot
  real(wp), dimension (:), allocatable :: mwe
  real(wp), dimension (:), allocatable :: drwe

  ! -- EWALD--
  real(wp), dimension (:, :), allocatable :: ewald
  real(wp), dimension (:, :, :), allocatable :: dewald
  real(wp), dimension (:, :), allocatable :: fewald
  real(wp), dimension (:, :, :, :), allocatable :: ewaldsr
  real(wp), dimension (:, :, :, :), allocatable :: dip
  real(wp), dimension (:, :, :, :, :), allocatable :: dipp
 
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
  real(wp), dimension (:, :, :, :), allocatable :: sVNL
  real(wp), dimension (:, :, :, :, :), allocatable :: spVNL
  real(wp), dimension (:, :, :, :, :), allocatable :: sp_mat
  real(wp), dimension (:, :, :, :, :), allocatable :: tp_mat
  real(wp), dimension (:, :, :), allocatable :: dipcm
  real(wp), dimension (:, :, :, :), allocatable :: dippcm
  real(wp), dimension (:, :, :, :, :, :), allocatable :: dippc
  real(wp), dimension(3,3) :: eps2  
  real(wp), dimension (:, :, :, :), allocatable :: vnl
  integer, dimension (:), allocatable :: neighPP_comn
  integer, dimension (:, :), allocatable :: neighPP_comm
  integer, dimension (:, :, :), allocatable :: neighPP_comj
  integer, dimension (:, :, :), allocatable :: neighPP_comb 
  integer,allocatable   :: neighj_aux(:,:)
  integer :: neighPP_max 
  integer :: num_neig_maxtot
 
  !CHARGES
  real(wp), dimension (:, :), allocatable :: Qin
  real(wp), dimension (:), allocatable :: Qinmixer
  real(wp), dimension (:, :), allocatable :: Qout
  real(wp), dimension (:), allocatable :: Qoutmixer
  real(wp), dimension (:), allocatable :: dq
  real(wp), dimension(:), allocatable  :: Q_partial
  real(wp), dimension (:, :), allocatable :: Qin_es
  real(wp), dimension (:), allocatable :: QLowdin_TOT
  real(wp), dimension (:), allocatable :: QMulliken_TOT
  real(wp), allocatable, dimension(:,:) :: qaux
  real(wp)    dip_x, dipQout_x, dipTot_x, dipProy_x, dipIntra_x, dip_res_x, dipQin_x, dipRes_x
  real(wp)    dip_y, dipQout_y, dipTot_y, dipProy_y, dipIntra_y, dip_res_y, dipQin_y, dipRes_y
  real(wp)    dip_z, dipQout_z, dipTot_z, dipProy_z, dipIntra_z, dip_res_z, dipQin_z, dipRes_z
  real(wp)    dip_tot, dip_proy, dipQin_tot, dipTot_tot, dipIntra_tot, dipQout_tot, dip_res_tot, dipRes_tot 
  real(wp), dimension (:), allocatable :: dq_DP

  !interaccions
  real(wp), dimension (:, :, :, :), allocatable :: vxc
  real(wp), dimension (:, :, :, :), allocatable :: vxc_ca
  real(wp), dimension (:, :, :, :), allocatable :: rho
  real(wp), dimension (:, :, :, :), allocatable :: rho_off
  real(wp), dimension (:, :, :, :), allocatable :: rhoij_off
  real(wp), dimension (:, :, :, :), allocatable :: s_mat 
  real(wp), dimension (:, :, :, :), allocatable :: sm_mat
  real(wp), dimension (:, :, :, :, :), allocatable :: spm_mat
  real(wp), dimension (:, :, :), allocatable :: rho_on
  real(wp), dimension (:, :, :), allocatable :: arho_on
  real(wp), dimension (:, :, :), allocatable :: rhoi_on
  real(wp), dimension (:, :, :), allocatable :: arhoi_on
  real(wp), dimension (:, :, :, :, :), allocatable :: arhop_on
  real(wp), dimension (:, :, :, :, :), allocatable :: rhop_on
  real(wp), dimension (:, :, :, :), allocatable :: arhoij_off
  real(wp), dimension (:, :, :, :), allocatable :: arho_off
  real(wp), dimension (:, :, :, :, :), allocatable :: arhopij_off 
  real(wp), dimension (:, :, :, :, :), allocatable :: arhop_off
  real(wp), dimension (:, :, :, :, :), allocatable :: rhop_off
  real(wp), dimension (:, :, :, :, :), allocatable :: rhopij_off
  real(wp), dimension (:, :, :, :), allocatable :: vca
  real(wp), dimension (:, :, :, :), allocatable :: ewaldlr
  real(wp), dimension (:, :, :, :), allocatable :: h_mat
  real(wp), dimension (:, :, :, :), allocatable :: t_mat
  real(wp), dimension (:, :, :, :), allocatable :: vna
  real(wp), dimension (:, :, :, :), allocatable :: ewaldqmmm
  real(wp), dimension (:, :, :, :, :), allocatable :: dipc
  integer, dimension(:,:), allocatable :: muR
  integer, dimension(:,:), allocatable :: nuR
  integer, dimension(:,:), allocatable :: alphaR
  integer, dimension(:,:), allocatable :: betaR
  real(wp), dimension(:,:,:,:), allocatable   :: hr_box

  
  real(wp), allocatable, dimension(:,:) :: Fv   
  real(wp), allocatable, dimension(:,:) :: Xv    
  real(wp), allocatable, dimension(:,:) :: delF 
  real(wp), allocatable, dimension(:,:) :: delX 
  real(wp), allocatable, dimension(:)   :: r2_sav
  real(wp), allocatable, dimension(:)   :: wi
  real(wp), allocatable, dimension(:)   :: x_best
 
  real(wp), allocatable, dimension(:,:) :: RJac
  real(wp), allocatable, dimension(:,:) :: betaInvH
  real(wp), allocatable, dimension(:,:) :: gamaH
  real(wp), dimension (:,:), allocatable :: xl

  ! Lattice vectors
  real(wp), dimension (3) :: a1vec
  real(wp), dimension (3) :: a2vec
  real(wp), dimension (3) :: a3vec

  ! allocate_f
  real(wp), dimension (:, :, :), allocatable :: fotnl
  real(wp), dimension (:, :, :), allocatable :: fanl
  real(wp), dimension (:, :, :), allocatable :: fotna
  real(wp), dimension (:, :, :), allocatable :: fana
  real(wp), dimension (:, :, :), allocatable :: faxc 
  real(wp), dimension (:, :, :), allocatable :: faxc_ca 
  real(wp), dimension (:, :, :), allocatable :: dxcdcc
  real(wp), dimension (:, :), allocatable :: ft
  real(wp), dimension (:, :), allocatable :: dusr
  real(wp), dimension (:, :), allocatable :: special_k_orig
  real(wp), dimension (:), allocatable :: weight_k_orig
  real(wp), dimension (:, :), allocatable :: scale_k
  real(wp), dimension (:, :, :), allocatable :: fotxc 
  real(wp), dimension (:, :, :), allocatable :: fotxc_ca
  real(wp), dimension (:, :, :), allocatable :: faca 
  real(wp), dimension (:, :, :), allocatable :: fotca
  real(wp), dimension (:, :), allocatable :: f3naa, f3nab, f3nac
  real(wp), dimension (:, :), allocatable :: f3nla, f3nlb, f3nlc
  real(wp), dimension (:, :), allocatable :: f3caa, f3cab, f3cac
  real(wp), dimension (:, :), allocatable :: flrew
  real(wp), dimension (:, :), allocatable :: f3xca_ca, f3xcb_ca, f3xcc_ca
  real(wp), dimension (:, :), allocatable :: f3xca, f3xcb, f3xcc
  real(wp), dimension (:, :), allocatable :: flrew_qmmm
  real(wp), dimension (:, :), allocatable :: fro
  real(wp), dimension (:, :), allocatable :: ftot
  real(wp), dimension (:, :), allocatable :: dxcv


  integer                        :: qmmm_qm_mm_pairs     !Number of pairs per QM atom. - length of pair_list. 
  real(wp), dimension(:,:), allocatable :: qmmm_dxyzcl  !Used to store the forces generated by qm_mm before adding them to the main f array.
  real(wp), dimension(:,:),  allocatable ::  qmmm_qm_xcrd         !Contains imaged mm coordinates and scaled mm charges.
  real(wp), dimension(:), allocatable :: qmmm_scf_mchg
  real(wp), dimension(:), allocatable :: qmmm_Qneutral_TOT


  real(wp) amat (3, 3, 5) 

end module
