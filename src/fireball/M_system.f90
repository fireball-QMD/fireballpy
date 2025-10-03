module M_system
  use, intrinsic :: iso_fortran_env, only: double => real64

  !========================
  integer :: iforce     = 1
  integer :: iqout      != 2 ! 1:Lowdin 2:Mulliken 3:NPA 4:M-dipole :7MD-pres..
  integer :: icluster   != 1 ! 
  integer :: idipole    != 1 
  integer :: igamma     != 1
  integer :: iqmmm      = 0
  integer :: ifixcharge
  integer :: qstate = 0
  !========================
  !               idipole icluster igamma 
  !molecule       1       1        1        !R non periodic + only gamma kpts
  !molecule_test  0       1        1        !R solo para comparar, hasta que idipole este en sistemas periodicos
  !periodic       0       0        0        !I periodic + kpoints
  !periodic_gamma 0       0        1        !R periodic + only gamma kpts
  !
  logical :: isgeneig

  integer :: errno = 0

  integer :: ialgmix  = 2 !1:anderson 2:johnson 3:custom
  real(double), parameter :: xc_overtol = 5.0d-5
  real(double), parameter :: smt_elect = 0.8d0 ! Ewald and electrostatic
  integer, parameter :: ithetamax = 5
  integer :: idmix = 6
  integer :: max_scf_iterations = 200
  real(double) :: tempfe = 100.0d0
  real(double) :: bmix = 0.1d0
  integer :: Kbest
  real(double) :: sigma = 0.0d0
  real(double) :: sigmaold = 0.0d0
  real(double) :: sigmatol = 1.0d-8
  real(double) :: sigmabest
  real(double) :: w02 = 0.0d0
  logical :: scf_achieved = .false.
 
  integer :: natoms
  real(double), dimension (:, :), allocatable :: ratom
  integer, dimension (:), allocatable :: degelec
  integer, dimension (:), allocatable :: imass
  character (len = 2), dimension (:), allocatable:: symbol
  integer :: mbeta_max  
  integer :: neigh_max
  !integer :: ishiftO
  !real(double), dimension (3) :: shifter

  !--- mandar a M_fdata ? AQUI
  integer, dimension (:), allocatable :: getmssh
  integer, dimension (:), allocatable :: getlssh
  integer, dimension (:), allocatable :: getissh

  ! --- ETOT ----
  real(double) :: etot
  real(double) :: etotold, etotnew
  real(double) :: etotper
  real(double) :: atomic_energy, efermi
  real(double) :: uiiuee
  real(double) :: uxcdcc
  real(double) :: uxcdcc_ols
  real(double) :: etotxc_1c
  real(double) :: ebs
  real(double) :: eqmmm
  real(double) :: dc_v_intra_dip_1c
  integer :: Kscf  
  real(double) :: Uexc_1c
  real(double) ::  Umuxc_1c
  real(double), dimension (:, :, :, :), allocatable :: vxc_1c

  !Charges
  real(double), dimension(:), allocatable  :: Q0_TOT
  integer, dimension (:), allocatable :: nelectron

  !--diag--
  real(double), dimension (:, :, :), allocatable :: blowre
  real(double), dimension (:, :, :), allocatable :: bbnkre
  real(double), dimension (:, :, :), allocatable :: blowim
  real(double), dimension (:, :, :), allocatable :: bbnkim
  real(double), dimension (:, :), allocatable :: sm12_real
  complex(double), dimension (:, :, :), allocatable :: sm12_complex
  

  real(double), dimension (:, :), allocatable :: eigen_k
  real(double), dimension (:, :), allocatable :: special_k
  integer :: norbitals
  integer :: norbitals_new
  integer :: nkpoints  
  integer, dimension (:), allocatable :: getiatom
  integer, dimension (:,:),allocatable :: ioccupy_k !AQUI pensar allo en denmat
  integer, dimension (:), allocatable :: ioccupy   !AQUI
  real(double), dimension (:,:), allocatable :: foccupy !AQUI
  
  !scf
  real(double), dimension (:, :, :, :), allocatable :: cape
  real(double), dimension (:, :, :, :), allocatable :: rhoPP
  real(double) :: ztot
  real(double), dimension (:), allocatable :: weight_k
  integer ::  nssh_tot
  real(double), dimension (:), allocatable :: mwe
  real(double), dimension (:), allocatable :: drwe

  ! -- EWALD--
  real(double), dimension (:, :), allocatable :: ewald
  real(double), dimension (:, :, :), allocatable :: dewald
  real(double), dimension (:, :), allocatable :: fewald
  real(double), dimension (:, :, :, :), allocatable :: ewaldsr
  real(double), dimension (:, :, :, :), allocatable :: dip
  real(double), dimension (:, :, :, :, :), allocatable :: dipp
 
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
  real(double), dimension (:, :, :, :), allocatable :: sVNL
  real(double), dimension (:, :, :, :, :), allocatable :: spVNL
  real(double), dimension (:, :, :, :, :), allocatable :: sp_mat
  real(double), dimension (:, :, :, :, :), allocatable :: tp_mat
  real(double), dimension (:, :, :), allocatable :: dipcm
  real(double), dimension (:, :, :, :), allocatable :: dippcm
  real(double), dimension (:, :, :, :, :, :), allocatable :: dippc
  real(double), dimension(3,3) :: eps2  
  real(double), dimension (:, :, :, :), allocatable :: vnl
  integer, dimension (:), allocatable :: neighPP_comn
  integer, dimension (:, :), allocatable :: neighPP_comm
  integer, dimension (:, :, :), allocatable :: neighPP_comj
  integer, dimension (:, :, :), allocatable :: neighPP_comb 
  integer, dimension(:, :), allocatable   :: neighj_aux
  integer :: neighPP_max 
  integer :: num_neig_maxtot
 
  !CHARGES
  real(double), dimension (:, :), allocatable :: Qin
  real(double), dimension (:), allocatable :: Qinmixer
  real(double), dimension (:, :), allocatable :: Qout
  real(double), dimension (:), allocatable :: Qoutmixer
  real(double), dimension (:), allocatable :: dq
  real(double), dimension(:), allocatable  :: Q_partial
  real(double), dimension (:, :), allocatable :: Qin_es
  real(double), dimension (:), allocatable :: QLowdin_TOT
  real(double), dimension (:), allocatable :: QMulliken_TOT
  real(double), allocatable, dimension(:,:) :: qaux
  real(double) :: dip_x, dipQout_x, dipTot_x, dipProy_x, dipIntra_x, dip_res_x, dipQin_x, dipRes_x
  real(double) :: dip_y, dipQout_y, dipTot_y, dipProy_y, dipIntra_y, dip_res_y, dipQin_y, dipRes_y
  real(double) :: dip_z, dipQout_z, dipTot_z, dipProy_z, dipIntra_z, dip_res_z, dipQin_z, dipRes_z
  real(double) :: dip_tot, dip_proy, dipQin_tot, dipTot_tot, dipIntra_tot, dipQout_tot, dip_res_tot, dipRes_tot 
  real(double), dimension (:), allocatable :: dq_DP

  !interaccions
  real(double), dimension (:, :, :, :), allocatable :: vxc
  real(double), dimension (:, :, :, :), allocatable :: vxc_ca
  real(double), dimension (:, :, :, :), allocatable :: rho
  real(double), dimension (:, :, :, :), allocatable :: rho_off
  real(double), dimension (:, :, :, :), allocatable :: rhoij_off
  real(double), dimension (:, :, :, :), allocatable :: s_mat 
  real(double), dimension (:, :, :, :), allocatable :: sm_mat
  real(double), dimension (:, :, :, :, :), allocatable :: spm_mat
  real(double), dimension (:, :, :), allocatable :: rho_on
  real(double), dimension (:, :, :), allocatable :: arho_on
  real(double), dimension (:, :, :), allocatable :: rhoi_on
  real(double), dimension (:, :, :), allocatable :: arhoi_on
  real(double), dimension (:, :, :, :, :), allocatable :: arhop_on
  real(double), dimension (:, :, :, :, :), allocatable :: rhop_on
  real(double), dimension (:, :, :, :), allocatable :: arhoij_off
  real(double), dimension (:, :, :, :), allocatable :: arho_off
  real(double), dimension (:, :, :, :, :), allocatable :: arhopij_off 
  real(double), dimension (:, :, :, :, :), allocatable :: arhop_off
  real(double), dimension (:, :, :, :, :), allocatable :: rhop_off
  real(double), dimension (:, :, :, :, :), allocatable :: rhopij_off
  real(double), dimension (:, :, :, :), allocatable :: vca
  real(double), dimension (:, :, :, :), allocatable :: ewaldlr
  real(double), dimension (:, :, :, :), allocatable :: h_mat
  real(double), dimension (:, :, :, :), allocatable :: t_mat
  real(double), dimension (:, :, :, :), allocatable :: vna
  real(double), dimension (:, :, :, :), allocatable :: ewaldqmmm
  real(double), dimension (:, :, :, :, :), allocatable :: dipc
  integer, dimension(:,:), allocatable :: muR
  integer, dimension(:,:), allocatable :: nuR
  integer, dimension(:,:), allocatable :: alphaR
  integer, dimension(:,:), allocatable :: betaR
  real(double), dimension(:,:,:,:), allocatable   :: hr_box

  
  real(double), allocatable, dimension(:,:) :: Fv   
  real(double), allocatable, dimension(:,:) :: Xv    
  real(double), allocatable, dimension(:,:) :: delF 
  real(double), allocatable, dimension(:,:) :: delX 
  real(double), allocatable, dimension(:)   :: r2_sav
  real(double), allocatable, dimension(:)   :: wi
  real(double), allocatable, dimension(:)   :: x_best
 
  real(double), allocatable, dimension(:,:) :: RJac
  real(double), allocatable, dimension(:,:) :: betaInvH
  real(double), allocatable, dimension(:,:) :: gamaH
  real(double), dimension (:,:), allocatable :: xl

  ! Lattice vectors
  real(double), dimension (3) :: a1vec
  real(double), dimension (3) :: a2vec
  real(double), dimension (3) :: a3vec

  ! allocate_f
  real(double), dimension (:, :, :), allocatable :: fotnl
  real(double), dimension (:, :, :), allocatable :: fanl
  real(double), dimension (:, :, :), allocatable :: fotna
  real(double), dimension (:, :, :), allocatable :: fana
  real(double), dimension (:, :, :), allocatable :: faxc 
  real(double), dimension (:, :, :), allocatable :: faxc_ca 
  real(double), dimension (:, :, :), allocatable :: dxcdcc
  real(double), dimension (:, :), allocatable :: ft
  real(double), dimension (:, :), allocatable :: dusr
  real(double), dimension (:, :, :), allocatable :: fotxc 
  real(double), dimension (:, :, :), allocatable :: fotxc_ca
  real(double), dimension (:, :, :), allocatable :: faca 
  real(double), dimension (:, :, :), allocatable :: fotca
  real(double), dimension (:, :), allocatable :: f3naa, f3nab, f3nac
  real(double), dimension (:, :), allocatable :: f3nla, f3nlb, f3nlc
  real(double), dimension (:, :), allocatable :: f3caa, f3cab, f3cac
  real(double), dimension (:, :), allocatable :: flrew
  real(double), dimension (:, :), allocatable :: f3xca_ca, f3xcb_ca, f3xcc_ca
  real(double), dimension (:, :), allocatable :: f3xca, f3xcb, f3xcc
  real(double), dimension (:, :), allocatable :: flrew_qmmm
  real(double), dimension (:, :), allocatable :: fro
  real(double), dimension (:, :), allocatable :: ftot
  real(double), dimension (:, :), allocatable :: dxcv


  !integer                        :: qmmm_qm_mm_pairs     !Number of pairs per QM atom. - length of pair_list. 
  integer                        :: qmmm_qm_natoms
  real(double), dimension(:,:), allocatable :: qmmm_dxyzcl  !Used to store the forces generated by qm_mm before adding them to the main f array.
  real(double), dimension(:,:),  allocatable ::  qmmm_qm_xcrd         !Contains imaged mm coordinates and scaled mm charges.
  real(double), dimension(:), allocatable :: qmmm_scf_mchg
  real(double), dimension(:), allocatable :: qmmm_Qneutral_TOT
  real(double) :: qmmm_rc1, qmmm_rc2, qmmm_width


  real(double), dimension(3, 3, 5) :: amat
  
  real(double), dimension (:,:,:,:,:,:), allocatable :: g_h
  real(double), dimension (:,:,:,:,:,:), allocatable :: g_xc  
  real(double), dimension (:,:,:,:,:), allocatable :: f_xc
  real(double), dimension (:,:,:,:), allocatable :: exc_aa
  real(double), dimension (:,:,:,:), allocatable :: vxc_aa
  integer, dimension (:), allocatable :: fix_shell_charge 
end module M_system
