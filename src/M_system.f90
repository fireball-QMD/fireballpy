module M_system
  use iso_c_binding

  !========================
  integer(c_long) :: iforce     = 1
  integer(c_long) :: iqout      != 2 ! 1:Lowdin 2:Mulliken 3:NPA 4:M-dipole :7MD-pres..
  integer(c_long) :: icluster   != 1 ! 
  integer(c_long) :: idipole    != 1 
  integer(c_long) :: igamma     != 1
  integer(c_long) :: iqmmm      = 0
  integer(c_long) :: ifixcharge = 0
  !========================
  !               idipole icluster igamma 
  !molecule       1       1        1        !R non periodic + only gamma kpts
  !molecule_test  0       1        1        !R solo para comparar, hasta que idipole este en sistemas periodicos
  !periodic       0       0        0        !I periodic + kpoints
  !periodic_gamma 0       0        1        !R periodic + only gamma kpts

  integer(c_long) :: ialgmix  = 1 !1:anderson 2:johnson 3:custom
  real(c_double), parameter :: xc_overtol = 5.0d-5
  real(c_double), parameter :: smt_elect = 0.8d0 ! Ewald and electrostatic
  integer(c_long), parameter :: ithetamax = 5
  integer(c_long) :: idmix = 6
  integer(c_long) :: max_scf_iterations = 200
  real(c_double) :: tempfe = 100.0d0
  real(c_double) :: bmix = 0.1d0
  real(c_double) :: sigma = 0.0d0
  real(c_double) :: sigmaold = 0.0d0
  real(c_double) :: sigmatol = 1.0E-8
  real(c_double) :: sigmabest
  real(c_double) :: w02 = 0.0d0
  logical :: scf_achieved = .false.
 
  integer(c_long) :: natoms
  real(c_double), dimension (:, :), allocatable :: ratom
  integer(c_long), dimension (:), allocatable :: degelec
  integer(c_long), dimension (:), allocatable :: imass
  character (len = 2), dimension (:), allocatable:: symbol
  integer(c_long) :: mbeta_max  
  integer(c_long) :: neigh_max
  integer(c_long) :: ishiftO
  real(c_double), dimension (3) :: shifter

  !--- mandar a M_fdata ? AQUI
  integer(c_long), dimension (:), allocatable :: getmssh
  integer(c_long), dimension (:), allocatable :: getlssh
  integer(c_long), dimension (:), allocatable :: getissh

  ! --- ETOT ----
  real(c_double) :: etot
  real(c_double) :: etotold, etotnew
  real(c_double) :: etotper
  real(c_double) :: atomic_energy, efermi
  real(c_double) :: uiiuee
  real(c_double) :: uxcdcc
  real(c_double) :: uxcdcc_ols
  real(c_double) :: etotxc_1c
  real(c_double) :: ebs
  real(c_double) :: eqmmm
  real(c_double) :: dc_v_intra_dip_1c
  integer(c_long) :: Kscf  
  real(c_double) :: Uexc_1c
  real(c_double) ::  Umuxc_1c
  real(c_double), dimension (:, :, :, :), allocatable :: vxc_1c

  !Charges
  real(c_double), dimension(:), allocatable  :: Q0_TOT
  integer(c_long), dimension (:), allocatable :: nelectron

  !--diag--
  real(c_double), dimension (:, :, :), allocatable :: blowre
  real(c_double), dimension (:, :, :), allocatable :: bbnkre
  real(c_double), dimension (:, :, :), allocatable :: blowim
  real(c_double), dimension (:, :, :), allocatable :: bbnkim
  real(c_double), dimension (:, :), allocatable :: sm12_real
  complex(c_double_complex), dimension (:, :, :), allocatable :: sm12_complex
  

  real(c_double), dimension (:, :), allocatable :: eigen_k
  real(c_double), dimension (:, :), allocatable :: special_k
  integer(c_long) :: norbitals
  integer(c_long) :: norbitals_new
  integer(c_long) :: nkpoints  
  integer(c_long), dimension (:), allocatable :: getiatom
  integer(c_long), dimension (:,:),allocatable :: ioccupy_k !AQUI pensar allo en denmat
  integer(c_long), dimension (:), allocatable :: ioccupy   !AQUI
  real(c_double), dimension (:,:), allocatable :: foccupy !AQUI
  
  !scf
  real(c_double), dimension (:, :, :, :), allocatable :: cape
  real(c_double), dimension (:, :, :, :), allocatable :: rhoPP
  real(c_double) :: ztot
  real(c_double), dimension (:), allocatable :: weight_k
  integer(c_long) ::  nssh_tot
  real(c_double), dimension (:), allocatable :: mwe
  real(c_double), dimension (:), allocatable :: drwe

  ! -- EWALD--
  real(c_double), dimension (:, :), allocatable :: ewald
  real(c_double), dimension (:, :, :), allocatable :: dewald
  real(c_double), dimension (:, :), allocatable :: fewald
  real(c_double), dimension (:, :, :, :), allocatable :: ewaldsr
  real(c_double), dimension (:, :, :, :), allocatable :: dip
  real(c_double), dimension (:, :, :, :, :), allocatable :: dipp
 
  ! --- NEIGH ----
  integer(c_long), dimension (:, :), allocatable :: neigh_b  
  integer(c_long), dimension (:, :), allocatable :: neigh_j
  integer(c_long), dimension (:), allocatable :: neighn
  integer(c_long), dimension (:, :, :), allocatable :: neigh_comb
  integer(c_long), dimension (:, :, :), allocatable :: neigh_comj
  integer(c_long), dimension (:, :, :), allocatable :: neigh_com_ng
  integer(c_long), dimension (:, :), allocatable :: neigh_comm
  integer(c_long), dimension (:), allocatable :: neigh_comn 
  integer(c_long), dimension (:,:), allocatable :: neigh_back
  integer(c_long), dimension (:), allocatable :: neigh_self
  integer(c_long), dimension (:, :), allocatable :: nPP_b
  integer(c_long), dimension (:, :), allocatable :: nPP_j
  integer(c_long), dimension (:, :), allocatable :: nPP_map
  integer(c_long), dimension (:), allocatable :: nPPn
  integer(c_long), dimension (:), allocatable :: nPP_self
  integer(c_long), dimension (:, :), allocatable :: nPPx_b
  integer(c_long), dimension (:, :), allocatable :: nPPx_j
  integer(c_long), dimension (:, :), allocatable :: nPPx_map
  integer(c_long), dimension (:, :), allocatable :: nPPx_point
  integer(c_long), dimension (:), allocatable :: nPPxn
  integer(c_long), dimension (:), allocatable :: nPPx_self
  integer(c_long), dimension(:), allocatable :: neigh_pair_a1
  integer(c_long), dimension(:), allocatable :: neigh_pair_a2
  integer(c_long), dimension(:), allocatable :: neigh_pair_n1
  integer(c_long), dimension(:), allocatable :: neigh_pair_n2
  integer(c_long), dimension(:,:), allocatable    :: neighj_tot
  integer(c_long), dimension(:,:), allocatable    :: neighb_tot
  integer(c_long), dimension(:), allocatable      :: neighn_tot
  integer(c_long) :: numorb_max
  integer(c_long), dimension (:), allocatable :: neighPP_self
  integer(c_long), dimension (:), allocatable :: neighPPn
  integer(c_long), dimension (:, :), allocatable :: neighPP_b
  integer(c_long), dimension (:, :), allocatable :: neighPP_j
  integer(c_long) :: tot_pairs
  real(c_double), dimension (:, :, :, :), allocatable :: sVNL
  real(c_double), dimension (:, :, :, :, :), allocatable :: spVNL
  real(c_double), dimension (:, :, :, :, :), allocatable :: sp_mat
  real(c_double), dimension (:, :, :, :, :), allocatable :: tp_mat
  real(c_double), dimension (:, :, :), allocatable :: dipcm
  real(c_double), dimension (:, :, :, :), allocatable :: dippcm
  real(c_double), dimension (:, :, :, :, :, :), allocatable :: dippc
  real(c_double), dimension(3,3) :: eps2  
  real(c_double), dimension (:, :, :, :), allocatable :: vnl
  integer(c_long), dimension (:), allocatable :: neighPP_comn
  integer(c_long), dimension (:, :), allocatable :: neighPP_comm
  integer(c_long), dimension (:, :, :), allocatable :: neighPP_comj
  integer(c_long), dimension (:, :, :), allocatable :: neighPP_comb 
  integer(c_long), dimension(:, :), allocatable   :: neighj_aux
  integer(c_long) :: neighPP_max 
  integer(c_long) :: num_neig_maxtot
 
  !CHARGES
  real(c_double), dimension (:, :), allocatable :: Qin
  real(c_double), dimension (:), allocatable :: Qinmixer
  real(c_double), dimension (:, :), allocatable :: Qout
  real(c_double), dimension (:), allocatable :: Qoutmixer
  real(c_double), dimension (:), allocatable :: dq
  real(c_double), dimension(:), allocatable  :: Q_partial
  real(c_double), dimension (:, :), allocatable :: Qin_es
  real(c_double), dimension (:), allocatable :: QLowdin_TOT
  real(c_double), dimension (:), allocatable :: QMulliken_TOT
  real(c_double), allocatable, dimension(:,:) :: qaux
  real(c_double) :: dip_x, dipQout_x, dipTot_x, dipProy_x, dipIntra_x, dip_res_x, dipQin_x, dipRes_x
  real(c_double) :: dip_y, dipQout_y, dipTot_y, dipProy_y, dipIntra_y, dip_res_y, dipQin_y, dipRes_y
  real(c_double) :: dip_z, dipQout_z, dipTot_z, dipProy_z, dipIntra_z, dip_res_z, dipQin_z, dipRes_z
  real(c_double) :: dip_tot, dip_proy, dipQin_tot, dipTot_tot, dipIntra_tot, dipQout_tot, dip_res_tot, dipRes_tot 
  real(c_double), dimension (:), allocatable :: dq_DP

  !interaccions
  real(c_double), dimension (:, :, :, :), allocatable :: vxc
  real(c_double), dimension (:, :, :, :), allocatable :: vxc_ca
  real(c_double), dimension (:, :, :, :), allocatable :: rho
  real(c_double), dimension (:, :, :, :), allocatable :: rho_off
  real(c_double), dimension (:, :, :, :), allocatable :: rhoij_off
  real(c_double), dimension (:, :, :, :), allocatable :: s_mat 
  real(c_double), dimension (:, :, :, :), allocatable :: sm_mat
  real(c_double), dimension (:, :, :, :, :), allocatable :: spm_mat
  real(c_double), dimension (:, :, :), allocatable :: rho_on
  real(c_double), dimension (:, :, :), allocatable :: arho_on
  real(c_double), dimension (:, :, :), allocatable :: rhoi_on
  real(c_double), dimension (:, :, :), allocatable :: arhoi_on
  real(c_double), dimension (:, :, :, :, :), allocatable :: arhop_on
  real(c_double), dimension (:, :, :, :, :), allocatable :: rhop_on
  real(c_double), dimension (:, :, :, :), allocatable :: arhoij_off
  real(c_double), dimension (:, :, :, :), allocatable :: arho_off
  real(c_double), dimension (:, :, :, :, :), allocatable :: arhopij_off 
  real(c_double), dimension (:, :, :, :, :), allocatable :: arhop_off
  real(c_double), dimension (:, :, :, :, :), allocatable :: rhop_off
  real(c_double), dimension (:, :, :, :, :), allocatable :: rhopij_off
  real(c_double), dimension (:, :, :, :), allocatable :: vca
  real(c_double), dimension (:, :, :, :), allocatable :: ewaldlr
  real(c_double), dimension (:, :, :, :), allocatable :: h_mat
  real(c_double), dimension (:, :, :, :), allocatable :: t_mat
  real(c_double), dimension (:, :, :, :), allocatable :: vna
  real(c_double), dimension (:, :, :, :), allocatable :: ewaldqmmm
  real(c_double), dimension (:, :, :, :, :), allocatable :: dipc
  integer(c_long), dimension(:,:), allocatable :: muR
  integer(c_long), dimension(:,:), allocatable :: nuR
  integer(c_long), dimension(:,:), allocatable :: alphaR
  integer(c_long), dimension(:,:), allocatable :: betaR
  real(c_double), dimension(:,:,:,:), allocatable   :: hr_box

  
  real(c_double), allocatable, dimension(:,:) :: Fv   
  real(c_double), allocatable, dimension(:,:) :: Xv    
  real(c_double), allocatable, dimension(:,:) :: delF 
  real(c_double), allocatable, dimension(:,:) :: delX 
  real(c_double), allocatable, dimension(:)   :: r2_sav
  real(c_double), allocatable, dimension(:)   :: wi
  real(c_double), allocatable, dimension(:)   :: x_best
 
  real(c_double), allocatable, dimension(:,:) :: RJac
  real(c_double), allocatable, dimension(:,:) :: betaInvH
  real(c_double), allocatable, dimension(:,:) :: gamaH
  real(c_double), dimension (:,:), allocatable :: xl

  ! Lattice vectors
  real(c_double), dimension (3) :: a1vec
  real(c_double), dimension (3) :: a2vec
  real(c_double), dimension (3) :: a3vec

  ! allocate_f
  real(c_double), dimension (:, :, :), allocatable :: fotnl
  real(c_double), dimension (:, :, :), allocatable :: fanl
  real(c_double), dimension (:, :, :), allocatable :: fotna
  real(c_double), dimension (:, :, :), allocatable :: fana
  real(c_double), dimension (:, :, :), allocatable :: faxc 
  real(c_double), dimension (:, :, :), allocatable :: faxc_ca 
  real(c_double), dimension (:, :, :), allocatable :: dxcdcc
  real(c_double), dimension (:, :), allocatable :: ft
  real(c_double), dimension (:, :), allocatable :: dusr
  real(c_double), dimension (:, :), allocatable :: special_k_orig
  real(c_double), dimension (:), allocatable :: weight_k_orig
  real(c_double), dimension (:, :), allocatable :: scale_k
  real(c_double), dimension (:, :, :), allocatable :: fotxc 
  real(c_double), dimension (:, :, :), allocatable :: fotxc_ca
  real(c_double), dimension (:, :, :), allocatable :: faca 
  real(c_double), dimension (:, :, :), allocatable :: fotca
  real(c_double), dimension (:, :), allocatable :: f3naa, f3nab, f3nac
  real(c_double), dimension (:, :), allocatable :: f3nla, f3nlb, f3nlc
  real(c_double), dimension (:, :), allocatable :: f3caa, f3cab, f3cac
  real(c_double), dimension (:, :), allocatable :: flrew
  real(c_double), dimension (:, :), allocatable :: f3xca_ca, f3xcb_ca, f3xcc_ca
  real(c_double), dimension (:, :), allocatable :: f3xca, f3xcb, f3xcc
  real(c_double), dimension (:, :), allocatable :: flrew_qmmm
  real(c_double), dimension (:, :), allocatable :: fro
  real(c_double), dimension (:, :), allocatable :: ftot
  real(c_double), dimension (:, :), allocatable :: dxcv


  integer(c_long)                        :: qmmm_qm_mm_pairs     !Number of pairs per QM atom. - length of pair_list. 
  real(c_double), dimension(:,:), allocatable :: qmmm_dxyzcl  !Used to store the forces generated by qm_mm before adding them to the main f array.
  real(c_double), dimension(:,:),  allocatable ::  qmmm_qm_xcrd         !Contains imaged mm coordinates and scaled mm charges.
  real(c_double), dimension(:), allocatable :: qmmm_scf_mchg
  real(c_double), dimension(:), allocatable :: qmmm_Qneutral_TOT


  real(c_double), dimension(3, 3, 5) :: amat

end module
