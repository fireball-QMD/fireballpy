module M_system
  use, intrinsic :: iso_fortran_env, only: dp => real64

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
  real(dp), parameter :: xc_overtol = 5.0e-5_dp
  real(dp), parameter :: xc_overtolG = 5.0e-5_dp
  real(dp), parameter :: smt_elect = 0.8_dp ! Ewald and electrostatic
  integer, parameter :: ithetamax = 5
  integer :: idmix = 6
  integer :: max_scf_iterations = 200
  real(dp) :: tempfe = 100.0_dp
  real(dp) :: bmix = 0.1_dp
  integer :: Kbest
  real(dp) :: sigma = 0.0_dp
  real(dp) :: sigmaold = 0.0_dp
  real(dp) :: sigmatol = 1.0e-8_dp
  real(dp) :: sigmabest
  real(dp) :: w02 = 0.0_dp
  logical :: scf_achieved = .false.
 
  integer :: natoms
  real(dp), dimension (:, :), allocatable :: ratom
  integer, dimension (:), allocatable :: degelec
  integer, dimension (:), allocatable :: imass
  character (len = 2), dimension (:), allocatable:: symbol
  integer :: mbeta_max  
  integer :: neigh_max
  !integer :: ishiftO
  !real(dp), dimension (3) :: shifter

  !--- mandar a M_fdata ? AQUI
  integer, dimension (:), allocatable :: getmssh
  integer, dimension (:), allocatable :: getlssh
  integer, dimension (:), allocatable :: getissh
  integer, dimension (:), allocatable :: get_orb_ofshell
  integer, dimension (:), allocatable :: get_iatom_ofshell
  integer, dimension (:), allocatable :: get_issh_ofshell
  integer, dimension (:), allocatable :: get_l_ofshell 
  integer, dimension (:), allocatable :: get_shell_oforb
  integer, dimension (:,:), allocatable :: get_shell_ofatom_imu
  integer, dimension (:,:), allocatable :: get_shell_ofatom_issh

  ! --- ETOT ----
  real(dp) :: etot
  real(dp) :: etotold, etotnew
  real(dp) :: etotper
  real(dp) :: atomic_energy, efermi
  real(dp) :: uiiuee
  real(dp) :: uxcdcc
  real(dp) :: uxcdcc_zw
  real(dp) :: uxcdcc_ols
  real(dp) :: etotxc_1c
  real(dp) :: ebs
  real(dp) :: eqmmm
  real(dp) :: dc_v_intra_dip_1c
  integer :: Kscf  
  real(dp) :: Uexc_1c
  real(dp) ::  Umuxc_1c
  real(dp), dimension (:, :, :, :), allocatable :: vxc_1c

  !Charges
  real(dp), dimension(:), allocatable  :: Q0_TOT
  integer, dimension (:), allocatable :: nelectron

  !--diag--
  real(dp), dimension (:, :, :), allocatable :: blowre
  real(dp), dimension (:, :, :), allocatable :: bbnkre
  real(dp), dimension (:, :, :), allocatable :: blowim
  real(dp), dimension (:, :, :), allocatable :: bbnkim
  real(dp), dimension (:, :), allocatable :: sm12_real
  complex(dp), dimension (:, :, :), allocatable :: sm12_complex
  

  real(dp), dimension (:, :), allocatable :: eigen_k
  real(dp), dimension (:, :), allocatable :: special_k
  integer :: norbitals
  integer :: norbitals_new
  integer :: nkpoints  
  integer, dimension (:), allocatable :: getiatom
  integer, dimension (:,:),allocatable :: ioccupy_k !AQUI pensar allo en denmat
  integer, dimension (:), allocatable :: ioccupy   !AQUI
  real(dp), dimension (:,:), allocatable :: foccupy !AQUI
  
  !scf
  real(dp), dimension (:, :, :, :), allocatable :: cape
  real(dp), dimension (:, :, :, :), allocatable :: rhoPP
  real(dp) :: ztot
  real(dp), dimension (:), allocatable :: weight_k
  integer ::  nssh_tot
  real(dp), dimension (:), allocatable :: mwe
  real(dp), dimension (:), allocatable :: drwe

  ! -- EWALD--
  real(dp), dimension (:, :), allocatable :: ewald
  real(dp), dimension (:, :, :), allocatable :: dewald
  real(dp), dimension (:, :), allocatable :: fewald
  real(dp), dimension (:, :, :, :), allocatable :: ewaldsr
  real(dp), dimension (:, :, :, :), allocatable :: dip
  real(dp), dimension (:, :, :, :, :), allocatable :: dipp
 
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
  real(dp), dimension (:, :, :, :), allocatable :: sVNL
  real(dp), dimension (:, :, :, :, :), allocatable :: spVNL
  real(dp), dimension (:, :, :, :, :), allocatable :: sp_mat
  real(dp), dimension (:, :, :, :, :), allocatable :: tp_mat
  real(dp), dimension (:, :, :), allocatable :: dipcm
  real(dp), dimension (:, :, :, :), allocatable :: dippcm
  real(dp), dimension (:, :, :, :, :, :), allocatable :: dippc
  real(dp), dimension(3,3) :: eps2  
  real(dp), dimension (:, :, :, :), allocatable :: vnl
  integer, dimension (:), allocatable :: neighPP_comn
  integer, dimension (:, :), allocatable :: neighPP_comm
  integer, dimension (:, :, :), allocatable :: neighPP_comj
  integer, dimension (:, :, :), allocatable :: neighPP_comb 
  integer, dimension(:, :), allocatable   :: neighj_aux
  integer :: neighPP_max 
  integer :: num_neig_maxtot
 
  !CHARGES
  real(dp), dimension (:, :), allocatable :: Qin
  real(dp), dimension (:), allocatable :: Qinmixer
  real(dp), dimension (:, :), allocatable :: Qout
  real(dp), dimension (:), allocatable :: Qoutmixer
  real(dp), dimension (:), allocatable :: dq
  real(dp), dimension(:), allocatable  :: Q_partial
  real(dp), dimension (:, :), allocatable :: Qin_es
  real(dp), dimension (:), allocatable :: QLowdin_TOT
  real(dp), dimension (:), allocatable :: QMulliken_TOT
  real(dp), allocatable, dimension(:,:) :: qaux
  real(dp) :: dip_x, dipQout_x, dipTot_x, dipProy_x, dipIntra_x, dip_res_x, dipQin_x, dipRes_x
  real(dp) :: dip_y, dipQout_y, dipTot_y, dipProy_y, dipIntra_y, dip_res_y, dipQin_y, dipRes_y
  real(dp) :: dip_z, dipQout_z, dipTot_z, dipProy_z, dipIntra_z, dip_res_z, dipQin_z, dipRes_z
  real(dp) :: dip_tot, dip_proy, dipQin_tot, dipTot_tot, dipIntra_tot, dipQout_tot, dip_res_tot, dipRes_tot 
  real(dp), dimension (:), allocatable :: dq_DP

  !interaccions
  real(dp), dimension (:, :, :, :), allocatable :: vxc
  real(dp), dimension (:, :, :, :), allocatable :: vxc_ca
  real(dp), dimension (:, :, :, :), allocatable :: rho
  real(dp), dimension (:, :, :, :), allocatable :: rho_off
  real(dp), dimension (:, :, :, :), allocatable :: rhoij_off
  real(dp), dimension (:, :, :, :), allocatable :: s_mat 
  real(dp), dimension (:, :, :, :), allocatable :: sm_mat
  real(dp), dimension (:, :, :, :, :), allocatable :: spm_mat
  real(dp), dimension (:, :, :), allocatable :: rho_on
  real(dp), dimension (:, :, :), allocatable :: arho_on
  real(dp), dimension (:, :, :), allocatable :: rhoi_on
  real(dp), dimension (:, :, :), allocatable :: arhoi_on
  real(dp), dimension (:, :, :, :, :), allocatable :: arhop_on
  real(dp), dimension (:, :, :, :, :), allocatable :: rhop_on
  real(dp), dimension (:, :, :, :), allocatable :: arhoij_off
  real(dp), dimension (:, :, :, :), allocatable :: arho_off
  real(dp), dimension (:, :, :, :, :), allocatable :: arhopij_off 
  real(dp), dimension (:, :, :, :, :), allocatable :: arhop_off
  real(dp), dimension (:, :, :, :, :), allocatable :: rhop_off
  real(dp), dimension (:, :, :, :, :), allocatable :: rhopij_off
  real(dp), dimension (:, :, :, :), allocatable :: vca
  real(dp), dimension (:, :, :, :), allocatable :: ewaldlr
  real(dp), dimension (:, :, :, :), allocatable :: h_mat
  real(dp), dimension (:, :, :, :), allocatable :: h_mat0   ! h^0_{mu nu}: H a cargas neutras (Kscf=1, dQ=0)
  real(dp), dimension (:, :, :, :), allocatable :: t_mat
  real(dp), dimension (:, :, :, :), allocatable :: vna
  real(dp), dimension (:, :, :, :), allocatable :: ewaldqmmm
  real(dp), dimension (:, :, :, :, :), allocatable :: dipc
  integer, dimension(:,:), allocatable :: muR
  integer, dimension(:,:), allocatable :: nuR
  integer, dimension(:,:), allocatable :: alphaR
  integer, dimension(:,:), allocatable :: betaR
  real(dp), dimension(:,:,:,:), allocatable   :: hr_box

  
  real(dp), allocatable, dimension(:,:) :: Fv   
  real(dp), allocatable, dimension(:,:) :: Xv    
  real(dp), allocatable, dimension(:,:) :: delF 
  real(dp), allocatable, dimension(:,:) :: delX 
  real(dp), allocatable, dimension(:)   :: r2_sav
  real(dp), allocatable, dimension(:)   :: wi
  real(dp), allocatable, dimension(:)   :: x_best
 
  real(dp), allocatable, dimension(:,:) :: RJac
  real(dp), allocatable, dimension(:,:) :: betaInvH
  real(dp), allocatable, dimension(:,:) :: gamaH
  real(dp), dimension (:,:), allocatable :: xl

  ! Lattice vectors
  real(dp), dimension (3) :: a1vec
  real(dp), dimension (3) :: a2vec
  real(dp), dimension (3) :: a3vec

  ! allocate_f
  real(dp), dimension (:, :, :), allocatable :: fotnl
  real(dp), dimension (:, :, :), allocatable :: fanl
  real(dp), dimension (:, :, :), allocatable :: fotna
  real(dp), dimension (:, :, :), allocatable :: fana
  real(dp), dimension (:, :, :), allocatable :: faxc 
  real(dp), dimension (:, :, :), allocatable :: faxc_ca 
  real(dp), dimension (:, :, :), allocatable :: dxcdcc
  real(dp), dimension (:, :, :), allocatable :: dxcdcc_zw
  real(dp), dimension (:, :), allocatable :: ft
  real(dp), dimension (:, :), allocatable :: dusr
  real(dp), dimension (:, :, :), allocatable :: fotxc_ca
  real(dp), dimension (:, :, :), allocatable :: faca 
  real(dp), dimension (:, :, :), allocatable :: fotca
  real(dp), dimension (:, :), allocatable :: f3naa, f3nab, f3nac
  real(dp), dimension (:, :), allocatable :: f3nla, f3nlb, f3nlc
  real(dp), dimension (:, :), allocatable :: f3caa, f3cab, f3cac
  real(dp), dimension (:, :), allocatable :: flrew
  real(dp), dimension (:, :), allocatable :: f3xca_ca, f3xcb_ca, f3xcc_ca
  real(dp), dimension (:, :), allocatable :: f3xca, f3xcb, f3xcc
  real(dp), dimension (:, :), allocatable :: flrew_qmmm
  real(dp), dimension (:, :), allocatable :: fro
  real(dp), dimension (:, :), allocatable :: ftot
  real(dp), dimension (:, :), allocatable :: dxcv


  !integer                        :: qmmm_qm_mm_pairs     !Number of pairs per QM atom. - length of pair_list. 
  integer                        :: qmmm_qm_natoms
  real(dp), dimension(:,:), allocatable :: qmmm_dxyzcl  !Used to store the forces generated by qm_mm before adding them to the main f array.
  real(dp), dimension(:,:),  allocatable ::  qmmm_qm_xcrd         !Contains imaged mm coordinates and scaled mm charges.
  real(dp), dimension(:), allocatable :: qmmm_scf_mchg
  real(dp), dimension(:), allocatable :: qmmm_Qneutral_TOT
  real(dp) :: qmmm_rc1, qmmm_rc2, qmmm_width


  real(dp), dimension(3, 3, 5) :: amat
  
  real(dp), dimension (:,:,:,:,:), allocatable :: g_h
  real(dp), dimension (:,:,:,:,:), allocatable :: g_xc
  ! N^gamma_munu = <phi_mu|(phi^S_gamma)^2|phi_nu>, gamma: global shell index
  real(dp), dimension (:,:,:,:,:), allocatable :: den_or
  ! M^gamma_alphabeta = <phi^S_alpha|(phi^S_gamma)^2|phi^S_beta>, gamma: global shell index
  real(dp), dimension (:,:,:,:,:), allocatable :: den_sh
  ! Speed up SCF loop (seccion 4 de doc/XC_basic_claude_fable.tex): tablas Q-independientes
  ! adicionales rellenadas en Kscf=1 para reconstruir densidades/vxc sin interpolar en Kscf>1.
  ! iscf_fast=0 recupera el camino clasico (reensamblado con interpolaciones en cada paso).
  integer :: iscf_fast = 1
  ! pieza 2c (L+R) de den_or/den_sh: reconstruye rhoij_off/arhoij_off (gamma: indice global)
  real(dp), dimension (:,:,:,:,:), allocatable :: denij_or
  real(dp), dimension (:,:,:,:,:), allocatable :: denij_sh
  ! pieza self (atomo i sobre si mismo, indice de shell local): reconstruye rhoi_on/arhoi_on
  real(dp), dimension (:,:,:,:), allocatable :: deni_or
  real(dp), dimension (:,:,:,:), allocatable :: deni_sh
  ! termino constante VXC_0 (isorp=0) off-site de assemble_xc_2c, congelado en Kscf=1
  real(dp), dimension (:,:,:,:), allocatable :: vxc_ca0
  ! integrales de Coulomb 2c por par de vecinos (assemble_usr), congeladas en Kscf=1
  real(dp), dimension (:,:,:,:), allocatable :: coulomb_mat
  real(dp), dimension (:,:,:,:,:), allocatable :: f_xc
  real(dp), dimension (:,:,:,:), allocatable :: exc_aa
  real(dp), dimension (:,:,:,:), allocatable :: vxc_aa
  real(dp), dimension (:,:), allocatable :: g_h_shell
  real(dp), dimension (:,:), allocatable :: g_xc_shell  
  real(dp), dimension (:,:), allocatable :: f_xc_shell
  real(dp), dimension (:), allocatable :: exc_aa_shell
  real(dp), dimension (:), allocatable :: vxc_aa_shell
  integer, dimension (:), allocatable :: fix_shell_charge
  ! Seleccion de shells fijadas a carga neutra en stationary_charges (se fija desde Python
  ! con set_fix_shells): 0 = ninguna, 1 = solo shells d, 2 = shells con Qneutral=0 (default),
  ! 3 = mascara de usuario en fix_shell_user (indice de shell global: atomo1-shells, atomo2-...)
  integer :: ifix_shells = 2
  integer, dimension (:), allocatable :: fix_shell_user
  integer, dimension(:,:), allocatable :: orb2shell
end module M_system
