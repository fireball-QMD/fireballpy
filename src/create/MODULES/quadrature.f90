        module quadrature
         use precision

! quadratures for kinetic energy
         integer, parameter :: nqke = 400
         integer, parameter :: nrke = 240
         integer, parameter :: nddke = 107
         real(kind=long), parameter :: ecutke = 40.0d3

! Three-center interactions
! ****************************************************************************
! quadrature for theta grids in three-center interactions
         integer, parameter :: ntheta_max = 5      !  Gauss Legendre expansion

! quadratures for the bond-charge and neutral atom distances
         integer, parameter :: nbcba = 29
         integer, parameter :: nnaba = 29

! exchange-correlation: phi, r, and theta grids
         integer, parameter :: nphi2ba_xc = 15
         integer, parameter :: nrba_xc = 31
         integer, parameter :: nthetaba_xc = 15  

! neutral atom: phi, r, and theta grids
         integer, parameter :: nphi2ba_na = 15
         integer, parameter :: nrba_na = 31
         integer, parameter :: nthetaba_na = 15 
! ****************************************************************************
        end module

! overlap parameters.
!       integer nzs
!       parameter (nzs = 106)
!
!       integer nrhos
!       parameter (nrhos = 106)
!
!       integer ndds
!       parameter (ndds = 107)
!
! (non)-neutral potential atom/ontop parameters
!       integer nznao
!       parameter (nznao = 106)
!
!       integer nrhonao
!       parameter (nrhonao = 106)
!
!       integer nddnao
!       parameter (nddnao = 107)
!
! (non)-neutral potential atom/atom parameters
!       integer nznat
!       parameter (nznat = 106)
!
!       integer nrhonat
!       parameter (nrhonat = 106)
!
!       integer nddnat
!       parameter (nddnat = 107)
!
! non-local parameters.
!       integer nznl
!       parameter (nznl = 146)
!
!       integer nrhonl
!       parameter (nrhonl = 146)
!
!       integer nddnl
!       parameter (nddnl = 107)
!
! exchange-correlation potential (ontop) parameters
!       integer nzxco
!       parameter (nzxco = 106)
!
!       integer nrhoxco
!       parameter (nrhoxco = 106)
!
!       integer nddxco
!       parameter (nddxco = 107)
!
! exchange-correlation potential (atom) parameters
!       integer nzxca
!       parameter (nzxca = 106)
!
!       integer nrhoxca
!       parameter (nrhoxca = 106)
!
!       integer nddxca
!       parameter (nddxca = 107)
!
! exact exchange potential (ontop) parameters
!       integer nzexo
!       parameter (nzexo = 106)
!
!       integer nrhoexo
!       parameter (nrhoexo = 106)
!
!       integer nddexo
!       parameter (nddexo = 107)
!
! exact exchange potential (atom) parameters
!       integer nzexa
!       parameter (nzexa = 106)
!
!       integer nrhoexa
!       parameter (nrhoexa = 106)
!
!       integer nddexa
!       parameter (nddexa = 107)
!
! dipole parameters. Also used for coulomb interactions.
!       integer nzd
!       parameter (nzd = 156)
!
!       integer nrhod
!       parameter (nrhod = 156)
!
!       integer nddd
!       parameter (nddd = 107)
!
! Extended Hubbard parameters
!       integer nzeh
!       parameter (nzeh = 106)
!
!       integer nrhoeh
!       parameter (nrhoeh = 106)
!
!       integer nddeh
!       parameter (nddeh = 107)
!
! When the density is calculated the grid size is dependent normally upon
! the grid size of the wavefunctions. This can be rather small most of the
! time. The factor ixcgridfactor allows the user to increase the grid size
! to something more manageable, making the code run faster.
!       integer ixcgridfactor
!       parameter (ixcgridfactor = 15)
! ======================= Last line of quadrature.inc ========
