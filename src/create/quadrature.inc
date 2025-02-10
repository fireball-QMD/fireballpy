! File quadrature.inc with all of the quadrature and grid information in it.
! =========================================================
! two-center parameters for computing matrix elements.
 
! overlap parameters.  (Must be accurate, because they determine othogonality)
        integer nzs
        parameter (nzs = 106)
 
        integer nrhos
        parameter (nrhos = 106)
 
        integer ndds
        parameter (ndds = 107)
 
! (non)-neutral potential atom/ontop parameters
        integer nznao
        parameter (nznao = 106)
 
        integer nrhonao
        parameter (nrhonao = 106)
 
        integer nddnao
        parameter (nddnao = 107)
 
! (non)-neutral potential atom/atom parameters
        integer nznat
        parameter (nznat = 106)
 
        integer nrhonat
        parameter (nrhonat = 106)
 
        integer nddnat
        parameter (nddnat = 107)
 
! non-local parameters.
        integer nznl
		parameter (nznl = 106)
	 
		integer nrhonl
		parameter (nrhonl = 106)
	 
		integer nddnl
		parameter (nddnl = 107)
	 
	! exchange-correlation potential (ontop) parameters
		integer nzxco
		parameter (nzxco = 106)
	 
		integer nrhoxco
		parameter (nrhoxco = 106)
	 
		integer nddxco
		parameter (nddxco = 107)
	 
	! exchange-correlation potential (atom) parameters
		integer nzxca
		parameter (nzxca = 106)
	 
		integer nrhoxca
		parameter (nrhoxca = 106)
	 
		integer nddxca
		parameter (nddxca = 107)
	 
	! exact exchange potential (ontop) parameters
		integer nzexo
		parameter (nzexo = 106)
	 
		integer nrhoexo
		parameter (nrhoexo = 106)
	 
		integer nddexo
		parameter (nddexo = 107)
	 
	! exact exchange potential (atom) parameters
		integer nzexa
		parameter (nzexa = 106)
	 
		integer nrhoexa
		parameter (nrhoexa = 106)
	 
		integer nddexa
		parameter (nddexa = 107)
	 
	! dipole parameters. Also used for coulomb interactions.
		integer nzd
		parameter (nzd = 106)
	 
		integer nrhod
		parameter (nrhod = 106)
	 
		integer nddd
		parameter (nddd = 107)
	 
	! Extended Hubbard parameters
		integer nzeh
		parameter (nzeh = 106)
	 
		integer nrhoeh
		parameter (nrhoeh = 106)
 
        integer nddeh
        parameter (nddeh = 107)
 
! When the density is calculated the grid size is dependent normally upon
! the grid size of the wavefunctions. This can be rather small most of the
! time. The factor ixcgridfactor allows the user to increase the grid size
! to something more manageable, making the code run faster.
        integer ixcgridfactor
        parameter (ixcgridfactor = 15)
 
! quadrature for theta grids in three-center interactions
        integer, parameter :: ntheta_max = 5      !  Gauss Legendre expansion
