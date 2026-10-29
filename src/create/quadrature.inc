! File quadrature.inc with the quadrature and grid information in it.
!
! Lo demas que habia aqui (nzs, nrhos, ndds, nznao ... nddeh, ixcgridfactor:
! 31 parametros) no lo usaba nadie: las mallas de dos centros las fijan ahora
! los modulos de INTEGRALS. Eliminado el 2026-10-22.

! quadrature for theta grids in three-center interactions
        integer, parameter :: ntheta_max = 5      !  Gauss Legendre expansion
