# copyright info:
#
#                             @Copyright 2003
#                           Fireball Committee
# Brigham Young University - James P. Lewis, Chair
# Arizona State University - Otto F. Sankey
# Universidad de Madrid - Jose Ortega
# Lawrence Livermore National Laboratory - Kurt Glaesemann
# Universidad de Madrid - Pavel Jelenik
#
# Other contributors, past and present:
# Auburn University - Jian Jun Dong
# Arizona State University - Gary B. Adams
# Arizona State University - Kevin Schmidt
# Arizona State University - John Tomfohr
# Motorola, Physical Sciences Research Labs - Alex Demkov
# Motorola, Physical Sciences Research Labs - Jun Wang
# Ohio University - Dave Drabold
# University of Regensburg - Juergen Fritsch

# RESTRICTED RIGHTS LEGEND
#
# Use, duplication, or disclosure of this software and its documentation
# by the Government is subject to restrictions as set forth in
# subdivsion { (b) (3) (ii) } of the Rights in Technical Data and
# Computer Software clause at 52.227-7013.

#Pick your machine type: AIX, AIX_CHPC, ICEBOX, SGI
MACHINE = ILINUX
#Pick you parallel type: <blank for none>, MPI, MPICH
METHOD = DEBUG
#METHOD = OPT
#Define debug to turn on symbols and array checking: <blank for none>, TRUE
#DEBUG =      
USEBLAS = YES
#Set bogus default to check for unkown machine types
F77 = ERROR_IN_THE_MAKEFILE
METHOD = DOUBLE
#METHOD = GAMMA 
#METHOD = SCALAPACK
# Do you want to do thermodynamic integration (NO,YES)?  If so 
# # you'll need a compatible C compiler.
CCOMPILE = NO
# # Do you use TclMD option? (YES/NO)
# # Note: you must have TclMD library installed and Tcl/Tk to use this.
# # You can get Tcl/Tk from http://dev.scriptics.com/software/tcltk/
TCLMD = NO
#
# # Choose a profiler - none, VAMPIR
PROFILER = -pg
#


#Do we debug?
ifneq (,$(findstring TRUE, ${DEBUG}))
	OPTIMIZATION = -C -g -check all -no-ipo
endif

#Do we do MPI?
MPIFILE = MPI_stubs.o
ifneq (,$(findstring MPI, ${METHOD}))
	MPIFILE = MPI_stuff.o
endif


#include MACHINES/AIX
#include MACHINES/HUDSON
#include MACHINES/LINUX
#include MACHINES/ILINUX
#include MACHINES/ILINUX_beem
#include MACHINES/MARYLOUX
#include MACHINES/SCRI
#include MACHINES/SGI
#include MACHINES/TRU64

include MACHINES/${MACHINE}

#FFLAGS =  $(PARALLEL) $(OPTIMIZATION)

MODULES = precision.o coefficients.o constants.o dimensions.o quadrature.o \
	  x_exact.o

EXX = clebsch_gordon.o factorial.o setup_clm.o x_1c.o xatom_2c.o \
	xatom_2c_integral.o xatom_2c_rprime.o xontopl_2c.o \
	xontopl_2c_integral.o xontopl_2c_rprime.o xontopr_2c.o \
	xontopr_2c_integral.o xontopr_2c_rprime.o

INTERPOLATERS = interpolate2d.o psiofr.o rprimeofr.o vnnaofr.o vppofr.o buildspline_1d.o \
		dpsiofr.o

INTEGRALS = onecenterxc.o onecenternuxc.o onecenternuxcs.o onecentervdip.o radial_integral.o append_string.o \
	gaunt.o gauntReal.o mk_index.o \
	mk_indexPP.o mk_indexDipX.o mk_indexDipY.o rescaled_psi.o threecenter.o \
	gleg.o  threecenter_integral.o twocenter.o twocenter_integral.o \
	exc1crho.o nuxc1crho.o goverlap1c.o

KINETIC = kinetic.o Pintegral.o Tintegral.o bessel.o tree.o

READFILES = iofile2c.o iofile2c_x.o iofile2c_xcna.o iofile3c.o iofile.o \
	readpsi.o readtheory.o readvnn.o readvpp.o readcreate.o

XC = bpgcxc.o ceperley_alder.o corlyp1c.o corlyp2c.o density_calc.o dexc.o \
	dpotxc12.o dpotxc12s.o dvxc.o dvxc3c.o get_potxc1c.o get_potxc2c.o \
	gga91_sr.o ggacrad1c.o ggaxrad1c.o ggaxrad2c.o ggacrad2c.o ldaxc.o \
	lsdavwn.o pbe.o rho1c_store.o rho2c_store.o vxc.o

GAUSS = gauss_create.o gauss_subs.o gauss_density.o

MAIN = create.o 

OBJECTS = $(EXX) $(INTERPOLATERS) $(INTEGRALS) $(KINETIC) $(READFILES) $(XC) \
	$(MAIN) $(MPIFILE) $(GAUSS)

create.x: $(MODULES) $(OBJECTS) 
	$(F90) $(FFLAGS) -o $@ $(MODULES) $(FFLAGS) $(OBJECTS) $(LFLAGS) $(VISFLAGS) $(PARLFLAGS)

create.o: create.f parameters.inc exchange.inc quadrature.inc \
	wavefunctions.inc
	$(F77) $(FFLAGS) -c create.f

.PHONY: clean veryclean

clean:
	rm -f -r core *.o .nfs* rii_files create.x.ipa* .#* *.mod ldtmp*

veryclean: clean
	rm -f create.x

# ****************************************************************************
# Modules
# ****************************************************************************
precision.o : MODULES/precision.f90
	$(F90) $(FFLAGS) -c MODULES/precision.f90
coefficients.o : MODULES/coefficients.f90
	$(F90) $(FFLAGS) -c MODULES/coefficients.f90
constants.o : MODULES/constants.f90
	$(F90) $(FFLAGS) -c MODULES/constants.f90
dimensions.o : MODULES/dimensions.f90
	$(F90) $(FFLAGS) -c MODULES/dimensions.f90
quadrature.o : MODULES/quadrature.f90
	$(F90) $(FFLAGS) -c MODULES/quadrature.f90
x_exact.o : MODULES/x_exact.f90
	$(F90) $(FFLAGS) -c MODULES/x_exact.f90


# ****************************************************************************
# exact exchange objects
# ****************************************************************************
clebsch_gordon.o : EXX/clebsch_gordon.f90 $(MODULES)
	$(F90) $(FFLAGS) -c EXX/clebsch_gordon.f90
factorial.o : EXX/factorial.f90 $(MODULES)
	$(F90) $(FFLAGS) -c EXX/factorial.f90
setup_clm.o : EXX/setup_clm.f90 $(MODULES)
	$(F90) $(FFLAGS) -c EXX/setup_clm.f90
x_1c.o : EXX/x_1c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c EXX/x_1c.f90
xatom_2c.o : EXX/xatom_2c.f90 quadrature.inc $(MODULES)
	$(F90) $(FFLAGS) -c EXX/xatom_2c.f90
xatom_2c_integral.o : EXX/xatom_2c_integral.f90 $(MODULES)
	$(F90) $(FFLAGS) -c EXX/xatom_2c_integral.f90
xatom_2c_rprime.o : EXX/xatom_2c_rprime.f90 $(MODULES)
	$(F90) $(FFLAGS) -c EXX/xatom_2c_rprime.f90
xontopl_2c.o : EXX/xontopl_2c.f90 quadrature.inc $(MODULES)
	$(F90) $(FFLAGS) -c EXX/xontopl_2c.f90
xontopl_2c_integral.o : EXX/xontopl_2c_integral.f90 $(MODULES)
	$(F90) $(FFLAGS) -c EXX/xontopl_2c_integral.f90
xontopl_2c_rprime.o : EXX/xontopl_2c_rprime.f90 $(MODULES)
	$(F90) $(FFLAGS) -c EXX/xontopl_2c_rprime.f90
xontopr_2c.o : EXX/xontopr_2c.f90 quadrature.inc $(MODULES)
	$(F90) $(FFLAGS) -c EXX/xontopr_2c.f90
xontopr_2c_integral.o : EXX/xontopr_2c_integral.f90 $(MODULES)
	$(F90) $(FFLAGS) -c EXX/xontopr_2c_integral.f90
xontopr_2c_rprime.o : EXX/xontopr_2c_rprime.f90 $(MODULES)
	$(F90) $(FFLAGS) -c EXX/xontopr_2c_rprime.f90


# ****************************************************************************
# exchange-correlations objects
# ****************************************************************************
bpgcxc.o : XC/bpgcxc.f $(MODULES)
	$(F77) $(FFLAGS) -c XC/bpgcxc.f
ceperley_alder.o : XC/ceperley_alder.f90 $(MODULES)
	$(F90) $(FFLAGS) -c XC/ceperley_alder.f90
corlyp1c.o : XC/corlyp1c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c XC/corlyp1c.f90
corlyp2c.o : XC/corlyp2c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c XC/corlyp2c.f90
density_calc.o : XC/density_calc.f parameters.inc exchange.inc $(MODULES)
	$(F77) $(FFLAGS) -c XC/density_calc.f
dexc.o : XC/dexc.f parameters.inc wavefunctions.inc $(MODULES)
	$(F77) $(FFLAGS) -c XC/dexc.f
dpotxc12.o : XC/dpotxc12.f parameters.inc wavefunctions.inc $(MODULES)
	$(F77) $(FFLAGS) -c XC/dpotxc12.f
dpotxc12s.o : XC/dpotxc12s.f $(MODULES)
	$(F77) $(FFLAGS) -c XC/dpotxc12s.f
dvxc.o : XC/dvxc.f parameters.inc wavefunctions.inc $(MODULES)
	$(F77) $(FFLAGS) -c XC/dvxc.f
dvxc3c.o : XC/dvxc3c.f parameters.inc wavefunctions.inc $(MODULES)
	$(F77) $(FFLAGS) -c XC/dvxc3c.f
get_potxc1c.o : XC/get_potxc1c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c XC/get_potxc1c.f90
get_potxc2c.o : XC/get_potxc2c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c XC/get_potxc2c.f90
gga91_sr.o : XC/gga91_sr.f $(MODULES)
	$(F77) $(FFLAGS) -c XC/gga91_sr.f
ggacrad1c.o : XC/ggacrad1c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c XC/ggacrad1c.f90
ggaxrad1c.o : XC/ggaxrad1c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c XC/ggaxrad1c.f90
ggaxrad2c.o : XC/ggaxrad2c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c XC/ggaxrad2c.f90
ggacrad2c.o : XC/ggacrad2c.f $(MODULES)
	$(F77) $(FFLAGS) -c XC/ggacrad2c.f
ldaxc.o : XC/ldaxc.f $(MODULES)
	$(F77) $(FFLAGS) -c XC/ldaxc.f
lsdavwn.o : XC/lsdavwn.f90 $(MODULES)
	$(F90) $(FFLAGS) -c XC/lsdavwn.f90
pbe.o : XC/pbe.f $(MODULES)
	$(F77) $(FFLAGS) -c XC/pbe.f
rho1c_store.o : XC/rho1c_store.f $(MODULES)
	$(F77) $(FFLAGS) -c XC/rho1c_store.f
rho2c_store.o : XC/rho2c_store.f parameters.inc exchange.inc $(MODULES) \
	quadrature.inc wavefunctions.inc
	$(F77) $(FFLAGS) -c XC/rho2c_store.f
vxc.o : XC/vxc.f parameters.inc wavefunctions.inc $(MODULES)
	$(F77) $(FFLAGS) -c XC/vxc.f


# ****************************************************************************
# integrals objects
# ****************************************************************************
onecenterxc.o : INTEGRALS/onecenterxc.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTEGRALS/onecenterxc.f90
gauntReal.o : INTEGRALS/gauntReal.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTEGRALS/gauntReal.f90
gaunt.o : INTEGRALS/gaunt.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTEGRALS/gaunt.f90
onecentervdip.o : INTEGRALS/onecentervdip.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTEGRALS/onecentervdip.f90
radial_integral.o : INTEGRALS/radial_integral.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTEGRALS/radial_integral.f90
append_string.o : INTEGRALS/append_string.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTEGRALS/append_string.f90
onecenternuxc.o : INTEGRALS/onecenternuxc.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTEGRALS/onecenternuxc.f90
onecenternuxcs.o : INTEGRALS/onecenternuxcs.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTEGRALS/onecenternuxcs.f90
mk_index.o : INTEGRALS/mk_index.f $(MODULES)
	$(F77) $(FFLAGS) -c INTEGRALS/mk_index.f
mk_indexPP.o : INTEGRALS/mk_indexPP.f $(MODULES)
	$(F77) $(FFLAGS) -c INTEGRALS/mk_indexPP.f
mk_indexDipX.o : INTEGRALS/mk_indexDipX.f $(MODULES)
	$(F77) $(FFLAGS) -c INTEGRALS/mk_indexDipX.f
mk_indexDipY.o : INTEGRALS/mk_indexDipY.f $(MODULES)
	$(F77) $(FFLAGS) -c INTEGRALS/mk_indexDipY.f
rescaled_psi.o : INTEGRALS/rescaled_psi.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTEGRALS/rescaled_psi.f90
threecenter.o : INTEGRALS/threecenter.f90 parameters.inc quadrature.inc $(MODULES)
	$(F90) $(FFLAGS) -c INTEGRALS/threecenter.f90
threecenter_integral.o : INTEGRALS/threecenter_integral.f parameters.inc $(MODULES)
	$(F77) $(FFLAGS) -c INTEGRALS/threecenter_integral.f
gleg.o : INTEGRALS/gleg.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTEGRALS/gleg.f90
twocenter.o : INTEGRALS/twocenter.f parameters.inc pseudopotentials.inc \
	quadrature.inc $(MODULES)
	$(F77) $(FFLAGS) -c INTEGRALS/twocenter.f
twocenter_integral.o : INTEGRALS/twocenter_integral.f $(MODULES)
	$(F77) $(FFLAGS) -c INTEGRALS/twocenter_integral.f
exc1crho.o : INTEGRALS/exc1crho.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTEGRALS/exc1crho.f90
nuxc1crho.o : INTEGRALS/nuxc1crho.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTEGRALS/nuxc1crho.f90
goverlap1c.o : INTEGRALS/goverlap1c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTEGRALS/goverlap1c.f90

# ****************************************************************************
# gaussian fits of psi, 3C-NA, and 3C-XC 
# ****************************************************************************
gauss_create.o : GAUSS/gauss_create.f $(MODULES)
	$(F77) $(FFLAGS) -c GAUSS/gauss_create.f

gauss_subs.o : GAUSS/gauss_subs.f $(MODULES)
	$(F77) $(FFLAGS) -c GAUSS/gauss_subs.f

gauss_density.o : GAUSS/gauss_density.f $(MODULES)
	$(F77) $(FFLAGS) -c GAUSS/gauss_density.f


# ****************************************************************************
# interpolaters objects
# ****************************************************************************
interpolate2d.o : INTERPOLATERS/interpolate2d.f $(MODULES)
	$(F77) $(FFLAGS) -c INTERPOLATERS/interpolate2d.f
psiofr.o : INTERPOLATERS/psiofr.f parameters.inc wavefunctions.inc $(MODULES)
	$(F77) $(FFLAGS) -c INTERPOLATERS/psiofr.f
dpsiofr.o : INTERPOLATERS/dpsiofr.f parameters.inc wavefunctions.inc $(MODULES)
	$(F77) $(FFLAGS) -c INTERPOLATERS/dpsiofr.f
rprimeofr.o : INTERPOLATERS/rprimeofr.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERPOLATERS/rprimeofr.f90
vnnaofr.o : INTERPOLATERS/vnnaofr.f parameters.inc vnonneutral.inc \
	wavefunctions.inc $(MODULES)
	$(F77) $(FFLAGS) -c INTERPOLATERS/vnnaofr.f
vppofr.o : INTERPOLATERS/vppofr.f parameters.inc pseudopotentials.inc $(MODULES)
	$(F77) $(FFLAGS) -c INTERPOLATERS/vppofr.f
buildspline_1d.o : INTERPOLATERS/buildspline_1d.f90
	$(F90) $(FFLAGS) -c INTERPOLATERS/buildspline_1d.f90


# ****************************************************************************
# kinetic objects
# ****************************************************************************
bessel.o : KINETIC/bessel.f $(MODULES)
	$(F77) $(FFLAGS) -c KINETIC/bessel.f
kinetic.o : KINETIC/kinetic.f90 parameters.inc $(MODULES)
	$(F90) $(FFLAGS) -c KINETIC/kinetic.f90
tree.o : KINETIC/tree.f $(MODULES)
	$(F77) $(FFLAGS) -c KINETIC/tree.f
Pintegral.o : KINETIC/Pintegral.f $(MODULES)
	$(F77) $(FFLAGS) -c KINETIC/Pintegral.f
Tintegral.o : KINETIC/Tintegral.f $(MODULES)
	$(F77) $(FFLAGS) -c KINETIC/Tintegral.f


# ****************************************************************************
# readfiles objects
# ****************************************************************************
iofile2c.o : READFILES/iofile2c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/iofile2c.f90
iofile2c_x.o : READFILES/iofile2c_x.f $(MODULES)
	$(F77) $(FFLAGS) -c READFILES/iofile2c_x.f
iofile2c_xcna.o : READFILES/iofile2c_xcna.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/iofile2c_xcna.f90
iofile3c.o : READFILES/iofile3c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/iofile3c.f90
iofile.o : READFILES/iofile.f $(MODULES)
	$(F77) $(FFLAGS) -c READFILES/iofile.f
readpsi.o : READFILES/readpsi.f parameters.inc wavefunctions.inc $(MODULES)
	$(F77) $(FFLAGS) -c READFILES/readpsi.f
readtheory.o : READFILES/readtheory.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/readtheory.f90
readvnn.o : READFILES/readvnn.f parameters.inc vnonneutral.inc $(MODULES)
	$(F77) $(FFLAGS) -c READFILES/readvnn.f
readvpp.o : READFILES/readvpp.f parameters.inc pseudopotentials.inc $(MODULES)
	$(F77) $(FFLAGS) -c READFILES/readvpp.f
readcreate.o : READFILES/readcreate.f parameters.inc exchange.inc $(MODULES)
	$(F77) $(FFLAGS) -c READFILES/readcreate.f

# ****************************************************************************
# MPI files
# ****************************************************************************
MPI_stubs.o : MPI_stubs.f90 $(MODULES)
	$(F90) $(FFLAGS) -c MPI_stubs.f90
MPI_stuff.o : MPI_stuff.f90 $(MODULES)
	$(F90) $(FFLAGS) -c MPI_stuff.f90
