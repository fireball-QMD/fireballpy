modules_begin_srcs = [
  'begin_input.f90',
  'constants.f90',
  'pp_storage.f90',
  'precision.f90',
]

xcb_srcs = [
  'bpgcxc.f90',
  'corlyp.f90',
  'Exx.f90',
  'get_uxc.f90',
  'gga91_sr.f90',
  'ggacrad.f90',
  'ggaxrad.f90',
  'ldaxc.f90',
  'lsdavwn.f90',
  'pbe.f90',
]

rcatms_srcs = [
  'clebsch_gordon.f90',
  'derf0.f90',
  'get_nlmesh_size.f90',
  'get_psi.f90',
  'get_vee.f90',
  'integrate_hpsi.f90',
  'pp.f90',
  'psirc.f90',
  'rcatms.f90',
  'rcatms_DMOL.f90',
  'rcatms_excite.f90',
  'rcatms_GSM.f90',
  'vshort.f90',
]

vnn_srcs= [
  'get_ewfmesh_size.f90',
  'get_wfmesh_size.f90',
  'psiofr.f90',
  'readpsi.f90',
  'readpsi_excite.f90',
  'vnn.f90',
  'vnn_excite.f90',
]

begin_srcs = []
foreach src : modules_begin_srcs
  begin_srcs += 'MODULES' / src
endforeach
foreach src : xcb_srcs
  begin_srcs += 'XC' / src
endforeach
foreach src : rcatms_srcs
  begin_srcs += 'RCATMS' / src
endforeach
foreach src : vnn_srcs
  begin_srcs += 'VNN' / src
endforeach

# Build begin library
begin_lib = static_library('begin', begin_srcs)
begin_dep = declare_dependency(link_with: begin_lib, sources: begin_srcs)
