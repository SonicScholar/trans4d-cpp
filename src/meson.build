libtrans4d_depends = ['Trans4d.cpp', 'Trans4dCommon.cpp', 'InitBd.cpp', 'InitEq.cpp', 'InitPs.cpp', 'InitVl.cpp']
trans4d_depends = ['Trans4dExample.cpp']

libtrans4d = shared_library('trans4d',
  sources: libtrans4d_depends,
  include_directories: inc,
  install : true)

trans4d = executable('trans4d',
  sources: trans4d_depends,
  include_directories: inc,
  link_with: libtrans4d,
  install: true)
