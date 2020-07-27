project('trans4d',
  'cpp',
  version: '0.0.1',
  license: 'MIT',
  default_options: ['cpp_std=c++11', 'warning_level=2'])

gtest_proj = subproject('gtest')
gtest_dep = gtest_proj.get_variable('gtest_main_dep')

inc = include_directories('include')

subdir('include')
subdir('src')
subdir('test')
subdir('doc')

pkg_mod = import('pkgconfig')
pkg_mod.generate(libraries: libtrans4d,
                 version: '1.0',
                 name: 'libtrans4d',
                 filebase: 'trans4d',
                 description: 'C++ Port of Fortran TRANS4D by Snay, Pearson, and Saleh')

