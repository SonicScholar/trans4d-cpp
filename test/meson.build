test_dep = dependency('gtest', main : true, required : false)
test_trans4d = executable('trans4d',
               'Trans4d_test.cpp',
               include_directories: inc,
               dependencies : [ test_dep, gtest_dep])
test('gtest test', test_trans4d)
