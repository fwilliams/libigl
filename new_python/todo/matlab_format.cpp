#include <npe.h>
#include <typedefs.h>






#include <igl/matlab_format.h>

const char* ds_matlab_format = R"igl_Qu8mg5v7(

Parameters
----------


Returns
-------


See also
--------


Notes
-----
None

Examples
--------

 This is a routine to print a matrix using format suitable for pasting into
   the matlab IDE
  
   Templates:
     DerivedM  e.g. derived from MatrixXd
   Input:
     input  some matrix to be formatted
     name  name of matrix
   Returns  Formatted matrix
  
   Example:
    M := [1 2 3;4 5 6];
   cout<<matlab_format(M)<<endl;
    Prints:
    [
      1 2 3
      4 5 6
    ];
   cout<<matlab_format(M,"M")<<endl;
    Prints:
    M = [
      1 2 3
      4 5 6
    ];
)igl_Qu8mg5v7";

npe_function(matlab_format)
npe_doc(ds_matlab_format)

npe_arg(m, Eigen::DenseBase<DerivedM> &)
npe_arg(name, std::string)


npe_begin_code()

  igl::matlab_format(m, name);
  return ;

npe_end_code()
#include <igl/matlab_format.h>

const char* ds_matlab_format = R"igl_Qu8mg5v7(

Parameters
----------


Returns
-------


See also
--------


Notes
-----
None

Examples
--------

 Same but for sparse matrices. Print IJV format into an auxiliary variable
   and then print a call to sparse which will construct the sparse matrix
   Example:
    S := [0 2 3;4 5 0];
   cout<<matlab_format(S,"S")<<endl;
    Prints:
    SIJV = [
      2 1 4
      1 2 2
      2 2 5
      1 3 3
    ];
    S = sparse(SIJV(:,1),SIJV(:,2),SIJV(:,3));
  
)igl_Qu8mg5v7";

npe_function(matlab_format)
npe_doc(ds_matlab_format)

npe_arg(s, Eigen::SparseMatrix<DerivedS> &)
npe_arg(name, std::string)


npe_begin_code()

  igl::matlab_format(s, name);
  return ;

npe_end_code()
#include <igl/matlab_format.h>

const char* ds_matlab_format = R"igl_Qu8mg5v7(
See matlab_format for the documentation.
)igl_Qu8mg5v7";

npe_function(matlab_format)
npe_doc(ds_matlab_format)

npe_arg(v, double)
npe_arg(name, std::string)


npe_begin_code()

  igl::matlab_format(v, name);
  return ;

npe_end_code()
#include <igl/matlab_format.h>

const char* ds_matlab_format = R"igl_Qu8mg5v7(
See matlab_format for the documentation.
)igl_Qu8mg5v7";

npe_function(matlab_format)
npe_doc(ds_matlab_format)

npe_arg(v, float)
npe_arg(name, std::string)


npe_begin_code()

  igl::matlab_format(v, name);
  return ;

npe_end_code()


