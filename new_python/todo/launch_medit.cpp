#include <npe.h>
#include <typedefs.h>
#include <igl/launch_medit.h>

const char* ds_launch_medit = R"igl_Qu8mg5v7(

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

 Writes the tetmesh in (V,T,F) to a temporary file, opens it with medit
   (forking with a system call) and returns
  
  
   Templates:
     DerivedV  real-value: i.e. from MatrixXd
     DerivedT  integer-value: i.e. from MatrixXi
     DerivedF  integer-value: i.e. from MatrixXi
   Inputs:
     V  double matrix of vertex positions  #V by 3
     T  #T list of tet indices into vertex positions
     F  #F list of face indices into vertex positions
     wait  whether to wait for medit process to finish before returning
   Returns returned value of system call (probably not useful if wait=false
   because of the fork)
)igl_Qu8mg5v7";

npe_function(launch_medit)
npe_doc(ds_launch_medit)

npe_arg(v, dense_f32, dense_f64)
npe_arg(t, dense_f32, dense_f64)
npe_arg(f, dense_i32, dense_i64)
npe_arg(wait, bool)


npe_begin_code()

  igl::launch_medit(v, t, f, wait);
  return ;

npe_end_code()


