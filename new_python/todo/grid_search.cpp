#include <npe.h>
#include <typedefs.h>






#include <igl/grid_search.h>

const char* ds_grid_search = R"igl_Qu8mg5v7(

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

 Solve the problem:
  
     minimize f(x)
     subject to lb ≤ x ≤ ub 
   
   by exhaustive grid search.
  
   Inputs:
     f  function to minimize
     LB  #X vector of finite lower bounds
     UB  #X vector of finite upper bounds
     I  #X vector of number of steps for each variable
   Outputs:
     X  #X optimal parameter vector
   Returns f(X)
  
)igl_Qu8mg5v7";

npe_function(grid_search)
npe_doc(ds_grid_search)

npe_arg(f, std::function<Scalar (DerivedX &)>)
npe_arg(lb, dense_f32, dense_f64)
npe_arg(ub, dense_f32, dense_f64)
npe_arg(i, dense_f32, dense_f64)


npe_begin_code()

  DerivedX & x;
  igl::grid_search(f, lb, ub, i, x);
  return npe::move(x);

npe_end_code()


