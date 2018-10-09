#include <npe.h>
#include <typedefs.h>






#include <igl/parallel_for.h>

const char* ds_parallel_for = R"igl_Qu8mg5v7(

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

 PARALLEL_FOR Functional implementation of a basic, open-mp style, parallel
   for loop. If the inner block of a for-loop can be rewritten/encapsulated in
   a single (anonymous/lambda) function call `func` so that the serial code
   looks like:
   
       for(int i = 0;i<loop_size;i++)
       {
         func(i);
       }
  
   then `parallel_for(loop_size,func,min_parallel)` will use as many threads as
   available on the current hardware to parallelize this for loop so long as
   loop_size<min_parallel, otherwise it will just use a serial for loop.
  
   Inputs:
     loop_size  number of iterations. I.e. for(int i = 0;i<loop_size;i++) ...
     func  function handle taking iteration index as only argument to compute
       inner block of for loop I.e. for(int i ...){ func(i); }
     min_parallel  min size of loop_size such that parallel (non-serial)
       thread pooling should be attempted {0}
   Returns true iff thread pool was invoked
)igl_Qu8mg5v7";

npe_function(parallel_for)
npe_doc(ds_parallel_for)

npe_arg(loop_size, Index)
npe_arg(func, FunctionType &)
npe_arg(min_parallel, size_t)


npe_begin_code()

  igl::parallel_for(loop_size, func, min_parallel);
  return ;

npe_end_code()
#include <igl/parallel_for.h>

const char* ds_parallel_for = R"igl_Qu8mg5v7(

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

 PARALLEL_FOR Functional implementation of an open-mp style, parallel for
   loop with accumulation. For example, serial code separated into n chunks
   (each to be parallelized with a thread) might look like:
       
       Eigen::VectorXd S;
       const auto & prep_func = [&S](int n){ S = Eigen:VectorXd::Zero(n); };
       const auto & func = [&X,&S](int i, int t){ S(t) += X(i); };
       const auto & accum_func = [&S,&sum](int t){ sum += S(t); };
       prep_func(n);
       for(int i = 0;i<loop_size;i++)
       {
         func(i,i%n);
       }
       double sum = 0;
       for(int t = 0;t<n;t++)
       {
         accum_func(t);
       }
   
   Inputs:
     loop_size  number of iterations. I.e. for(int i = 0;i<loop_size;i++) ...
     prep_func function handle taking n >= number of threads as only
       argument 
     func  function handle taking iteration index i and thread id t as only
       arguments to compute inner block of for loop I.e. 
       for(int i ...){ func(i,t); }
     accum_func  function handle taking thread index as only argument, to be
       called after all calls of func, e.g., for serial accumulation across
       all n (potential) threads, see n in description of prep_func.
     min_parallel  min size of loop_size such that parallel (non-serial)
       thread pooling should be attempted {0}
   Returns true iff thread pool was invoked
)igl_Qu8mg5v7";

npe_function(parallel_for)
npe_doc(ds_parallel_for)

npe_arg(loop_size, Index)
npe_arg(prep_func, PrepFunctionType &)
npe_arg(func, FunctionType &)
npe_arg(accum_func, AccumFunctionType &)
npe_arg(min_parallel, size_t)


npe_begin_code()

  igl::parallel_for(loop_size, prep_func, func, accum_func, min_parallel);
  return ;

npe_end_code()
#include <igl/parallel_for.h>

const char* ds_parallel_for = R"igl_Qu8mg5v7(

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

 PARALLEL_FOR Functional implementation of a basic, open-mp style, parallel
   for loop. If the inner block of a for-loop can be rewritten/encapsulated in
   a single (anonymous/lambda) function call `func` so that the serial code
   looks like:
   
       for(int i = 0;i<loop_size;i++)
       {
         func(i);
       }
  
   then `parallel_for(loop_size,func,min_parallel)` will use as many threads as
   available on the current hardware to parallelize this for loop so long as
   loop_size<min_parallel, otherwise it will just use a serial for loop.
  
   Inputs:
     loop_size  number of iterations. I.e. for(int i = 0;i<loop_size;i++) ...
     func  function handle taking iteration index as only argument to compute
       inner block of for loop I.e. for(int i ...){ func(i); }
     min_parallel  min size of loop_size such that parallel (non-serial)
       thread pooling should be attempted {0}
   Returns true iff thread pool was invoked
)igl_Qu8mg5v7";

npe_function(parallel_for)
npe_doc(ds_parallel_for)

npe_arg(loop_size, Index)
npe_arg(func, FunctionType &)
npe_arg(min_parallel, size_t)


npe_begin_code()

  igl::parallel_for(loop_size, func, min_parallel);
  return ;

npe_end_code()
#include <igl/parallel_for.h>

const char* ds_parallel_for = R"igl_Qu8mg5v7(

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

 PARALLEL_FOR Functional implementation of an open-mp style, parallel for
   loop with accumulation. For example, serial code separated into n chunks
   (each to be parallelized with a thread) might look like:
       
       Eigen::VectorXd S;
       const auto & prep_func = [&S](int n){ S = Eigen:VectorXd::Zero(n); };
       const auto & func = [&X,&S](int i, int t){ S(t) += X(i); };
       const auto & accum_func = [&S,&sum](int t){ sum += S(t); };
       prep_func(n);
       for(int i = 0;i<loop_size;i++)
       {
         func(i,i%n);
       }
       double sum = 0;
       for(int t = 0;t<n;t++)
       {
         accum_func(t);
       }
   
   Inputs:
     loop_size  number of iterations. I.e. for(int i = 0;i<loop_size;i++) ...
     prep_func function handle taking n >= number of threads as only
       argument 
     func  function handle taking iteration index i and thread id t as only
       arguments to compute inner block of for loop I.e. 
       for(int i ...){ func(i,t); }
     accum_func  function handle taking thread index as only argument, to be
       called after all calls of func, e.g., for serial accumulation across
       all n (potential) threads, see n in description of prep_func.
     min_parallel  min size of loop_size such that parallel (non-serial)
       thread pooling should be attempted {0}
   Returns true iff thread pool was invoked
)igl_Qu8mg5v7";

npe_function(parallel_for)
npe_doc(ds_parallel_for)

npe_arg(loop_size, Index)
npe_arg(prep_func, PreFunctionType &)
npe_arg(func, FunctionType &)
npe_arg(accum_func, AccumFunctionType &)
npe_arg(min_parallel, size_t)


npe_begin_code()

  igl::parallel_for(loop_size, prep_func, func, accum_func, min_parallel);
  return ;

npe_end_code()


