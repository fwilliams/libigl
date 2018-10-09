#include <npe.h>
#include <typedefs.h>






#include <igl/validate_arg_scalar.h>

const char* ds_validate_arg_scalar = R"igl_Qu8mg5v7(

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

 Throw an error if arg i+1 is not a scalar
    
     Inputs:
       i  index of current argument
       nrhs  total number of arguments
       prhs  pointer to arguments array
       name   name of current argument
)igl_Qu8mg5v7";

npe_function(validate_arg_scalar)
npe_doc(ds_validate_arg_scalar)

npe_arg(i, int)
npe_arg(nrhs, int)
npe_arg(prhs, int *[])
npe_arg(name, char *)


npe_begin_code()

  igl::  matlab::validate_arg_scalar(i, nrhs, prhs, name);
  return ;

npe_end_code()
#include <igl/validate_arg_logical.h>

const char* ds_validate_arg_logical = R"igl_Qu8mg5v7(
See validate_arg_logical for the documentation.
)igl_Qu8mg5v7";

npe_function(validate_arg_logical)
npe_doc(ds_validate_arg_logical)

npe_arg(i, int)
npe_arg(nrhs, int)
npe_arg(prhs, int *[])
npe_arg(name, char *)


npe_begin_code()

  igl::  matlab::validate_arg_logical(i, nrhs, prhs, name);
  return ;

npe_end_code()
#include <igl/validate_arg_char.h>

const char* ds_validate_arg_char = R"igl_Qu8mg5v7(
See validate_arg_char for the documentation.
)igl_Qu8mg5v7";

npe_function(validate_arg_char)
npe_doc(ds_validate_arg_char)

npe_arg(i, int)
npe_arg(nrhs, int)
npe_arg(prhs, int *[])
npe_arg(name, char *)


npe_begin_code()

  igl::  matlab::validate_arg_char(i, nrhs, prhs, name);
  return ;

npe_end_code()
#include <igl/validate_arg_double.h>

const char* ds_validate_arg_double = R"igl_Qu8mg5v7(
See validate_arg_double for the documentation.
)igl_Qu8mg5v7";

npe_function(validate_arg_double)
npe_doc(ds_validate_arg_double)

npe_arg(i, int)
npe_arg(nrhs, int)
npe_arg(prhs, int *[])
npe_arg(name, char *)


npe_begin_code()

  igl::  matlab::validate_arg_double(i, nrhs, prhs, name);
  return ;

npe_end_code()
#include <igl/validate_arg_function_handle.h>

const char* ds_validate_arg_function_handle = R"igl_Qu8mg5v7(
See validate_arg_function_handle for the documentation.
)igl_Qu8mg5v7";

npe_function(validate_arg_function_handle)
npe_doc(ds_validate_arg_function_handle)

npe_arg(i, int)
npe_arg(nrhs, int)
npe_arg(prhs, int *[])
npe_arg(name, char *)


npe_begin_code()

  igl::  matlab::validate_arg_function_handle(i, nrhs, prhs, name);
  return ;

npe_end_code()


