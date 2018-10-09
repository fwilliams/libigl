#include <npe.h>
#include <typedefs.h>






#include <igl/mexErrMsgTxt.h>

const char* ds_mex_err_msg_txt = R"igl_Qu8mg5v7(

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

 Wrapper for mexErrMsgTxt that only calls error if test fails
)igl_Qu8mg5v7";

npe_function(mex_err_msg_txt)
npe_doc(ds_mex_err_msg_txt)

npe_arg(message, char *)


npe_begin_code()

  bool test;
  igl::  matlab::mexErrMsgTxt(message, test);
  return npe::move(test);

npe_end_code()


