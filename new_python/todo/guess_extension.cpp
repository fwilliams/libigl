#include <npe.h>
#include <typedefs.h>






#include <igl/guess_extension.h>

const char* ds_guess_extension = R"igl_Qu8mg5v7(

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

 Given a file pointer at the beginning of a "mesh" file, try to guess the
   extension of the file format it comes from. The file pointer is rewound on
   return.
  
   Inputs:
     fp  file pointer (see output)
   Outputs:
     fp  file pointer rewound 
     guess  extension as string. One of "mesh",{"obj"},"off","ply","stl", or
       "wrl"
  
)igl_Qu8mg5v7";

npe_function(guess_extension)
npe_doc(ds_guess_extension)



npe_begin_code()

  FILE * fp;
  std::string & guess;
  igl::guess_extension(fp, guess);
  return std::make_tuple(npe::move(fp), npe::move(guess));

npe_end_code()
#include <igl/guess_extension.h>

const char* ds_guess_extension = R"igl_Qu8mg5v7(
See guess_extension for the documentation.
)igl_Qu8mg5v7";

npe_function(guess_extension)
npe_doc(ds_guess_extension)



npe_begin_code()

  FILE * fp;
  igl::guess_extension(fp);
  return npe::move(fp);

npe_end_code()


