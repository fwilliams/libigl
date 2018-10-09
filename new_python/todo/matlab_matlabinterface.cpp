#include <npe.h>
#include <typedefs.h>






#include <igl/mlinit.h>

const char* ds_mlinit = R"igl_Qu8mg5v7(

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

 Init the MATLAB engine
     (no need to call it directly since it is automatically invoked by any other command)
)igl_Qu8mg5v7";

npe_function(mlinit)
npe_doc(ds_mlinit)



npe_begin_code()

  int ** engine;
  igl::  matlab::mlinit(engine);
  return npe::move(engine);

npe_end_code()
#include <igl/mlclose.h>

const char* ds_mlclose = R"igl_Qu8mg5v7(

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

 Closes the MATLAB engine
)igl_Qu8mg5v7";

npe_function(mlclose)
npe_doc(ds_mlclose)



npe_begin_code()

  int ** engine;
  igl::  matlab::mlclose(engine);
  return npe::move(engine);

npe_end_code()
#include <igl/mlsetmatrix.h>

const char* ds_mlsetmatrix = R"igl_Qu8mg5v7(

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

 Send a matrix to MATLAB
)igl_Qu8mg5v7";

npe_function(mlsetmatrix)
npe_doc(ds_mlsetmatrix)

npe_arg(m, Eigen::MatrixXd &)


npe_begin_code()

  int ** engine;
  std::string name;
  igl::  matlab::mlsetmatrix(m, engine, name);
  return std::make_tuple(npe::move(engine), npe::move(name));

npe_end_code()
#include <igl/mlsetmatrix.h>

const char* ds_mlsetmatrix = R"igl_Qu8mg5v7(

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

 Send a matrix to MATLAB
)igl_Qu8mg5v7";

npe_function(mlsetmatrix)
npe_doc(ds_mlsetmatrix)

npe_arg(m, Eigen::MatrixXf &)


npe_begin_code()

  int ** engine;
  std::string name;
  igl::  matlab::mlsetmatrix(m, engine, name);
  return std::make_tuple(npe::move(engine), npe::move(name));

npe_end_code()
#include <igl/mlsetmatrix.h>

const char* ds_mlsetmatrix = R"igl_Qu8mg5v7(

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

 Send a matrix to MATLAB
)igl_Qu8mg5v7";

npe_function(mlsetmatrix)
npe_doc(ds_mlsetmatrix)

npe_arg(m, Eigen::MatrixXi &)


npe_begin_code()

  int ** engine;
  std::string name;
  igl::  matlab::mlsetmatrix(m, engine, name);
  return std::make_tuple(npe::move(engine), npe::move(name));

npe_end_code()
#include <igl/mlsetmatrix.h>

const char* ds_mlsetmatrix = R"igl_Qu8mg5v7(

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

 Send a matrix to MATLAB
)igl_Qu8mg5v7";

npe_function(mlsetmatrix)
npe_doc(ds_mlsetmatrix)

npe_arg(m, int &)


npe_begin_code()

  int ** mlengine;
  std::string name;
  igl::  matlab::mlsetmatrix(m, mlengine, name);
  return std::make_tuple(npe::move(mlengine), npe::move(name));

npe_end_code()
#include <igl/mlgetmatrix.h>

const char* ds_mlgetmatrix = R"igl_Qu8mg5v7(

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

 Receive a matrix from MATLAB
)igl_Qu8mg5v7";

npe_function(mlgetmatrix)
npe_doc(ds_mlgetmatrix)



npe_begin_code()

  int ** engine;
  std::string name;
  Eigen::MatrixXd & m;
  igl::  matlab::mlgetmatrix(engine, name, m);
  return std::make_tuple(npe::move(engine), npe::move(name), npe::move(m));

npe_end_code()
#include <igl/mlgetmatrix.h>

const char* ds_mlgetmatrix = R"igl_Qu8mg5v7(

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

 Receive a matrix from MATLAB
)igl_Qu8mg5v7";

npe_function(mlgetmatrix)
npe_doc(ds_mlgetmatrix)



npe_begin_code()

  int ** engine;
  std::string name;
  Eigen::MatrixXf & m;
  igl::  matlab::mlgetmatrix(engine, name, m);
  return std::make_tuple(npe::move(engine), npe::move(name), npe::move(m));

npe_end_code()
#include <igl/mlgetmatrix.h>

const char* ds_mlgetmatrix = R"igl_Qu8mg5v7(

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

 Receive a matrix from MATLAB
)igl_Qu8mg5v7";

npe_function(mlgetmatrix)
npe_doc(ds_mlgetmatrix)



npe_begin_code()

  int ** engine;
  std::string name;
  Eigen::MatrixXi & m;
  igl::  matlab::mlgetmatrix(engine, name, m);
  return std::make_tuple(npe::move(engine), npe::move(name), npe::move(m));

npe_end_code()
#include <igl/mlgetmatrix.h>

const char* ds_mlgetmatrix = R"igl_Qu8mg5v7(

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

 Receive a matrix from MATLAB
)igl_Qu8mg5v7";

npe_function(mlgetmatrix)
npe_doc(ds_mlgetmatrix)



npe_begin_code()

  int ** mlengine;
  std::string name;
  int & m;
  igl::  matlab::mlgetmatrix(mlengine, name, m);
  return std::make_tuple(npe::move(mlengine), npe::move(name), npe::move(m));

npe_end_code()
#include <igl/mlsetscalar.h>

const char* ds_mlsetscalar = R"igl_Qu8mg5v7(

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

 Send a single scalar to MATLAB
)igl_Qu8mg5v7";

npe_function(mlsetscalar)
npe_doc(ds_mlsetscalar)



npe_begin_code()

  int ** engine;
  std::string name;
  double s;
  igl::  matlab::mlsetscalar(engine, name, s);
  return std::make_tuple(npe::move(engine), npe::move(name), npe::move(s));

npe_end_code()
#include <igl/mlgetscalar.h>

const char* ds_mlgetscalar = R"igl_Qu8mg5v7(

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

 Receive a single scalar from MATLAB
)igl_Qu8mg5v7";

npe_function(mlgetscalar)
npe_doc(ds_mlgetscalar)



npe_begin_code()

  int ** engine;
  std::string name;
  igl::  matlab::mlgetscalar(engine, name);
  return std::make_tuple(npe::move(engine), npe::move(name));

npe_end_code()
#include <igl/mleval.h>

const char* ds_mleval = R"igl_Qu8mg5v7(

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

 Execute arbitrary MATLAB code and return the MATLAB output
)igl_Qu8mg5v7";

npe_function(mleval)
npe_doc(ds_mleval)



npe_begin_code()

  int ** engine;
  std::string code;
  igl::  matlab::mleval(engine, code);
  return std::make_tuple(npe::move(engine), npe::move(code));

npe_end_code()
#include <igl/mlsetmatrix.h>

const char* ds_mlsetmatrix = R"igl_Qu8mg5v7(

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

 Send a sparse matrix to MATLAB
)igl_Qu8mg5v7";

npe_function(mlsetmatrix)
npe_doc(ds_mlsetmatrix)

npe_arg(m, int &)


npe_begin_code()

  int ** mlengine;
  std::string name;
  igl::  matlab::mlsetmatrix(m, mlengine, name);
  return std::make_tuple(npe::move(mlengine), npe::move(name));

npe_end_code()


