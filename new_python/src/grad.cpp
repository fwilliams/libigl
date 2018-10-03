#include <tuple>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <npe.h>
#include <typedefs.h>

#include <igl/grad.h>

const char* ds_grad = R"igl_Qu8mg5v7(
Compute the numerical gradient operator.

Parameters
----------
v : #v by 3 list of mesh vertex positions
f : #f by 3 list of mesh face indices [or a #faces by 4 list of tetrahedral indices]
uniform : boolean (default false). Use a uniform mesh instead of the vertices v

Returns
-------
g : #faces * dim by #v gradient operator
dtype : data-type of the returned objects, optional. Default is `float64`.
(All integer return types are `int32` by default.)

See also
--------
None

Notes
-----
Gradient of a scalar function defined on piecewise linear elements (mesh)
is constant on each triangle [tetrahedron] i,j,k:
grad(Xijk) = (Xj-Xi) * (Vi - Vk)^R90 / 2A + (Xk-Xi) * (Vj - Vi)^R90 / 2A
where Xi is the scalar value at vertex i, Vi is the 3D position of vertex
i, and A is the area of triangle (i,j,k). ^R90 represent a rotation of 
90 degrees.

Examples
--------
# Mesh in (v, f)
>>> g = grad(v, f)
)igl_Qu8mg5v7";

npe_function(grad)
npe_doc(ds_grad)

npe_arg(v, dense_f64, dense_f32)
npe_arg(f, dense_i32)
npe_default_arg(dtype, npe::dtype, "float64")
npe_default_arg(uniform, bool, false)

npe_begin_code()

if (dtype.type() == npe::type_f32) {
    Eigen::SparseMatrix<std::float_t> g;
    igl::grad(v, f, g, uniform);
    return npe::move(g);
} else {
    Eigen::SparseMatrix<std::double_t> g;
    igl::grad(v, f, g, uniform);
    return npe::move(g);
}


npe_end_code()


