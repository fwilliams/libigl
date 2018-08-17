#include <tuple>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <npe.h>
#include <igl/readOFF.h>

npe_function(readOFF)
npe_arg(str, std::string)
npe_default_arg(read_normals, bool, true)

npe_begin_code()
using namespace std;

Eigen::Matrix<std::double_t, Eigen::Dynamic, Eigen::Dynamic> v;
Eigen::Matrix<std::int32_t, Eigen::Dynamic, Eigen::Dynamic> f;
Eigen::Matrix<std::double_t, Eigen::Dynamic, Eigen::Dynamic> n;
bool ret;
//bool read_normals = false;

if (read_normals) {
    ret = igl::readOFF(str, v, f, n);
}
else {
    ret = igl::readOFF(str, v, f);
}
if (!ret) {
    throw std::invalid_argument("File '" + str + "' not found.");
}
return std::make_tuple(NPE_MOVE_DENSE(v), NPE_MOVE_DENSE(f), NPE_MOVE_DENSE(n));

npe_end_code()
