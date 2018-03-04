#ifndef IGL_MARCHING_TETS_H
#define IGL_MARCHING_TETS_H

#include "igl_inline.h"
#include <Eigen/Core>

namespace igl {
template <typename DerivedTV,
          typename DerivedTT,
          typename Derivedisovalues,
          typename DerivedoutV,
          typename DerivedoutF>
IGL_INLINE void marching_tets(
        const Eigen::PlainObjectBase<DerivedTV>& TV,
        const Eigen::PlainObjectBase<DerivedTT>& TT,
        const Eigen::PlainObjectBase<Derivedisovalues>& isovals,
        double isovalue,
        Eigen::PlainObjectBase<DerivedoutV>& outV,
        Eigen::PlainObjectBase<DerivedoutF>& outF);
}

#ifndef IGL_STATIC_LIBRARY
#  include "marching_tets.cpp"
#endif

#endif // IGL_MARCHING_TETS_H
