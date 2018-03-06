// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2018 Francis Williams <francis@fwilliams.info>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IGL_MARCHING_TETS_H
#define IGL_MARCHING_TETS_H

#include "igl_inline.h"
#include <Eigen/Core>

namespace igl {
  // marching_tets( TV, TT, isovals, isovalue, outV, outF _
  //
  // performs the marching tetrahedra algorithm on a tet mesh defined by TV and TT
  // with scalar values defined at each vertex in TV. The output is a triangle
  // mesh approximating the isosurface
  // coresponding to the value isovalue.
  //
  // Input:
  //  TV  #tet_vertices x 3 array -- The vertices of the tetrahedral mesh
  //  TT  #tets x 4 array --  The indexes of each tet in the tetrahedral mesh
  //  isovals  #tet_vertices x 1 array -- The values defined on each tet vertex
  //  isovalue  scalar -- The isovalue of the level set we want to compute
  //
  // Output:
  //  outV  #V x 3 array -- The vertices of the output level surface mesh
  //  outF  #F x 3 array -- The face indexes of the output level surface mesh
  //
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
