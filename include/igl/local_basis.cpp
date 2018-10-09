// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "local_basis.h"

#include <sstream>
#include <string>
#include <fstream>

#include <vector>
#include <Eigen/Geometry>


template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::local_basis(
  const Eigen::MatrixBase<DerivedV>& V,
  const Eigen::MatrixBase<DerivedF>& F,
  Eigen::PlainObjectBase<DerivedV>& B1,
  Eigen::PlainObjectBase<DerivedV>& B2,
  Eigen::PlainObjectBase<DerivedV>& B3
  )
{
  using namespace Eigen;
  using namespace std;
  B1.resize(F.rows(),3);
  B2.resize(F.rows(),3);
  B3.resize(F.rows(),3);

  for (unsigned i=0;i<F.rows();++i)
  {
      Eigen::Matrix<typename DerivedV::Scalar, 1, 3> v1 = (V.row(F(i,1)) - V.row(F(i,0))).normalized();
      Eigen::Matrix<typename DerivedV::Scalar, 1, 3> t = V.row(F(i,2)) - V.row(F(i,0));
      Eigen::Matrix<typename DerivedV::Scalar, 1, 3> v3 = v1.cross(t).normalized();
      Eigen::Matrix<typename DerivedV::Scalar, 1, 3> v2 = v1.cross(v3).normalized();

      B1.row(i) = v1;
      B2.row(i) = -v2;
      B3.row(i) = v3;
  }
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
// generated by autoexplicit.sh
template void igl::local_basis<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> >&);
template void igl::local_basis<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
#endif
