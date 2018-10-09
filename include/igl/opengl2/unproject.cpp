// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "unproject.h"
#include "model_proj_viewport.h"
#include "../unproject.h"
#include "gl.h"

#include <Eigen/Dense>
#include <Eigen/LU>

IGL_INLINE void igl::opengl2::unproject(
  const double winX,
  const double winY,
  const double winZ,
  double* objX,
  double* objY,
  double* objZ)
{
  Eigen::Vector3d obj;
  igl::opengl2::unproject(Eigen::Vector3d(winX,winY,winZ),obj);
  *objX = obj(0);
  *objY = obj(1);
  *objZ = obj(2);
}

template <typename Derivedwin, typename Derivedobj>
IGL_INLINE void igl::opengl2::unproject(
  const Eigen::MatrixBase<Derivedwin> & win,
  Eigen::PlainObjectBase<Derivedobj> & obj)
{
  const auto ret = igl::opengl2::unproject(win);
  obj = ret.template cast<typename Derivedobj::Scalar>();
}


template <typename Derivedwin>
IGL_INLINE Derivedwin igl::opengl2::unproject(
  const Eigen::MatrixBase<Derivedwin> & win)
{
  using namespace Eigen;
  typedef typename Derivedwin::Scalar Scalar;
  Matrix4d MV,P;
  Vector4d VPd;
  model_proj_viewport(MV,P,VPd);
  Vector3d wind = win.template cast<double>();
  Vector3d objd = igl::unproject(wind,MV,P,VPd);
  return objd.template cast<Scalar>();
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
// generated by autoexplicit.sh
template void igl::opengl2::unproject<Eigen::Matrix<double, 3, 1, 0, 3, 1>, Eigen::Matrix<float, 3, 1, 0, 3, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, 3, 1, 0, 3, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<float, 3, 1, 0, 3, 1> >&);
template void igl::opengl2::unproject<Eigen::Matrix<float, 3, 1, 0, 3, 1>, Eigen::Matrix<float, 3, 1, 0, 3, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<float, 3, 1, 0, 3, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<float, 3, 1, 0, 3, 1> >&);
template Eigen::Matrix<float, 3, 1, 0, 3, 1> igl::opengl2::unproject<Eigen::Matrix<float, 3, 1, 0, 3, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<float, 3, 1, 0, 3, 1> > const&);
template Eigen::Matrix<double, 3, 1, 0, 3, 1> igl::opengl2::unproject<Eigen::Matrix<double, 3, 1, 0, 3, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, 3, 1, 0, 3, 1> > const&);
template void igl::opengl2::unproject<Eigen::Matrix<double, 3, 1, 0, 3, 1>, Eigen::Matrix<double, 3, 1, 0, 3, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, 3, 1, 0, 3, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 3, 1, 0, 3, 1> >&);
#endif
