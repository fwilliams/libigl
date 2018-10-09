// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "piecewise_constant_winding_number.h"
#include "unique_edge_map.h"
#include "PI.h"

template <
  typename DerivedF,
  typename DeriveduE,
  typename uE2EType>
IGL_INLINE bool igl::piecewise_constant_winding_number(
  const Eigen::MatrixBase<DerivedF>& F,
  const Eigen::MatrixBase<DeriveduE>& uE,
  const std::vector<std::vector<uE2EType> >& uE2E)
{
  const size_t num_faces = F.rows();
  const size_t num_edges = uE.rows();
  const auto edge_index_to_face_index = [&](size_t ei)
  {
    return ei % num_faces;
  };
  const auto is_consistent = [&](size_t fid, size_t s, size_t d)
  {
    if ((size_t)F(fid, 0) == s && (size_t)F(fid, 1) == d) return true;
    if ((size_t)F(fid, 1) == s && (size_t)F(fid, 2) == d) return true;
    if ((size_t)F(fid, 2) == s && (size_t)F(fid, 0) == d) return true;

    if ((size_t)F(fid, 0) == d && (size_t)F(fid, 1) == s) return false;
    if ((size_t)F(fid, 1) == d && (size_t)F(fid, 2) == s) return false;
    if ((size_t)F(fid, 2) == d && (size_t)F(fid, 0) == s) return false;
    throw "Invalid face!!";
  };
  for (size_t i=0; i<num_edges; i++)
  {
    const size_t s = uE(i,0);
    const size_t d = uE(i,1);
    int count=0;
    for (const auto& ei : uE2E[i])
    {
      const size_t fid = edge_index_to_face_index(ei);
      if (is_consistent(fid, s, d))
      {
        count++;
      }
      else
      {
        count--;
      }
    }
    if (count != 0)
    {
      return false;
    }
  }
  return true;
}
template <typename DerivedF>
IGL_INLINE bool igl::piecewise_constant_winding_number(
  const Eigen::MatrixBase<DerivedF>& F)
{
  Eigen::Matrix<typename DerivedF::Scalar,Eigen::Dynamic,2> E, uE;
  Eigen::Matrix<typename DerivedF::Scalar,Eigen::Dynamic,1> EMAP;
  std::vector<std::vector<size_t> > uE2E;
  unique_edge_map(F, E, uE, EMAP, uE2E);
  return piecewise_constant_winding_number(F,uE,uE2E);
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
// generated by autoexplicit.sh
template bool igl::piecewise_constant_winding_number<Eigen::Matrix<int, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, unsigned long>(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 1, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, std::vector<std::vector<unsigned long, std::allocator<unsigned long> >, std::allocator<std::vector<unsigned long, std::allocator<unsigned long> > > > const&);
template bool igl::piecewise_constant_winding_number<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, unsigned long>(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, std::vector<std::vector<unsigned long, std::allocator<unsigned long> >, std::allocator<std::vector<unsigned long, std::allocator<unsigned long> > > > const&);
#ifdef WIN32
template bool igl::piecewise_constant_winding_number<class Eigen::Matrix<int, -1, -1, 0, -1, -1>, class Eigen::Matrix<int, -1, -1, 0, -1, -1>, unsigned __int64>(class Eigen::PlainObjectBase<class Eigen::Matrix<int, -1, -1, 0, -1, -1>> const &, class Eigen::PlainObjectBase<class Eigen::Matrix<int, -1, -1, 0, -1, -1>> const &, class std::vector<class std::vector<unsigned __int64, class std::allocator<unsigned __int64>>, class std::allocator<class std::vector<unsigned __int64, class std::allocator<unsigned __int64>>>> const &);
template bool igl::piecewise_constant_winding_number<class Eigen::Matrix<int, -1, 3, 1, -1, 3>, class Eigen::Matrix<int, -1, -1, 0, -1, -1>, unsigned __int64>(class Eigen::PlainObjectBase<class Eigen::Matrix<int, -1, 3, 1, -1, 3>> const &, class Eigen::PlainObjectBase<class Eigen::Matrix<int, -1, -1, 0, -1, -1>> const &, class std::vector<class std::vector<unsigned __int64, class std::allocator<unsigned __int64>>, class std::allocator<class std::vector<unsigned __int64, class std::allocator<unsigned __int64>>>> const &);
#endif
#endif
