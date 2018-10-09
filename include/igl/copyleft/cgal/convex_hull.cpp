// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2017 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "convex_hull.h"
#include "../../ismember.h"
#include "polyhedron_to_mesh.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_3.h>
#include <vector>

template <
  typename DerivedV,
  typename DerivedW,
  typename DerivedG>
IGL_INLINE void igl::copyleft::cgal::convex_hull(
  const Eigen::MatrixBase<DerivedV> & V,
  Eigen::PlainObjectBase<DerivedW> & W,
  Eigen::PlainObjectBase<DerivedG> & G)
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel      K;
  typedef K::Point_3                                              Point_3;
  //typedef CGAL::Delaunay_triangulation_3<K>                       Delaunay;
  //typedef Delaunay::Vertex_handle                                 Vertex_handle;
  //typedef CGAL::Surface_mesh<Point_3>                             Surface_mesh;
  typedef CGAL::Polyhedron_3<K>                             Polyhedron_3;
  std::vector<Point_3> points(V.rows());
  for(int i = 0;i<V.rows();i++)
  {
    points[i] = Point_3(V(i,0),V(i,1),V(i,2));
  }
  Polyhedron_3 poly;
  CGAL::convex_hull_3(points.begin(),points.end(),poly);
  assert(poly.is_pure_triangle() && "Assuming CGAL outputs a triangle mesh");
  polyhedron_to_mesh(poly,W,G);
}

template <
  typename DerivedV,
  typename DerivedF>
IGL_INLINE void igl::copyleft::cgal::convex_hull(
  const Eigen::MatrixBase<DerivedV> & V,
  Eigen::PlainObjectBase<DerivedF> & F)
{
  Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, Eigen::Dynamic> W;
  Eigen::Matrix<typename DerivedF::Scalar, Eigen::Dynamic, Eigen::Dynamic> G;
  convex_hull(V,W,G);
  // This is a lazy way to reindex into the original mesh
  Eigen::Matrix<bool,Eigen::Dynamic,1> I;
  Eigen::VectorXi J;
  igl::ismember_rows(W,V,I,J);
  assert(I.all() && "Should find all W in V");
  F.resizeLike(G);
  for(int f = 0;f<G.rows();f++)
  {
    for(int c = 0;c<3;c++)
    {
      F(f,c) = J(G(f,c));
    }
  }
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
// generated by autoexplicit.sh
template void igl::copyleft::cgal::convex_hull<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
#endif
