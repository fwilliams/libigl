#include "marching_tets.h"

#include <unordered_map>
#include <vector>
#include <utility>
#include <cstdint>
#include <iostream>


inline std::int64_t make_edge_key(const std::pair<std::int32_t, std::int32_t>& p) {
    std::int64_t ret = 0;
    ret |= p.first;
    ret |= static_cast<std::int64_t>(p.second) << 32;
    return ret;
}

template <typename DerivedTV,
          typename DerivedTT,
          typename Derivedisovalues,
          typename DerivedoutV,
          typename DerivedoutF>
void igl::marching_tets(
        const Eigen::PlainObjectBase<DerivedTV>& TV,
        const Eigen::PlainObjectBase<DerivedTT>& TT,
        const Eigen::PlainObjectBase<Derivedisovalues>& isovals,
        double isovalue,
        Eigen::PlainObjectBase<DerivedoutV>& outV,
        Eigen::PlainObjectBase<DerivedoutF>& outF) {
    using namespace std;

    const int mt_cell_lookup[16][4] = {
        { -1, -1, -1, -1 }, // 0000 0a
        {  0,  2,  1, -1 }, // 0001 1a
        {  0,  3,  4, -1 }, // 0010 1b
        {  2,  1,  3,  4 }, // 0011 2b
        {  5,  3,  1, -1 }, // 0100 1d
        {  0,  2,  5,  3 }, // 0101 2a
        {  0,  1,  5,  4 }, // 0110 2d
        {  2,  5,  4, -1 }, // 0111 3a
        {  4,  5,  2, -1 }, // 1000 1c
        {  0,  4,  5,  1 }, // 1001 2c
        {  0,  3,  5,  2 }, // 1010 2f
        {  1,  3,  5, -1 }, // 1011 3d
        {  4,  3,  1,  2 }, // 1100 2e
        {  0,  4,  3, -1 }, // 1101 3c
        {  0,  1,  2, -1 }, // 1110 3b
        { -1, -1, -1, -1 }, // 1111 0b
    };

    const int mt_edge_lookup[6][2] = {
        {0, 1},
        {0, 2},
        {0, 3},
        {1, 2},
        {1, 3},
        {2, 3},
    };

    vector<Eigen::RowVector3d> vertices;
    vector<Eigen::RowVector3i> faces;
    vector<pair<int, int>> edge_table;


    assert(TT.cols() == 4 && TT.rows() >= 1);
    assert(TV.cols() == 3 && TV.rows() >= 4);
    assert(isovals.cols() == 1);

    static_assert(sizeof(std::int64_t) == 2*sizeof(std::int32_t), "need to fit 2 ints into a size_t");

    // For each tet
    for (int i = 0; i < TT.rows(); i++) {
        uint8_t key = 0;
        for (int v = 0; v < 4; v++) {
            int vid = TT(i, v);
            uint8_t flag = isovals[vid] > isovalue;
            key |= flag << v;
        }

        // Insert the vertices if they don't exist
        int v_ids[4] = {-1, -1, -1, -1}; // This will contain the index in TV of each vertex in the tet
        for (int e = 0; mt_cell_lookup[key][e] != -1 && e < 4; e++) {
            const int v1id = mt_edge_lookup[mt_cell_lookup[key][e]][0];
            const int v2id = mt_edge_lookup[mt_cell_lookup[key][e]][1];
            const int v1i = TT(i, v1id), v2i = TT(i, v2id);
            const Eigen::RowVector3d v1 = TV.row(v1i);
            const Eigen::RowVector3d v2 = TV.row(v2i);
            const double a = fabs(isovals[v1i] - isovalue), b = fabs(isovals[v2i] - isovalue);
            const double w = a / (a+b);

            const int vertex_id = vertices.size();
            vertices.push_back((1-w)*v1 + w*v2); // Push back the midpoint
            if (v1i < v2i) {
                edge_table.push_back(make_pair(v1i, v2i));
            } else {
                edge_table.push_back(make_pair(v2i, v1i));
            }

            v_ids[e] = vertex_id;
        }

        if (v_ids[0] != -1) {
            bool is_quad = mt_cell_lookup[key][3] != -1;
            if (is_quad) {
                Eigen::RowVector3i f1(v_ids[0], v_ids[1], v_ids[3]);
                Eigen::RowVector3i f2(v_ids[1], v_ids[2], v_ids[3]);
                faces.push_back(f1);
                faces.push_back(f2);
            } else {
                Eigen::RowVector3i f(v_ids[0], v_ids[1], v_ids[2]);
                faces.push_back(f);
            }
        }
    }

    // Deduplicate vertices
    int num_unique = 0;
    outV.resize(vertices.size(), 3);
    outF.resize(faces.size(), 3);
    unordered_map<int64_t, int> emap;
    for (int i = 0; i < faces.size(); i++) {
        for (int v = 0; v < 3; v++) {
            const int vi = faces[i][v];
            const pair<int32_t, int32_t> edge = edge_table[vi];
            const int64_t key = make_edge_key(edge);
            auto it = emap.find(key);
            if (it == emap.end()) { // New unique vertex, insert it
                outV.row(num_unique) = vertices[vi];
                outF(i, v) = num_unique;
                emap.insert(make_pair(key, num_unique));
                num_unique += 1;
            } else {
                outF(i, v) = it->second;
            }
        }
    }
    outV.conservativeResize(num_unique, 3);
}

#ifdef IGL_STATIC_LIBRARY
template void igl::marching_tets<typename Eigen::Matrix<double, -1, -1, 0, -1, -1>,
                                 typename Eigen::Matrix<int, -1, -1, 0, -1, -1>,
                                 typename Eigen::Matrix<double, -1, 1, 0, -1, 1>,
                                 typename Eigen::Matrix<double, -1, -1, 0, -1, -1>,
                                 typename Eigen::Matrix<int, -1, -1, 0, -1, -1>>(
        Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>> const&,
        Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1>> const&,
        const Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1>>&,
        double isovalue,
        Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>&,
        Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1>>&);

#endif // IGL_STATIC_LIBRARY
