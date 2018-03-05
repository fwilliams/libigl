#include <igl/opengl/glfw/Viewer.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/readOFF.h>
#include <igl/colormap.h>
#include <igl/marching_tets.h>

#include "tutorial_shared_path.h"

// Input polygon
Eigen::MatrixXd V;
Eigen::MatrixXi F;

// Tetrahedralized interior
Eigen::MatrixXd TV;
Eigen::MatrixXi TT;
Eigen::MatrixXi TF;

// Isovalues stored at each tet vertex
Eigen::VectorXd isovalues;

// Vertices and faces of the level set computed by marching tets
Eigen::MatrixXd isoV;
Eigen::MatrixXi isoF;


// Visualize the tet mesh as a wireframe
void visualize_tet_wireframe(igl::opengl::glfw::Viewer& viewer, const Eigen::MatrixXd& TV, const Eigen::MatrixXi& TT, const Eigen::VectorXd& isovals) {
    std::vector<std::pair<int, int>> edges;
    for (int i = 0; i < TT.rows(); i++) {
        int tf1 = TT(i, 0);
        int tf2 = TT(i, 1);
        int tf3 = TT(i, 2);
        int tf4 = TT(i, 2);
        edges.push_back(std::make_pair(tf1, tf2));
        edges.push_back(std::make_pair(tf1, tf3));
        edges.push_back(std::make_pair(tf1, tf4));
        edges.push_back(std::make_pair(tf2, tf3));
        edges.push_back(std::make_pair(tf2, tf4));
        edges.push_back(std::make_pair(tf3, tf4));
    }

    Eigen::MatrixXd v1(edges.size(), 3), v2(edges.size(), 3);
    for (int i = 0; i < edges.size(); i++) {
        v1.row(i) = TV.row(edges[i].first);
        v2.row(i) = TV.row(edges[i].second);
    }

    Eigen::MatrixXd C;
    Eigen::VectorXd isovals_normalized = (isovals - isovals.minCoeff() * Eigen::VectorXd::Ones(isovals.size())) / (isovals.maxCoeff() - isovals.minCoeff());
    igl::colormap(igl::COLOR_MAP_TYPE_MAGMA, isovals_normalized, false, C);
    viewer.data().point_size = 5.0;
    viewer.data().add_points(TV, C);
    viewer.data().add_edges(v1, v2, Eigen::RowVector3d(0.1, 0.1, 0.1));
}


int main(int argc, char *argv[]) {
    using namespace Eigen;
    using namespace std;

    // Load a surface mesh
    igl::readOFF(TUTORIAL_SHARED_PATH "/fertility.off",V,F);

    // Tetrahedralize the interior
    igl::copyleft::tetgen::tetrahedralize(V,F,"pq1.414Y", TV,TT,TF);

    // Make the isovalues of each tet vertex, be the distance form the origin
    isovalues.resize(TV.rows());
    for (int i = 0; i < isovalues.size(); i++) {
        isovalues[i] = TV.row(i).norm();
    }

    // We'll compute the level set with this value
    double level_set_value = 45.0;

    igl::marching_tets(TV, TT, isovalues, level_set_value, isoV, isoF);

    // Plot the generated mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(isoV, isoF);
    visualize_tet_wireframe(viewer, TV, TT, isovalues);
    viewer.launch();
}
