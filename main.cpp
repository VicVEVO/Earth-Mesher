#include "build/_deps/libigl-src/include/igl/opengl/glfw/Viewer.cpp"

void createSphere(int N, Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
    V.resize((N + 1) * (N + 1), 3);
    F.resize(N * N * 2, 3);

    for (int i = 0; i <= N; ++i) {
        for (int j = 0; j <= N; ++j) {
            double theta = 2.0 * M_PI * static_cast<double>(i) / static_cast<double>(N);
            double phi = M_PI * static_cast<double>(j) / static_cast<double>(N);

            V.row(i * (N + 1) + j) << std::sin(phi) * std::cos(theta), std::sin(phi) * std::sin(theta), std::cos(phi);
        }
    }

    for (int i = 0; i < N; ++i) {
        for (int     j = 0; j < N; ++j) {
            int index = i * N + j;
            F.row(2 * index) << i * (N + 1) + j, (i + 1) * (N + 1) + j, i * (N + 1) + j + 1;
            F.row(2 * index + 1) << (i + 1) * (N + 1) + j, (i + 1) * (N + 1) + j + 1, i * (N + 1) + j + 1;
        }
    }
}

int main(int argc, char* argv[]) {
    // Sphere creation
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    createSphere(10, V, F);

    // Mesh plotting
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);
    viewer.data().set_face_based(true);
    viewer.launch();
}
