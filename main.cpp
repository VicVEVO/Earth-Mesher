#include <igl/opengl/glfw/Viewer.h>

Eigen::MatrixXf matrixFromCSV(const std::string& csvFile) {
    std::ifstream csvData(csvFile);

    // Determining the number of rows and columns of the matrix
    std::string matrixSizeString;
    std::getline(csvData, matrixSizeString);
    std::istringstream matrixSizeStream(matrixSizeString);
    int rows, cols;
    matrixSizeStream >> rows >> cols;


    // Creation of the matrix
    Eigen::MatrixXf matrix(rows, cols);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::string matrixEntry;
            std::getline(csvData, matrixEntry, ',');
            matrix(i, j) = std::stof(matrixEntry);
        }
    }
    csvData.close();
    return matrix;
}

double radiansToDegrees(double radians) {
    return radians * 180.0 / M_PI;
}

// Function to create a sphere with vertices and faces based on altitude, latitude, and longitude
void createSphere(const Eigen::VectorXd& altitudes, const Eigen::VectorXd& latitudes, const Eigen::VectorXd& longitudes, Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
    // Determine the number of vertices and faces
    int N = altitudes.cols();
    int numVertices = (N + 1) * (N + 1);
    int numFaces = 2 * N * N;

    // Initialize the vertices and faces matrices
    V.resize(numVertices, 3);
    F.resize(numFaces, 3);

    // Fill in the vertices based on altitude, latitude, and longitude
    for (int i = 0; i <= N; ++i) {
        for (int j = 0; j <= N; ++j) {
            double altitude = altitudes[i];
            double latitude = latitudes[i];
            double longitude = longitudes[j];

            // Calculation of the spherical coordinates of latitude and longitude
            double theta = longitude;
            double phi = M_PI / 2 - latitude;

            // Add the altitude
            V.row(i * (N + 1) + j) << std::sin(phi) * std::cos(theta) * altitude, std::sin(phi) * std::sin(theta) * altitude, std::cos(phi) * altitude;
        }
    }

    // Fill in the faces
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
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
    std::cout << "bruh" << std::endl;
    std::string csvFile = "/home/vic_pabo/Documents/Earth-Mesher/submodules/NC-Converter/data/csv_files/TP_GPN_2PfP314_006_20010323_223727_20010323_233339.csv";
    std::cout << csvFile <<std::endl;
    Eigen::MatrixXf dataMat = matrixFromCSV(csvFile);
    std::cout << dataMat << std::endl;
    /*
    createSphere(10, V, F);
    // Mesh plotting
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);
    viewer.data().set_face_based(true);
    viewer.launch();
    */
}
