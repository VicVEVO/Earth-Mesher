/**
 *  *description du projet*
 *
 * A optimiser
 *
 */
#include <igl/opengl/glfw/Viewer.h>

Eigen::MatrixXf getMatrixFromCSV(const std::string& csvFile) {
    std::ifstream csvData;
    csvData.open(csvFile);

    std::string lineRead;
    std::getline(csvData, lineRead);
    int nbFields = int(std::count(lineRead.begin(), lineRead.end(), ',')) + 1;
    int nbMeasurements = 0;
    Eigen::MatrixXf matrix(nbMeasurements, nbFields);

    int indexDataField; // number of the current column studied in lineRead
    bool nullData = false; // Is the data equal to zero at a certain point in time?
    float precision = 1e-6; // precision of floats

    // Read each line
    while (std::getline(csvData, lineRead)) {
        // Extend the matrix
        nbMeasurements++;
        matrix.conservativeResize(nbMeasurements,nbFields);
        indexDataField = 0;

        // Read the values in the line
        size_t pos = 0;
        std::string token;
        while ((pos = lineRead.find(',')) != std::string::npos) {
            token = lineRead.substr(0, pos);
            matrix(nbMeasurements-1,indexDataField++) = std::stof(token);
            lineRead.erase(0, pos + 1);
        }
        token = lineRead.substr(0, pos);
        matrix(nbMeasurements-1,indexDataField) = std::stof(token);

        // Verification if null values
        if (abs(matrix(nbMeasurements-1,indexDataField)) < precision) {
            nullData = true;
            for (int i=0 ; i < nbFields; i++) {
                if (abs(matrix(nbMeasurements-1,i)) > precision) {
                    nullData = false;
                }
            }
            if (nullData) {
                matrix.conservativeResize(nbMeasurements-1,nbFields);
                return matrix;
            }
        }
    };
    return matrix;
}

float radiansToDegrees(float radians) {
    return radians * 180.0 / M_PI;
}

// Function to create a sphere with vertices and faces based on altitude, latitude, and longitude
void createSphere(const Eigen::VectorXf& altitudes, const Eigen::VectorXf& latitudes, const Eigen::VectorXf& longitudes, Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
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
            float altitude = altitudes[i];
            float latitude = latitudes[i];
            float longitude = longitudes[j];

            // Calculation of the spherical coordinates of latitude and longitude
            float theta = longitude;
            float phi = M_PI / 2 - latitude;

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
    std::string csvFile = "/home/vic_pabo/Documents/Earth-Mesher/submodules/NC-Converter/data/csv_files/TP_GPN_2PfP314_006_20010323_223727_20010323_233339.csv";
    std::cout << csvFile <<std::endl;
    Eigen::MatrixXf dataMat = getMatrixFromCSV(csvFile);

    Eigen::VectorXf altitudes = dataMat.col(0);
    Eigen::VectorXf latitudes = dataMat.col(1);
    Eigen::VectorXf longitudes = dataMat.col(2);

    createSphere(altitudes, latitudes, longitudes, V, F);
    // Mesh plotting
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);
    viewer.data().set_face_based(true);
    viewer.launch();

}
