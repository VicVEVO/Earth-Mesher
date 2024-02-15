/**
 *  Project description (bruh)
 *
 */
#include <igl/opengl/glfw/Viewer.h>
#include "errorCode.h"
#include "errorMessages.h"


// Function to convert a .csv File into a Eigen Matrix
Eigen::MatrixXf getMatrixFromCSV(const std::string& csvFile) {
    std::ifstream csvData;
    csvData.open(csvFile);

    std::string lineRead;
    std::getline(csvData, lineRead);

    // Checking if the csvFile actually exists
    if (lineRead.empty()){
        throw ErrorCode::FileNotFound;
    }

    int nbFields = int(std::count(lineRead.begin(), lineRead.end(), ',')) + 1;
    if (nbFields < 3){
        throw ErrorCode::NotEnoughFields;
    }
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

// Function to create a sphere with vertices and faces based on altitude, latitude, and longitude
void setPoints(const int n, const Eigen::VectorXf& altitudes, const Eigen::VectorXf& latitudes, const Eigen::VectorXf& longitudes, Eigen::MatrixXd& V, Eigen::MatrixXd& C) {

    // Fill the matrices with data
    for (int i = 0; i < n; ++i) {
        double theta = longitudes(i) + M_PI/2;
        double phi = latitudes(i);

        double x = sin(phi) * cos(theta);
        double y = sin(phi) * sin(theta);
        double z = cos(phi);
        V.row(i) = (altitudes(i)/1356038) * Eigen::Vector3d(x, y, z);

        // Set colors for visualization
        C(i, 0) = 255;
        C(i, 1) = 255;
        C(i, 2) = 255;
    }
}



/*
void color(float altitude, float x, float y, float z, int R1, int G1, int B1, int R2, int G2, int B2, int R3, int G3, int B3, int& R, int& G, int& B) {
    if (altitude <= x) {
        R = R1;
        G = G1;
        B = B1;
    } else if (altitude <= y) {
        float t = (altitude - x) / (y - x);
        R = R1 + t * (R2 - R1);
        G = G1 + t * (G2 - G1);
        B = B1 + t * (B2 - B1);
    } else if (altitude <= z) {
        float t = (altitude - y) / (z - y);
        R = R2 + t * (R3 - R2);
        G = G2 + t * (G3 - G2);
        B = B2 + t * (B3 - B2);
    } else {
        R = R3;
        G = G3;
        B = B3;
    }
}
*/


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

    // Populate the face matrix
    for (int i = 0; i < N; ++i) {
        for (int     j = 0; j < N; ++j) {
            int index = i * N + j;
            F.row(2 * index) << i * (N + 1) + j, (i + 1) * (N + 1) + j, i * (N + 1) + j + 1;
            F.row(2 * index + 1) << (i + 1) * (N + 1) + j, (i + 1) * (N + 1) + j + 1, i * (N + 1) + j + 1;
        }
    }
}



void createLines(const Eigen::VectorXf& altitudes, const Eigen::VectorXf& latitudes, const Eigen::VectorXf& longitudes, Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
    // Get the number of vertices
    int N = altitudes.rows();

    // Resize the vertex and face matrices
    V.resize(N, 3);
    F.resize(N - 2, 3);

    // Populate the vertex matrix
    for (int i = 0; i < N; i++) {
        double theta = longitudes(i) + M_PI/2;
        double phi = latitudes(i);

        double x = sin(phi) * cos(theta);
        double y = sin(phi) * sin(theta);
        double z = cos(phi);


        // Store the spherical coordinates in the vertex matrix
        V(i, 0) = x;
        V(i, 1) = y;
        V(i, 2) = z;
    }

    // Populate the face matrix
    for (int i = 0; i < N - 2; ++i) {
        // Find the three vertices for this triangle
        int idx1 = i;
        int idx2 = (i + 1) % (N - 1);
        int idx3 = (i + 2) % (N - 1) + 1;

        // Compute the latitudes and longitudes of the three vertices
        double lat1 = latitudes(idx1);
        double lat2 = latitudes(idx2);
        double lat3 = latitudes(idx3);
        double lon1 = longitudes(idx1);
        double lon2 = longitudes(idx2);
        double lon3 = longitudes(idx3);

        // Find the nearest vertex for each edge
        int j1, j2, j3;
        double min_dist, dist;

        // Edge 1-2
        j1 = idx2;
        j2 = idx1;
        min_dist = std::numeric_limits<double>::max();
        for (int j = 0; j < N; ++j) {
            if (j == idx1 || j == idx2) continue;
            dist = std::abs(latitudes(j) - lat1) + std::abs(longitudes(j) - lon1);
            if (dist < min_dist) {
                min_dist = dist;
                j1 = j;
            }
        }

        // Edge 2-3
        j2 = idx3;
        min_dist = std::numeric_limits<double>::max();
        for (int j = 0; j < N; ++j) {
            if (j == idx2 || j == idx3) continue;
            dist = std::abs(latitudes(j) - lat2) + std::abs(longitudes(j) - lon2);
            if (dist < min_dist) {
                min_dist = dist;
                j2 = j;
            }
        }

        // Edge 3-1
        j3 = idx1;
        min_dist = std::numeric_limits<double>::max();
        for (int j = 0; j < N; ++j) {
            if (j == idx3 || j == idx1) continue;
            dist = std::abs(latitudes(j) - lat3) + std::abs(longitudes(j) - lon3);
            if (dist < min_dist) {
                min_dist = dist;
                j3 = j;
            }
        }

        // Populate the face matrix
        F.row(i) << idx1, j1, j3;
    }
}

// Display procedures

void viewPoints(const Eigen::MatrixXd& V, const Eigen::MatrixXd& C) {
    // Set up the viewer
    igl::opengl::glfw::Viewer viewer;

    // Set the viewer data
    viewer.data().set_points(V, C);
    viewer.data().point_size = 5;

    // Set the background color
    viewer.core().background_color = Eigen::Vector4f(0, 0, 0, 1);

    // Launch the viewer
    viewer.launch(false,"Example window with lines");
}


void viewMesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) {
    // Set up the viewer
    igl::opengl::glfw::Viewer viewer;

    // Set the viewer data
    viewer.data().set_mesh(V, F);

    // Set the background color
    viewer.core().background_color = Eigen::Vector4f(0, 0, 0, 1);

    // Set the colors
    Eigen::MatrixXd C(V.rows(), 3);
    C.setConstant(255);
    viewer.data().set_colors(C);

    // Disable the lines
    viewer.data().show_lines = false;

    // Launch the viewer
    viewer.launch(false,"Example window with lines");
}


int main(const int argc, const char *argv[]) {
    // To unmute
    try {
        if (argc == 1 | argc == 2) {
            throw ErrorCode::NotEnoughArguments;
        }

        // Get the plot type and CSV file name from the command-line arguments
        std::string plotType = argv[1];
        std::string csvFile = argv[2];

        // Read the matrix from the CSV file
        Eigen::MatrixXf dataMat = getMatrixFromCSV(csvFile);

        // Extract the altitudes, latitudes, and longitudes from the matrix
        Eigen::VectorXf altitudes = dataMat.col(0);
        Eigen::VectorXf latitudes = dataMat.col(1);
        Eigen::VectorXf longitudes = dataMat.col(2);

        // Extract the size of the data set
        int n = altitudes.size();

        if (plotType == "-P") {
            // Create matrices to store the coordinates and colors of the points
            Eigen::MatrixXd V(n, 3);
            Eigen::MatrixXd C(n, 3);

            setPoints(n, altitudes, latitudes, longitudes, V, C);

            // Display data
            viewPoints(V, C);
        } else {
            Eigen::MatrixXd V;
            Eigen::MatrixXi F;

            //createLines(altitudes, latitudes, longitudes, V, F);

            createSphere(100,V,F);

            viewMesh(V,F);

            // Display data
        }
    }

    catch(const ErrorCode& errorCode) {
        std::cerr << errorMessages.at(errorCode) << std::endl;
    }
}

