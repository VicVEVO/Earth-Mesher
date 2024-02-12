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

// Function to create a sphere with vertices and faces based on altitude, latitude, and longitude
void displayPoints(const Eigen::VectorXf& altitudes, const Eigen::VectorXf& latitudes, const Eigen::VectorXf& longitudes) {
    int n = altitudes.size();

    // Create matrices to store the coordinates and colors of the points
    Eigen::MatrixXd V(n, 3);
    Eigen::MatrixXd C(n, 3);

    // Fill the matrices with random data
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

    // Set up the viewer
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_points(V, C);
    viewer.core().background_color = Eigen::Vector4f(0, 0, 0, 1);
    viewer.data().point_size = 5;

    // Launch the viewer
    viewer.launch(false,"Test");
}

int main(const int argc, const char *argv[]) {
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

        displayPoints(altitudes, latitudes, longitudes);
    }
    catch(const ErrorCode& errorCode) {
        std::cerr << errorMessages.at(errorCode) << std::endl;
    }
}
