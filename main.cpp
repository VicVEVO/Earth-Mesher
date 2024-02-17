/**
 *  Project description (bruh)
 *
 */
#include <igl/opengl/glfw/Viewer.h>
#include "errorCode.h"
#include "errorMessages.h"
#include "constants.h"


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
    if (nbFields < NUMBER_FIELDS_MIN){
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
        C(i, ALTITUDE_FIELD_INDEX) = MAX_COLOR;
        C(i, LATITUDE_FIELD_INDEX) = MAX_COLOR;
        C(i, LONGITUDE_FIELD_INDEX) = MAX_COLOR;
    }
}

void color(const float altitude, const float x, const float y, const float z, const int R1, const int G1, const int B1, const int R2, const int G2, const int B2, const int R3, const int G3, const int B3, int& R, int& G, int& B) {
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

int indicePlusProche(float latitude, float longitude, Eigen::VectorXf& latitudes, Eigen::VectorXf& longitudes) {
    float distanceMin = (latitudes(0) - latitude) * (latitudes(0) - latitude) + (longitudes(0) - longitude)*(longitudes(0) - longitude);
    int indiceMin = 0;
    for (int i = 1; i < latitudes.size(); i++) {
        if ((latitudes(i) - latitude) * (latitudes(i) - latitude) + (longitudes(i) - longitude)*(longitudes(i) - longitude) < distanceMin ) {
            distanceMin = (latitudes(i) - latitude) * (latitudes(i) - latitude) + (longitudes(i) - longitude)*(longitudes(i) - longitude);
            indiceMin = i;
        }
    }
    return indiceMin;
}

void createColorSphere(int N, Eigen::VectorXf& altitudes, Eigen::VectorXf& latitudes, Eigen::VectorXf& longitudes, Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& C) {
    // Initialisation of the colors of the triangles
    int R,G,B;

    // Initialisation of the level of a triangle
    float level;

    //std::cout << latitudes << std::endl; entre 0 et 360
    //std::cout << longitudes << std::endl; entre -90 et 90

    // Levels normalization
    Eigen::VectorXf levelsNormalized = altitudes.array() - altitudes.minCoeff();
    float levelsNormalizedMin = levelsNormalized.minCoeff();
    float levelsNormalizedMax = levelsNormalized.maxCoeff();
    levelsNormalized /= levelsNormalizedMax;

    V.resize((N + 1) * (N + 1), NUMBER_FIELDS);
    C.resize((N + 1) * (N + 1), NUMBER_FIELDS);
    F.resize(N * N * 2, NUMBER_FIELDS);

    for (int i = 0; i <= N; ++i) {
        for (int j = 0; j <= N; ++j) {
            double theta = 2.0 * M_PI * static_cast<double>(i) / static_cast<double>(N); // entre 0 et 2 PI
            double phi = M_PI * static_cast<double>(j) / static_cast<double>(N); // entre 0 et PI

            V.row(i * (N + 1) + j) << std::sin(phi) * std::cos(theta), std::sin(phi) * std::sin(theta), std::cos(phi);

            // Find the level corresponding to i*(N+1)*j
            level = indicePlusProche(360*theta/(2.0*M_PI), 90 + 180*phi/M_PI, latitudes, longitudes);
            //std::cout << indicePlusProche(360*theta/(2.0*M_PI), 90 + 180*phi/M_PI, latitudes, longitudes) << std::endl;
            //std::cout << latitudes.size() << std::endl;
            //std::cout << theta << " " << phi << std::endl;

            // Colorize the triangle
            color(level, levelsNormalizedMin, (levelsNormalizedMin+levelsNormalizedMax)/2, levelsNormalizedMax, MAX_COLOR, MIN_COLOR, MIN_COLOR, MIN_COLOR, MAX_COLOR, MIN_COLOR, MIN_COLOR, MIN_COLOR, MAX_COLOR, R, G, B);
            C.row(i * (N + 1) + j) << R, G, B;

            /*
            if (phi > 0.1 * M_PI && phi < 0.4 * M_PI && theta > 0.1 * 2.0 * M_PI && theta < 0.4 * 2.0 * M_PI) {
                C.row(i * (N + 1) + j) << MAX_COLOR, MAX_COLOR, MAX_COLOR;
            } else {
                C.row(i * (N + 1) + j) << MAX_COLOR, MIN_COLOR, MIN_COLOR;
            }
             */
        }
    }

    // Populate the face matrix
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            int index = i * N + j;
            F.row(2 * index) << i * (N + 1) + j, (i + 1) * (N + 1) + j, i * (N + 1) + j + 1;
            F.row(2 * index + 1) << (i + 1) * (N + 1) + j, (i + 1) * (N + 1) + j + 1, i * (N + 1) + j + 1;
        }
    }
}


void createLines(const Eigen::VectorXf& altitudes, const Eigen::VectorXf& latitudes, const Eigen::VectorXf& longitudes, Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
    // Get the number of vertices
    int N = altitudes.rows();

    // Populate the vertex matrix
    for (int i = 0; i < N; i++) {
        double theta = longitudes(i) + M_PI/2;
        double phi = latitudes(i);

        double x = sin(phi) * cos(theta);
        double y = sin(phi) * sin(theta);
        double z = cos(phi);

        // Store the spherical coordinates in the vertex matrix
        V(i, ALTITUDE_FIELD_INDEX) = x;
        V(i, LATITUDE_FIELD_INDEX) = y;
        V(i, LONGITUDE_FIELD_INDEX) = z;
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
    viewer.data().point_size = POINT_SIZE;

    // Set the background color
    viewer.core().background_color = Eigen::Vector4f(MIN_COLOR, MIN_COLOR, MIN_COLOR, ALPHA_VALUE);

    // Launch the viewer
    viewer.launch(false,"Example window with lines");
}


void viewMeshColor(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::MatrixXd& C) {

    // Set up the viewer
    igl::opengl::glfw::Viewer viewer;

    // Set the viewer data
    viewer.data().set_mesh(V, F);

    // Set the background color
    viewer.core().background_color = Eigen::Vector4f(MIN_COLOR, MIN_COLOR, MIN_COLOR, ALPHA_VALUE);

    // Set the colors
    viewer.data().set_colors(C);

    // Disable the lines view by default
    viewer.data().show_lines = false;

    // Launch the viewer
    viewer.launch(false,"Example window with lines");
}

// Main program

int main(const int argc, const char *argv[]) {
    try {
        if (argc < NUMBER_ARGS_MIN) {
            throw ErrorCode::NotEnoughArguments;
        }

        // Get the plot type and CSV file name from the command-line arguments
        std::string plotType = argv[PLOT_TYPE_INDEX];
        std::string csvFile = argv[CSV_FILE_INDEX];

        // Read the matrix from the CSV file
        Eigen::MatrixXf dataMat = getMatrixFromCSV(csvFile);

        // Extract the altitudes, latitudes, and longitudes from the matrix
        Eigen::VectorXf altitudes = dataMat.col(ALTITUDE_FIELD_INDEX);
        Eigen::VectorXf latitudes = dataMat.col(LATITUDE_FIELD_INDEX);
        Eigen::VectorXf longitudes = dataMat.col(LONGITUDE_FIELD_INDEX);

        // Extract the size of the data set
        int n = altitudes.size();

        // Create matrices to store the coordinates and colors of the points
        Eigen::MatrixXd V(n, NUMBER_FIELDS);
        Eigen::MatrixXd C(n, NUMBER_FIELDS);

        if (plotType == "-P") {
            setPoints(n, altitudes, latitudes, longitudes, V, C);

            // Display data
            viewPoints(V, C);

        } else {
            // Create a matrix to store triangle connectives
            Eigen::MatrixXi F(n-2,NUMBER_FIELDS);

            if (plotType == "-L") {
                createLines(altitudes, latitudes, longitudes, V, F);
                C.setConstant(MAX_COLOR);
            } else {
                createColorSphere(SPHERE_RESOLUTION, altitudes, latitudes, longitudes, V, F, C);
            }

            // Display data
            viewMeshColor(V,F,C);
        }
    }

    catch(const ErrorCode& errorCode) {
        std::cerr << errorMessages.at(errorCode) << std::endl;
    }
}
