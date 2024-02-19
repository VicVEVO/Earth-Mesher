/**
 *  Project description (bruh)
 *
 *  @TODO Voir unités des données d'entrée pour bien les représenter (aucun sens actuellement)
 *
 */
#include <igl/opengl/glfw/Viewer.h>
#include "errorCode.h"
#include "errorMessages.h"
#include "constants.h"

/**
 * Converts a .csv File into a Eigen Matrix
 *
 * @param csvFile The file -in a Comma-Separated Values format- to convert.
 *
 * @return The eigen matrix corresponding to the .csv file.
 */
Eigen::MatrixXf getMatrixFromCSV(const std::string& csvFile) {
    std::ifstream csvData;
    csvData.open(csvFile);

    std::string lineRead;
    std::getline(csvData, lineRead);

    // Checking if the csvFile actually exists
    if (lineRead.empty()){
        throw ErrorCode::FileNotFound;
    }

    // Checking if there are anough fields for the display
    int nbFields = int(std::count(lineRead.begin(), lineRead.end(), ',')) + 1;
    if (nbFields < NUMBER_FIELDS_MIN){
        throw ErrorCode::NotEnoughFields;
    }

    // Initialisation of the matrix
    int nbMeasurements = 0;
    Eigen::MatrixXf matrix(nbMeasurements, nbFields);

    // Number of the current column studied in lineRead
    int indexDataField;

    // Is the data equal to zero at a certain point in time?
    bool nullData = false;

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
        if (abs(matrix(nbMeasurements-1,indexDataField)) < FLOAT_PRECISION) {
            nullData = true;
            for (int i=0 ; i < nbFields; i++) {
                if (abs(matrix(nbMeasurements-1,i)) > FLOAT_PRECISION) {
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

/**
 * Creates a sphere with vertices and faces based on measures according latitude, and longitude
 *
 * @param n The
 * @param measures, latitudes, longitudes The
 * @param V, C The
 */
void setPoints(const int n, const Eigen::VectorXf& levels, const Eigen::VectorXf& latitudes, const Eigen::VectorXf& longitudes, Eigen::MatrixXd& V, Eigen::MatrixXd& C) {

    double levelMax = levels.maxCoeff();
    double correctedLatitude, correctedLongitude;

    // Fill the matrices with data
    for (int i = 0; i < n; ++i) {

        // Conversion in radians
        correctedLatitude = latitudes(i)*M_PI/180;
        correctedLongitude = longitudes(i)*M_PI/180;

        double x = cos(correctedLatitude) * cos(correctedLongitude);
        double y = cos(correctedLatitude) * sin(correctedLongitude);
        double z = sin(correctedLatitude);

        V.row(i) = Eigen::Vector3d(x, y, z);//(levels(i)/levelMax) * Eigen::Vector3d(x, y, z);

        // Set colors for visualization
        C(i, VIEW_FIELD_INDEX) = MAX_COLOR;
        C(i, LATITUDE_FIELD_INDEX) = MAX_COLOR;
        C(i, LONGITUDE_FIELD_INDEX) = MAX_COLOR;
    }
}

/**
 * Associates a color (R,G,B) to a measure according to the colors associated with the extreme and average values.
 *
 * @param measure The measure itself
 * @param minV, avgV, maxV The minimal, average and maximal values to identify the color of the measure
 * @param Rmin, Gmin, Bmin, Ravg, Gavg, Bavg, Rmax, Gmax, Bmax The colors associated with the boundary values
 */
void color(const float measure, const float minV, const float avgV, const float maxV, const int Rmin, const int Gmin, const int Bmin, const int Ravg, const int Gavg, const int Bavg, const int Rmax, const int Gmax, const int Bmax, int& R, int& G, int& B) {
    if (measure <= minV) {
        R = Rmin;
        G = Gmin;
        B = Bmin;
    } else if (measure <= avgV) {
        float t = (measure - minV) / (avgV - minV);
        R = Rmin + t * (Ravg - Rmin);
        G = Gmin + t * (Gavg - Gmin);
        B = Bmin + t * (Bavg - Bmin);
    } else if (measure <= maxV) {
        float t = (measure - avgV) / (maxV - avgV);
        R = Ravg + t * (Rmax - Ravg);
        G = Gavg + t * (Gmax - Gavg);
        B = Bavg + t * (Bmax - Bavg);
    } else {
        R = Rmax;
        G = Gmax;
        B = Bmax;
    }
}

/**
 * Indicates the nearest index for the (latitude,longitude) coordinates in a given sphere that fits with the latitudes and longitudes data.
 *
 * @param latitude, longitude The latitude and longitude to associate with an index.
 * @param latitudes, longitudes The measures of latitude and longitude.
 *
 * @return The index associated with the best fit
 */
int nearestIndexforCoords(float latitude, float longitude, Eigen::VectorXf& latitudes, Eigen::VectorXf& longitudes) {
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

/**
 * Creates a triangular mesh of a sphere with colored triangles.
 *
 * @param latitude, longitude The latitude and longitude to associate with an index.
 * @param latitudes, longitudes The measures of latitude and longitude.
 *
 * @return The index associated with the best fit.
 */
void createColorSphere(int N, Eigen::VectorXf& measures, Eigen::VectorXf& latitudes, Eigen::VectorXf& longitudes, Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& C) {
    // Triangle colors initialisation
    int R,G,B;

    // Triangle measure initialisation
    float measure;

    // Measures normalization
    Eigen::VectorXf measuresNormalized = measures.array() - measures.minCoeff();
    float measuresNormalizedMax = measuresNormalized.maxCoeff();
    measuresNormalized /= measuresNormalizedMax;
    measuresNormalizedMax = measuresNormalized.maxCoeff();
    float measuresNormalizedMin = measuresNormalized.minCoeff();

    // Resize the matrices
    V.resize((N + 1) * (N + 1), NUMBER_FIELDS);
    C.resize((N + 1) * (N + 1), NUMBER_FIELDS);
    F.resize(N * N * 2, NUMBER_FIELDS);

    for (int i = 0; i <= N; ++i) {
        for (int j = 0; j <= N; ++j) {
            double theta = 2.0 * M_PI * static_cast<double>(i) / static_cast<double>(N); // entre 0 et 2 PI : latitude
            double phi = M_PI * static_cast<double>(j) / static_cast<double>(N); // entre 0 et PI

            // Generate the triangle
            V.row(i * (N + 1) + j) << std::sin(phi) * std::cos(theta), std::sin(phi) * std::sin(theta), std::cos(phi);

            // Find the measure corresponding to i*(N+1)*j
            measure = measuresNormalized(nearestIndexforCoords(360*(0.5-theta/M_PI), 360*phi/M_PI, latitudes, longitudes));

            // Colorize the triangle
            color(measure, measuresNormalizedMin, 0.89, measuresNormalizedMax, MAX_COLOR, MIN_COLOR, MIN_COLOR, MIN_COLOR, MAX_COLOR, MIN_COLOR, MIN_COLOR, MIN_COLOR, MAX_COLOR, R, G, B);
            C.row(i * (N + 1) + j) << R, G, B;
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

/**
 * Creates the line mesh from the measures (points in 3D) by linking them.
 *
 * @param measures, latitudes, longitudes The data triplets that reconstruct 3D points.
 * @param V, F The matrices storing respectively the indices of the rows and the triangle connectivty.
 * @see https://libigl.github.io/tutorial/
 */
void createLines(const Eigen::VectorXf& measures, const Eigen::VectorXf& latitudes, const Eigen::VectorXf& longitudes, Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
    // Get the number of vertices
    int N = measures.rows();

    // Populate the vertex matrix
    for (int i = 0; i < N; i++) {
        // Conversion in radians
        double correctedLatitude = latitudes(i)*M_PI/180;
        double correctedLongitude = longitudes(i)*M_PI/180;

        double x = cos(correctedLatitude) * cos(correctedLongitude);
        double y = cos(correctedLatitude) * sin(correctedLongitude);
        double z = sin(correctedLatitude);

        // Store the spherical coordinates in the vertex matrix
        V(i, VIEW_FIELD_INDEX) = x;
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

/**
 * Display points in 3D.
 *
 * @param V The matrix storing the coordinates of the vertices.
 * @param C The matrix corresponding to the color of each point.
 */
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

/**
 * Display a figure in 3D.
 *
 * @param V, F The matrices storing respectively the indices of the rows and the triangle connectivty.
 * @param C The matrix corresponding to the color of each triangle.
 */
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

        // Extract the measures, latitudes, and longitudes from the matrix
        Eigen::VectorXf measures = dataMat.col(MEASURE_FIELD_INDEX);
        Eigen::VectorXf latitudes = dataMat.col(LATITUDE_FIELD_INDEX);
        Eigen::VectorXf longitudes = dataMat.col(LONGITUDE_FIELD_INDEX);

        // Extract the size of the data set
        int n = measures.size();

        // Create matrices to store the coordinates and colors of the points
        Eigen::MatrixXd V(n, NUMBER_FIELDS);
        Eigen::MatrixXd C(n, NUMBER_FIELDS);

        if (plotType == "-P") {
            setPoints(n, measures, latitudes, longitudes, V, C);

            // Display data
            viewPoints(V, C);

        } else {
            // Create a matrix to store triangle connectives
            Eigen::MatrixXi F(n-2,NUMBER_FIELDS);

            if (plotType == "-L") {
                createLines(measures, latitudes, longitudes, V, F);
                C.setConstant(MAX_COLOR);
            } else {
                createColorSphere(SPHERE_RESOLUTION, measures, latitudes, longitudes, V, F, C);
            }

            // Display data
            viewMeshColor(V,F,C);
        }
    }

    catch(const ErrorCode& errorCode) {
        std::cerr << errorMessages.at(errorCode) << std::endl;
    }
}
