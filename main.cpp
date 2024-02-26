/**
 *  Project description.
 *
 *  @TODO The description
 *
 */
#include <igl/opengl/glfw/Viewer.h>
#include "errorCode.h"
#include "errorMessages.h"
#include "constants.h"

#define TO_RAD (M_PI/180.0)

/**
 * Converts a .csv File into a Eigen Matrix.
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
    /* Checking if the csvFile actually exists */
    if (lineRead.empty()){
        throw ErrorCode::FileNotFound;
    }
    /* Checking if there are enough fields for the display */
    int nbFields = int(std::count(lineRead.begin(), lineRead.end(), ',')) + 1;
    if (nbFields < NUMBER_FIELDS_MIN){
        throw ErrorCode::NotEnoughFields;
    }
    /* Initialisation of the matrix */
    int nbMeasurements = 0;
    Eigen::MatrixXf matrix(nbMeasurements, nbFields);
    /* Number of the current column studied in lineRead */
    int indexDataField;
    /* Read each line */
    while (std::getline(csvData, lineRead)) {
        /* Extend the matrix */
        nbMeasurements++;
        matrix.conservativeResize(nbMeasurements,nbFields);
        indexDataField = 0;
        /* Read the values in the line */
        size_t pos = 0;
        std::string token;
        while ((pos = lineRead.find(',')) != std::string::npos) {
            token = lineRead.substr(0, pos);
            matrix(nbMeasurements-1,indexDataField++) = std::stof(token);
            lineRead.erase(0, pos + 1);
        }
        token = lineRead.substr(0, pos);
        matrix(nbMeasurements-1,indexDataField) = std::stof(token);
        // Verification if null values */
        if (abs(matrix(nbMeasurements-1,indexDataField)) < FLOAT_PRECISION) {
            // Is the data equal to zero at a certain point in time? */
            bool nullData = true;
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
 * Creates points with vertices and faces based on measures according latitude, and longitude.
 *
 * @param n The number of measurements.
 * @param measures, latitudes, longitudes The measures.
 * @param V, C The vertices and color matrices.
 */
void createPoints(const int n, const Eigen::VectorXf& levels, const Eigen::VectorXf& latitudes, const Eigen::VectorXf& longitudes, Eigen::MatrixXd& V, Eigen::MatrixXd& C) {
    /* Cartesian coordinates initialisation */
    double x, y, z;
    double levelMax = levels.maxCoeff();
    double correctedLatitude, correctedLongitude;
    /* Fill the matrices with data */
    for (int i = 0; i < n; ++i) {
        /* Conversion in radians */
        correctedLatitude = latitudes(i)*M_PI/180;
        correctedLongitude = longitudes(i)*M_PI/180;

        x = cos(correctedLatitude) * cos(correctedLongitude);
        y = cos(correctedLatitude) * sin(correctedLongitude);
        z = sin(correctedLatitude);

        V.row(i) = Eigen::Vector3d(x, y, z);

        /* Set colors for visualization */
        C(i, 0) = POINTS_COLOR_R;
        C(i, 1) = POINTS_COLOR_G;
        C(i, 2) = POINTS_COLOR_B;
    }
}

/**
 * Associates a color (R,G,B) to a measure according to the colors associated with the extreme and average values.
 *
 * @param measure The measure itself.
 * @param minV, avgV, maxV The minimal, average and maximal values to identify the color of the measure.
 * @param Rmin, Gmin, Bmin, Ravg, Gavg, Bavg, Rmax, Gmax, Bmax The colors associated with the boundary values.
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
 * Returns the distance between (lat1, lon1) and (lat2, lon2).
 *
 * @param lat1, lat2, lon1, lon2 The different altitudes and longitudes.
 */
float distance(float lat1, float lon1, float lat2, float lon2) {
    float d1 = abs(lat1 - lat2) + abs((lon1 - lon2));
    float d2 = abs(360-lat1 - lat2) + abs((lon1 - lon2));
    float d3 = abs(lat1 - lat2) + abs((90-lon1 - lon2));
    float d4 = abs(360-lat1 - lat2) + abs((90-lon1 - lon2));
    return std::min(std::min(d1,d2),std::min(d3,d4));
}


/**
 * Indicates the nearest index for the (latitude,longitude) coordinates in a given sphere that fits with the latitudes and longitudes data.
 *
 * @param latitude, longitude The latitude and longitude to associate with an index.
 * @param latitudes, longitudes The measures of latitude and longitude.
 *
 * @return The index associated with the best fit.
 */
int nearestIndexforCoords(const float latitude, const float longitude, const Eigen::VectorXf& latitudes, const Eigen::VectorXf& longitudes) {
    float distanceMin = distance(latitude, longitude, latitudes(0), longitudes(0));
    int minIndex = 0;
    for (int i = 1; i < latitudes.size(); i++) {
        float currentDistance = distance(latitude, longitude, latitudes(i), longitudes(i));
        if (currentDistance < distanceMin ) {
            distanceMin = currentDistance;
            minIndex = i;
        }
    }
    /* In that case, we consider the zone has not been observed. */
    if (distanceMin > DISTANCE_MIN_OBSERVATION) {
        return -1;
    }
    return minIndex;
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
    /* Triangle colors initialisation */
    int R, G, B;
    /* Measures normalization */
    Eigen::VectorXf measuresNormalized = measures.array() - measures.minCoeff();
    float measuresNormalizedMax = measuresNormalized.maxCoeff();
    measuresNormalized /= measuresNormalizedMax;
    measuresNormalizedMax = measuresNormalized.maxCoeff();
    float measuresNormalizedMin = measuresNormalized.minCoeff();
    /* Resize the matrices */
    V.resize((N + 1) * (N + 1), NUMBER_FIELDS);
    C.resize((N + 1) * (N + 1), NUMBER_FIELDS);
    F.resize(N * N * 2, NUMBER_FIELDS);

    for (int i = 0; i <= N; ++i) {
        for (int j = 0; j <= N; ++j) {
            double longitude = 2.0 * M_PI * static_cast<double>(j) / static_cast<double>(N);
            double latitude =  M_PI * static_cast<double>(i) / static_cast<double>(N) - M_PI/2;
            /* Generate the triangle */
            V.row(i * (N + 1) + j) << std::cos(latitude) * std::cos(longitude), std::cos(latitude) * std::sin(longitude), std::sin(latitude);
            latitude *= 180 / M_PI;
            longitude *= 180 / M_PI;
            /* Find the measure corresponding to i*(N+1)*j */
            int nearestIndexCoords = nearestIndexforCoords(latitude, longitude, latitudes, longitudes);

            if (nearestIndexCoords < 0) {
                R = R_UNEXPLORED; G = G_UNEXPLORED; B = B_UNEXPLORED; // better : C.row(i * (N + 1) + j) << 0.4,0.4,0.4;
            } else {
                color( measuresNormalized(nearestIndexCoords), measuresNormalizedMin, 0.7, measuresNormalizedMax, R_MIN, G_MIN, B_MIN,
                      R_AVG, G_AVG, B_AVG, R_MAX, G_MAX, B_MAX, R, G, B);
            }
            C.row(i * (N + 1) + j) << R, G, B;
            /* Face matrix */
            if (i < N && j < N) {
                int index = i * N + j;
                F.row(2 * index) << i * (N + 1) + j, (i + 1) * (N + 1) + j, i * (N + 1) + j + 1;
                F.row(2 * index + 1) << (i + 1) * (N + 1) + j, (i + 1) * (N + 1) + j + 1, i * (N + 1) + j + 1;
        }   }
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
    int nbVertices = measures.rows();
    /* Populate the vertex matrix */
    for (int i = 0; i < nbVertices; i++) {
        /* Conversion to radians */
        double correctedLatitude = latitudes(i)*M_PI/180;
        double correctedLongitude = longitudes(i)*M_PI/180;

        double x = cos(correctedLatitude) * cos(correctedLongitude);
        double y = cos(correctedLatitude) * sin(correctedLongitude);
        double z = sin(correctedLatitude);
        /* Store the spherical coordinates in the vertex matrix */
        V(i, VIEW_FIELD_INDEX) = x;
        V(i, LATITUDE_FIELD_INDEX) = y;
        V(i, LONGITUDE_FIELD_INDEX) = z;
    }
    /* Populate the face matrix */
    for (int i = 0; i < nbVertices - 2; ++i) {
        /* Find the three vertices for this triangle */
        int idx1 = i;
        int idx2 = (i + 1) % (nbVertices - 1);
        int idx3 = (i + 2) % (nbVertices - 1) + 1;

        /* Compute the latitudes and longitudes of the three vertices */
        double lat1 = latitudes(idx1);
        double lat2 = latitudes(idx2);
        double lat3 = latitudes(idx3);
        double lon1 = longitudes(idx1);
        double lon2 = longitudes(idx2);
        double lon3 = longitudes(idx3);

        /* Find the nearest vertex for each edge */
        int j1, j2, j3;
        double min_dist, dist;

        /* Edge 1-2 */
        j1 = idx2;
        j2 = idx1;
        min_dist = std::numeric_limits<double>::max();
        for (int j = 0; j < nbVertices; ++j) {
            if (j == idx1 || j == idx2) continue;
            dist = std::abs(latitudes(j) - lat1) + std::abs(longitudes(j) - lon1);
            if (dist < min_dist) {
                min_dist = dist;
                j1 = j;
            }
        }
        /* Edge 2-3 */
        j2 = idx3;
        min_dist = std::numeric_limits<double>::max();
        for (int j = 0; j < nbVertices; ++j) {
            if (j == idx2 || j == idx3) continue;
            dist = std::abs(latitudes(j) - lat2) + std::abs(longitudes(j) - lon2);
            if (dist < min_dist) {
                min_dist = dist;
                j2 = j;
            }
        }
        /* Edge 3-1 */
        j3 = idx1;
        min_dist = std::numeric_limits<double>::max();
        for (int j = 0; j < nbVertices; ++j) {
            if (j == idx3 || j == idx1) continue;
            dist = std::abs(latitudes(j) - lat3) + std::abs(longitudes(j) - lon3);
            if (dist < min_dist) {
                min_dist = dist;
                j3 = j;
            }
        }
        /* Populate the face matrix */
        F.row(i) << idx1, j1, j3;
    }
}

/* Display procedure */

/**
 * Launch the viewer.
 *
 * @param viewer The viewer.
 */
void viewData(igl::opengl::glfw::Viewer viewer){
    viewer.core().background_color = Eigen::Vector4f(BACKGROUND_COLOR_R, BACKGROUND_COLOR_G, BACKGROUND_COLOR_B, ALPHA_VALUE);
    viewer.launch(false,"Data visualization window");
}

/* Main program */

int main(const int argc, const char *argv[]) {
    try {
        if (argc < NUMBER_ARGS_MIN) {
            throw ErrorCode::NotEnoughArguments;
        }
        /* Get the plot type and CSV file name from the command-line arguments */
        std::string plotType = argv[PLOT_TYPE_INDEX];
        std::string csvFile = argv[CSV_FILE_INDEX];

        /* Read the matrix from the CSV file */
        Eigen::MatrixXf dataMat = getMatrixFromCSV(csvFile);

        /* Extract the measures, latitudes, and longitudes from the matrix */
        Eigen::VectorXf measures = dataMat.col(MEASURE_FIELD_INDEX);
        Eigen::VectorXf latitudes = dataMat.col(LATITUDE_FIELD_INDEX);
        Eigen::VectorXf longitudes = dataMat.col(LONGITUDE_FIELD_INDEX);

        /* Extract the size of the data set */
        int n = measures.size();

        /* Create matrices to store the coordinates and colors of the points */
        Eigen::MatrixXd V(n, NUMBER_FIELDS);
        Eigen::MatrixXd C(n, NUMBER_FIELDS);

        /* Viewer initialization */
        igl::opengl::glfw::Viewer viewer;

        if (plotType == "-P") {
            createPoints(n, measures, latitudes, longitudes, V, C);
            /* Set the viewer data */
            viewer.data().set_points(V, C);
            viewer.data().point_size = POINT_SIZE;

        } else {
            /* Create a matrix to store triangle connectives */
            Eigen::MatrixXi F(n-2,NUMBER_FIELDS);
            if (plotType == "-L") {
                createLines(measures, latitudes, longitudes, V, F);
                C.setConstant(MAX_COLOR);
            } else {
                createColorSphere(SPHERE_RESOLUTION, measures, latitudes, longitudes, V, F, C);
            }

            /* Set the viewer data */
            viewer.data().set_mesh(V, F);
            viewer.data().set_colors(C);
            viewer.data().show_lines = false;

            /* To uncomment if you want to show the satellite trajectory */
            /* Eigen::MatrixXd P(n, NUMBER_FIELDS);
            Eigen::MatrixXd C2(n, NUMBER_FIELDS);
            createPoints(n, measures, latitudes, longitudes, P, C2);
            viewer.data().add_points(P,Eigen::RowVector3d(POINTS_COLOR_R,POINTS_COLOR_G,POINTS_COLOR_B));
             */
        }
        /* Launch the viewer */
        viewData(viewer);
    }
    catch(const ErrorCode& errorCode) {
        std::cerr << errorMessages.at(errorCode) << std::endl;
    }
}
