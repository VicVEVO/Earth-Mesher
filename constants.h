#ifndef RENDER_CONSTANTS_H
#define RENDER_CONSTANTS_H

// Precision
const float FLOAT_PRECISION = 1e-6;

// Input settings
const int PLOT_TYPE_INDEX = 1;
const int CSV_FILE_INDEX = 2;
const int NUMBER_ARGS_MIN = 3;

const int NUMBER_FIELDS = 3;
const int NUMBER_FIELDS_MIN = 3;

// Data index
int MEASURE_FIELD_INDEX = 3;
const int VIEW_FIELD_INDEX = 0;
const int LATITUDE_FIELD_INDEX = 1; // Try not to change
const int LONGITUDE_FIELD_INDEX = 2; // Try not to change

// Plot configuration
const int POINT_SIZE = 5;
const int SPHERE_RESOLUTION = 100;
const int ALPHA_VALUE = 1;

const int MIN_COLOR = 0;
const int MAX_COLOR = 255;

#endif //RENDER_CONSTANTS_H
