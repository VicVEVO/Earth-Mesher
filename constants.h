#ifndef RENDER_CONSTANTS_H
#define RENDER_CONSTANTS_H

/* Precision */
const float FLOAT_PRECISION = 1e-6;

/* Input settings */
const int PLOT_TYPE_INDEX = 1;
const int CSV_FILE_INDEX = 2;
const int NUMBER_ARGS_MIN = 3;

const int NUMBER_FIELDS = 3;
const int NUMBER_FIELDS_MIN = 3;

/* Data index */
int MEASURE_FIELD_INDEX = 4;
const int VIEW_FIELD_INDEX = 0; // Try not to change
const int LATITUDE_FIELD_INDEX = 1; // Try not to change
const int LONGITUDE_FIELD_INDEX = 2; // Try not to change

/* Plot configuration */
const int POINT_SIZE = 5;
const int SPHERE_RESOLUTION = 100;
const int ALPHA_VALUE = 1;
const int DISTANCE_MIN_OBSERVATION = 10;

const int MIN_COLOR = 0;
const int MAX_COLOR = 255;

const int BACKGROUND_COLOR_R = 0;
const int BACKGROUND_COLOR_G = 0;
const int BACKGROUND_COLOR_B = 0;

const int POINTS_COLOR_R = 120;
const int POINTS_COLOR_G = 0;
const int POINTS_COLOR_B = 0;

/* Color gradient for the sphere display */
const int R_MIN = 0;
const int G_MIN = 0;
const int B_MIN = 1;

const int R_AVG = 0;
const int G_AVG = 0;
const int B_AVG = 2;

const int R_MAX = 0;
const int G_MAX = 0;
const int B_MAX = 4;

const int R_UNEXPLORED = 1;
const int G_UNEXPLORED = 1;
const int B_UNEXPLORED = 1;

#endif //RENDER_CONSTANTS_H
