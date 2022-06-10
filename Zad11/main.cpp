//
// Created by volie on 27.05.2022.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#ifndef M_PIl
#define M_PIl (3.14159265358979323846264338327950288)
#endif

using namespace std;

const double T_MAX = 2.0;
const double D = 1.0;
const double TAU = 0.1;
const double A = 10.0;
const double ALPHA = 0.0, BETA = 1.0, GAMMA = 0.0, PHI = 0.0, PSI = 1.0, THETA = 0.0;

double **createDoubleMatrix(int rows, int columns);
double *createDoubleVector(int columns);
void destroyMatrix(double **matrix, int rows, int columns);
void destroyVector(double *vector);

void LU_decomposition(double **matrixL, double **matrixU, int *indexes, int n);
double *L_matrixEquation(double **matrixL, double *vectorb, int *indexes, int n);
double *U_matrixEquation(double **matrixU, double *vectory, int *indexes, int n);
bool choice(double **matrixA, int *indexes, int x, int n);
void matrixUfromA(double **matrixA, double **matrixU, int n);
void indexesInput(int *vector, int columns);
void matrixZeros(double **matrix, int rows, int columns);


void Thomass(double *u, double *d, double *l, double *b, double *x, int n);
void get_exact(double start_x, double max_x, double h, double start_t, double max_t, double delta_t);
double exact_value(double x, double t);
double initial_condition(double x);

void ML_Thomass(double start_x, double max_x, double h, double start_t, double max_t, double delta_t);
void MCN_Thomass(double start_x, double max_x, double h, double start_t, double max_t, double delta_t);
void ML_LU(double start_x, double max_x, double h, double start_t, double max_t, double delta_t);
void MCN_LU(double start_x, double max_x, double h, double start_t, double max_t, double delta_t);

void ML_Thomass_error(double start_x, double max_x, double h, double start_t, double max_t, double delta_t);
void MCN_Thomass_error(double start_x, double max_x, double h, double start_t, double max_t, double delta_t);
void ML_LU_error(double start_x, double max_x, double h, double start_t, double max_t, double delta_t);
void MCN_LU_error(double start_x, double max_x, double h, double start_t, double max_t, double delta_t);

double ML_Thomass_maxt(double start_x, double max_x, double h, double start_t, double max_t, double delta_t);
double MCN_Thomass_maxt(double start_x, double max_x, double h, double start_t, double max_t, double delta_t);
double ML_LU_maxt(double start_x, double max_x, double h, double start_t, double max_t, double delta_t);
double MCN_LU_maxt(double start_x, double max_x, double h, double start_t, double max_t, double delta_t);


int main() {
    double start_x = -A, start_t = 0.0, max_x = A, max_t = T_MAX, h = 0.0625, delta_t = h * h, max_t_error;

//    get_exact(start_x, max_x, h, start_t, max_t, delta_t);
//    ML_Thomass(start_x, max_x, h, start_t, max_t, delta_t);
//    MCN_Thomass(start_x, max_x, h, start_t, max_t, delta_t);
//    ML_LU(start_x, max_x, h, start_t, max_t, delta_t);
//    MCN_LU(start_x, max_x, h, start_t, max_t, delta_t);
//    ML_Thomass_error(start_x, max_x, h, start_t, max_t, delta_t);
//    MCN_Thomass_error(start_x, max_x, h, start_t, max_t, delta_t);
//    ML_LU_error(start_x, max_x, h, start_t, max_t, delta_t);
//    MCN_LU_error(start_x, max_x, h, start_t, max_t, delta_t);

//    ofstream ml_thomas_maxt;
//    ml_thomas_maxt.open("ml_thomas_maxt.txt");
//    while(h > 2e-3) {
//        max_t_error = ML_Thomass_maxt(start_x, max_x, h, start_t, max_t, delta_t);
//        ml_thomas_maxt << log10(h) << " " << log10(max_t_error)<< endl;
//        h /= 2;
//        delta_t = h * h;
//    }

//    ofstream mcn_thomas_maxt;
//    mcn_thomas_maxt.open("mcn_thomas_maxt.txt");
//    while(h > 2e-3) {
//        max_t_error = MCN_Thomass_maxt(start_x, max_x, h, start_t, max_t, delta_t);
//        mcn_thomas_maxt << log10(h) << " " << log10(max_t_error)<< endl;
//        h /= 2;
//        delta_t = h * h;
//    }

//    ofstream ml_lu_maxt;
//    ml_lu_maxt.open("ml_lu_maxt.txt");
//    while(h > 5e-3) {
//        max_t_error = ML_LU_maxt(start_x, max_x, h, start_t, max_t, delta_t);
//        ml_lu_maxt << log10(h) << " " << log10(max_t_error)<< endl;
//        h /= 2;
//        delta_t = h * h;
//    }

//    ofstream mcn_lu_maxt;
//    mcn_lu_maxt.open("mcn_lu_maxt.txt");
//    while(h > 5e-3) {
//        max_t_error = ML_LU_maxt(start_x, max_x, h, start_t, max_t, delta_t);
//        mcn_lu_maxt << log10(h) << " " << log10(max_t_error)<< endl;
//        h /= 2;
//        delta_t = h * h;
//    }

    return 0;
}



void get_exact(double start_x, double max_x, double h, double start_t, double max_t, double delta_t) {
    ofstream exact_05;
    exact_05.open("exact_05.txt");
    ofstream exact_1;
    exact_1.open("exact_1.txt");
    ofstream exact_15;
    exact_15.open("exact_15.txt");
    double iter_x = start_x;
    while(iter_x <= max_x) {
        exact_05 << iter_x << " " << exact_value(iter_x, 0.5) << endl;
        iter_x += h;
    }
    iter_x = start_x;
    while(iter_x <= max_x) {
        exact_1 << iter_x << " " << exact_value(iter_x, 1) << endl;
        iter_x += h;
    }
    iter_x = start_x;
    while(iter_x <= max_x) {
        exact_15 << iter_x << " " << exact_value(iter_x, 1.5) << endl;
        iter_x += h;
    }
}

double exact_value(double x, double t) {
    return (1.0 / (2.0 * sqrt(M_PIl * D * (TAU + t)))) * exp(-(x*x) / (4 * D * (TAU + t)));
}

double initial_condition(double x) {
    return (1.0 / (2.0 * sqrt(M_PIl * D * TAU)) * exp(-(x*x) / (4 * D * TAU)));
}

void Thomass(double *u, double *d, double *l, double *b, double *x, int n) {
    double *temp_b = new double[n];
    double *temp_d = new double[n];

    temp_b[0] = b[0];
    temp_d[0] = d[0];
    for(int i = 1; i < n; i++) {
        temp_d[i] = d[i] - ((l[i] * u[i - 1]) / temp_d[i - 1]);
        temp_b[i] = b[i] - ((l[i] * temp_b[i - 1]) / temp_d[i - 1]);
    }

    x[n - 1] = temp_b[n - 1] / temp_d[n - 1];
    for(int i = n - 2; i >= 0; i--)
        x[i] = (temp_b[i] - (u[i] * x[i + 1])) / temp_d[i];

    delete[] temp_d;
    delete[] temp_b;
}


void ML_Thomass(double start_x, double max_x, double h, double start_t, double max_t, double delta_t) {
    int x_len = int((max_x - start_x) / h) + 1;
    double **matrix_A = new double *[3];

    for (int i = 0; i < 3; i++) {
        matrix_A[i] = new double[x_len];
    }
    double *vector_b = new double[x_len];
    double *vector_u = new double[x_len];

    for (int i = 0; i < x_len; i++) {
        vector_b[i] = initial_condition(start_x + h * i) * (-1.0);
    }

    vector_b[0] = 0.0;
    vector_b[x_len - 1] = 0.0;
    vector_u[0] = 0.0;
    vector_u[x_len - 1] = 0.0;

    matrix_A[0][0] = 0.0;
    matrix_A[1][0] = 1.0;
    matrix_A[2][0] = 0.0;
    matrix_A[0][x_len - 1] = 0.0;
    matrix_A[1][x_len - 1] = 1.0;
    matrix_A[2][x_len - 1] = 0.0;
    for (int i = 1; i < x_len - 1; i++) {
        matrix_A[0][i] = 1.0;
        matrix_A[1][i] = -3.0;
        matrix_A[2][i] = 1.0;
    }

    ofstream ml_thomas_05;
    ml_thomas_05.open("ml_thomas_05.txt");
    ofstream ml_thomas_1;
    ml_thomas_1.open("ml_thomas_1.txt");
    ofstream ml_thomas_15;
    ml_thomas_15.open("ml_thomas_15.txt");

    while(start_t < max_t) {
        Thomass(matrix_A[0], matrix_A[1], matrix_A[2], vector_b, vector_u, x_len);

        if(start_t == 0.5) {
            for(int i = 0; i < x_len; i++) {
                ml_thomas_05 << start_x + i * h << " " << vector_u[i] << endl;
            }
        }
        if(start_t == 1) {
            for(int i = 0; i < x_len; i++) {
                ml_thomas_1 << start_x + i * h << " " << vector_u[i] << endl;
            }
        }
        if(start_t == 1.5) {
            for(int i = 0; i < x_len; i++) {
                ml_thomas_15 << start_x + i * h << " " << vector_u[i] << endl;
            }
        }

        for (int i = 0; i < x_len; i++) {
            vector_b[i] = vector_u[i] * (-1.0);
        }

        vector_b[0] = 0.0;
        vector_b[x_len - 1] = 0.0;
        vector_u[0] = 0.0;
        vector_u[x_len - 1] = 0.0;
        start_t += delta_t;
    }

    destroyMatrix(matrix_A, 3, x_len);
    destroyVector(vector_u);
    destroyVector(vector_b);
}

void MCN_Thomass(double start_x, double max_x, double h, double start_t, double max_t, double delta_t) {
    int x_len = int((max_x - start_x) / h) + 1;
    double **matrix_A = new double *[3];

    for (int i = 0; i < 3; i++) {
        matrix_A[i] = new double[x_len];
    }
    double *vector_b = new double[x_len];
    double *vector_u = new double[x_len];

    for (int i = 0; i < x_len; i++) {
        vector_b[i] = initial_condition(start_x + h * i) * (-1.0);
    }

    vector_b[0] = 0.0;
    vector_b[x_len - 1] = 0.0;
    vector_u[0] = 0.0;
    vector_u[x_len - 1] = 0.0;

    matrix_A[0][0] = 0.0;
    matrix_A[1][0] = 1.0;
    matrix_A[2][0] = 0.0;
    matrix_A[0][x_len - 1] = 0.0;
    matrix_A[1][x_len - 1] = 1.0;
    matrix_A[2][x_len - 1] = 0.0;
    for (int i = 1; i < x_len - 1; i++) {
        matrix_A[0][i] = 0.5;
        matrix_A[1][i] = -2.0;
        matrix_A[2][i] = 0.5;
    }

    ofstream mcn_thomas_05;
    mcn_thomas_05.open("mcn_thomas_05.txt");
    ofstream mcn_thomas_1;
    mcn_thomas_1.open("mcn_thomas_1.txt");
    ofstream mcn_thomas_15;
    mcn_thomas_15.open("mcn_thomas_15.txt");

    while(start_t < max_t) {
        Thomass(matrix_A[0], matrix_A[1], matrix_A[2], vector_b, vector_u, x_len);

        if(start_t == 0.5) {
            for(int i = 0; i < x_len; i++) {
                mcn_thomas_05 << start_x + i * h << " " << vector_u[i] << endl;
            }
        }
        if(start_t == 1) {
            for(int i = 0; i < x_len; i++) {
                mcn_thomas_1 << start_x + i * h << " " << vector_u[i] << endl;
            }
        }
        if(start_t == 1.5) {
            for(int i = 0; i < x_len; i++) {
                mcn_thomas_15 << start_x + i * h << " " << vector_u[i] << endl;
            }
        }

        for (int i = 1; i < x_len - 1; i++) {
            vector_b[i] = (0.5 * vector_u[i - 1] + 0.5 * vector_u[i + 1]) * (-1.0);
        }

        vector_b[0] = 0.0;
        vector_b[x_len - 1] = 0.0;
        vector_u[0] = 0.0;
        vector_u[x_len - 1] = 0.0;
        start_t += delta_t;
    }

    destroyMatrix(matrix_A, 3, x_len);
    destroyVector(vector_u);
    destroyVector(vector_b);
}

void ML_LU(double start_x, double max_x, double h, double start_t, double max_t, double delta_t) {
    int x_len = int((max_x - start_x) / h) + 1;

    double **matrix_A = createDoubleMatrix(x_len, x_len);
    double **matrix_L = createDoubleMatrix(x_len, x_len);
    double **matrix_U = createDoubleMatrix(x_len, x_len);
    matrixZeros(matrix_A, x_len, x_len);
    matrixZeros(matrix_L, x_len, x_len);

    double *vector_b = new double[x_len];
    double *vector_u = new double[x_len];
    double *vector_y = new double[x_len];
    int *indexes = new int[x_len];
    indexesInput(indexes, x_len);

    for (int i = 0; i < x_len; i++) {
        vector_b[i] = initial_condition(start_x + h * i) * (-1.0);
    }

    vector_b[0] = 0.0;
    vector_b[x_len - 1] = 0.0;
    vector_u[0] = 0.0;
    vector_u[x_len - 1] = 0.0;

    matrix_A[0][0] = 1.0;
    matrix_A[x_len - 1][x_len - 1] = 1.0;
    for(int i = 1; i < x_len - 1; i++) {
        matrix_A[i][i - 1] = 1.0;
        matrix_A[i][i] = -3.0;
        matrix_A[i][i + 1] = 1.0;
    }
    matrixUfromA(matrix_A, matrix_U, x_len);

    ofstream ml_lu_05;
    ml_lu_05.open("ml_lu_05.txt");
    ofstream ml_lu_1;
    ml_lu_1.open("ml_lu_1.txt");
    ofstream ml_lu_15;
    ml_lu_15.open("ml_lu_15.txt");

    while(start_t < max_t) {
        LU_decomposition(matrix_L, matrix_U, indexes, x_len);
        vector_y = L_matrixEquation(matrix_L, vector_b, indexes, x_len);
        vector_u = U_matrixEquation(matrix_U, vector_y, indexes, x_len);

        if(start_t == 0.5) {
            for(int i = 0; i < x_len; i++) {
                ml_lu_05 << start_x + i * h << " " << vector_u[i] << endl;
            }
        }
        if(start_t == 1) {
            for(int i = 0; i < x_len; i++) {
                ml_lu_1 << start_x + i * h << " " << vector_u[i] << endl;
            }
        }
        if(start_t == 1.5) {
            for(int i = 0; i < x_len; i++) {
                ml_lu_15 << start_x + i * h << " " << vector_u[i] << endl;
            }
        }

        for (int i = 0; i < x_len; i++) {
            vector_b[i] = vector_u[i] * (-1.0);
        }

        vector_b[0] = 0.0;
        vector_b[x_len - 1] = 0.0;
        vector_u[0] = 0.0;
        vector_u[x_len - 1] = 0.0;

        indexesInput(indexes, x_len);

        matrixZeros(matrix_A, x_len, x_len);
        matrixZeros(matrix_L, x_len, x_len);
        matrix_A[0][0] = 1.0;
        matrix_A[x_len - 1][x_len - 1] = 1.0;
        for(int i = 1; i < x_len - 1; i++) {
            matrix_A[i][i - 1] = 1.0;
            matrix_A[i][i] = -3.0;
            matrix_A[i][i + 1] = 1.0;
        }
        matrixUfromA(matrix_A, matrix_U, x_len);

        start_t += delta_t;
    }

    destroyMatrix(matrix_A, x_len, x_len);
    destroyMatrix(matrix_L, x_len, x_len);
    destroyMatrix(matrix_U, x_len, x_len);
    destroyVector(vector_b);
    destroyVector(vector_u);
    destroyVector(vector_y);
    delete[] indexes;
}

void MCN_LU(double start_x, double max_x, double h, double start_t, double max_t, double delta_t) {
    int x_len = int((max_x - start_x) / h) + 1;

    double **matrix_A = createDoubleMatrix(x_len, x_len);
    double **matrix_L = createDoubleMatrix(x_len, x_len);
    double **matrix_U = createDoubleMatrix(x_len, x_len);
    matrixZeros(matrix_A, x_len, x_len);
    matrixZeros(matrix_L, x_len, x_len);

    double *vector_b = new double[x_len];
    double *vector_u = new double[x_len];
    double *vector_y = new double[x_len];
    int *indexes = new int[x_len];
    indexesInput(indexes, x_len);

    for (int i = 0; i < x_len; i++) {
        vector_b[i] = initial_condition(start_x + h * i) * (-1.0);
    }

    vector_b[0] = 0.0;
    vector_b[x_len - 1] = 0.0;
    vector_u[0] = 0.0;
    vector_u[x_len - 1] = 0.0;

    matrix_A[0][0] = 1.0;
    matrix_A[x_len - 1][x_len - 1] = 1.0;
    for(int i = 1; i < x_len - 1; i++) {
        matrix_A[i][i - 1] = 0.5;
        matrix_A[i][i] = -2.0;
        matrix_A[i][i + 1] = 0.5;
    }
    matrixUfromA(matrix_A, matrix_U, x_len);

    ofstream mcn_lu_05;
    mcn_lu_05.open("mcn_lu_05.txt");
    ofstream mcn_lu_1;
    mcn_lu_1.open("mcn_lu_1.txt");
    ofstream mcn_lu_15;
    mcn_lu_15.open("mcn_lu_15.txt");

    while(start_t < max_t) {
        LU_decomposition(matrix_L, matrix_U, indexes, x_len);
        vector_y = L_matrixEquation(matrix_L, vector_b, indexes, x_len);
        vector_u = U_matrixEquation(matrix_U, vector_y, indexes, x_len);

        if(start_t == 0.5) {
            for(int i = 0; i < x_len; i++) {
                mcn_lu_05 << start_x + i * h << " " << vector_u[i] << endl;
            }
        }
        if(start_t == 1) {
            for(int i = 0; i < x_len; i++) {
                mcn_lu_1 << start_x + i * h << " " << vector_u[i] << endl;
            }
        }
        if(start_t == 1.5) {
            for(int i = 0; i < x_len; i++) {
                mcn_lu_15 << start_x + i * h << " " << vector_u[i] << endl;
            }
        }

        for (int i = 0; i < x_len; i++) {
            vector_b[i] = (0.5 * vector_u[i - 1] + 0.5 * vector_u[i + 1]) * (-1.0);
        }

        vector_b[0] = 0.0;
        vector_b[x_len - 1] = 0.0;
        vector_u[0] = 0.0;
        vector_u[x_len - 1] = 0.0;

        indexesInput(indexes, x_len);

        matrixZeros(matrix_A, x_len, x_len);
        matrixZeros(matrix_L, x_len, x_len);
        matrix_A[0][0] = 1.0;
        matrix_A[x_len - 1][x_len - 1] = 1.0;
        for(int i = 1; i < x_len - 1; i++) {
            matrix_A[i][i - 1] = 0.5;
            matrix_A[i][i] = -2.0;
            matrix_A[i][i + 1] = 0.5;
        }
        matrixUfromA(matrix_A, matrix_U, x_len);

        start_t += delta_t;
    }

    destroyMatrix(matrix_A, x_len, x_len);
    destroyMatrix(matrix_L, x_len, x_len);
    destroyMatrix(matrix_U, x_len, x_len);
    destroyVector(vector_b);
    destroyVector(vector_u);
    destroyVector(vector_y);
    delete[] indexes;
}



void ML_Thomass_error(double start_x, double max_x, double h, double start_t, double max_t, double delta_t) {
    double max_error = 0.0;
    double error;
    int x_len = int((max_x - start_x) / h) + 1;
    double **matrix_A = new double *[3];

    for (int i = 0; i < 3; i++) {
        matrix_A[i] = new double[x_len];
    }
    double *vector_b = new double[x_len];
    double *vector_u = new double[x_len];

    for (int i = 0; i < x_len; i++) {
        vector_b[i] = initial_condition(start_x + h * i) * (-1.0);
    }

    vector_b[0] = 0.0;
    vector_b[x_len - 1] = 0.0;
    vector_u[0] = 0.0;
    vector_u[x_len - 1] = 0.0;

    matrix_A[0][0] = 0.0;
    matrix_A[1][0] = 1.0;
    matrix_A[2][0] = 0.0;
    matrix_A[0][x_len - 1] = 0.0;
    matrix_A[1][x_len - 1] = 1.0;
    matrix_A[2][x_len - 1] = 0.0;
    for (int i = 1; i < x_len - 1; i++) {
        matrix_A[0][i] = 1.0;
        matrix_A[1][i] = -3.0;
        matrix_A[2][i] = 1.0;
    }

    ofstream ml_thomas_error;
    ml_thomas_error.open("ml_thomas_error.txt");

    while(start_t < max_t) {
        Thomass(matrix_A[0], matrix_A[1], matrix_A[2], vector_b, vector_u, x_len);

        for(int i = 0; i < x_len; i++) {
            error = fabs(vector_u[i] - exact_value(start_x + i * h, start_t));
            if(error > max_error)
                max_error = error;
        }
        ml_thomas_error << start_t << " " << max_error << endl;
        max_error = 0.0;

        for (int i = 0; i < x_len; i++) {
            vector_b[i] = vector_u[i] * (-1.0);
        }

        vector_b[0] = 0.0;
        vector_b[x_len - 1] = 0.0;
        vector_u[0] = 0.0;
        vector_u[x_len - 1] = 0.0;
        start_t += delta_t;
    }

    destroyMatrix(matrix_A, 3, x_len);
    destroyVector(vector_u);
    destroyVector(vector_b);
}

void MCN_Thomass_error(double start_x, double max_x, double h, double start_t, double max_t, double delta_t) {
    double max_error = 0.0;
    double error;
    int x_len = int((max_x - start_x) / h) + 1;
    double **matrix_A = new double *[3];

    for (int i = 0; i < 3; i++) {
        matrix_A[i] = new double[x_len];
    }
    double *vector_b = new double[x_len];
    double *vector_u = new double[x_len];

    for (int i = 0; i < x_len; i++) {
        vector_b[i] = initial_condition(start_x + h * i) * (-1.0);
    }

    vector_b[0] = 0.0;
    vector_b[x_len - 1] = 0.0;
    vector_u[0] = 0.0;
    vector_u[x_len - 1] = 0.0;

    matrix_A[0][0] = 0.0;
    matrix_A[1][0] = 1.0;
    matrix_A[2][0] = 0.0;
    matrix_A[0][x_len - 1] = 0.0;
    matrix_A[1][x_len - 1] = 1.0;
    matrix_A[2][x_len - 1] = 0.0;
    for (int i = 1; i < x_len - 1; i++) {
        matrix_A[0][i] = 0.5;
        matrix_A[1][i] = -2.0;
        matrix_A[2][i] = 0.5;
    }

    ofstream mcn_thomas_error;
    mcn_thomas_error.open("mcn_thomas_error.txt");

    while(start_t < max_t) {
        Thomass(matrix_A[0], matrix_A[1], matrix_A[2], vector_b, vector_u, x_len);

        for(int i = 0; i < x_len; i++) {
            error = fabs(vector_u[i] - exact_value(start_x + i * h, start_t));
            if(error > max_error)
                max_error = error;
        }
        mcn_thomas_error << start_t << " " << max_error << endl;
        max_error = 0.0;

        for (int i = 1; i < x_len - 1; i++) {
            vector_b[i] = (0.5 * vector_u[i - 1] + 0.5 * vector_u[i + 1]) * (-1.0);
        }

        vector_b[0] = 0.0;
        vector_b[x_len - 1] = 0.0;
        vector_u[0] = 0.0;
        vector_u[x_len - 1] = 0.0;
        start_t += delta_t;
    }

    destroyMatrix(matrix_A, 3, x_len);
    destroyVector(vector_u);
    destroyVector(vector_b);
}

void ML_LU_error(double start_x, double max_x, double h, double start_t, double max_t, double delta_t) {
    double max_error = 0.0;
    double error;
    int x_len = int((max_x - start_x) / h) + 1;

    double **matrix_A = createDoubleMatrix(x_len, x_len);
    double **matrix_L = createDoubleMatrix(x_len, x_len);
    double **matrix_U = createDoubleMatrix(x_len, x_len);
    matrixZeros(matrix_A, x_len, x_len);
    matrixZeros(matrix_L, x_len, x_len);

    double *vector_b = new double[x_len];
    double *vector_u = new double[x_len];
    double *vector_y = new double[x_len];
    int *indexes = new int[x_len];
    indexesInput(indexes, x_len);

    for (int i = 0; i < x_len; i++) {
        vector_b[i] = initial_condition(start_x + h * i) * (-1.0);
    }

    vector_b[0] = 0.0;
    vector_b[x_len - 1] = 0.0;
    vector_u[0] = 0.0;
    vector_u[x_len - 1] = 0.0;

    matrix_A[0][0] = 1.0;
    matrix_A[x_len - 1][x_len - 1] = 1.0;
    for(int i = 1; i < x_len - 1; i++) {
        matrix_A[i][i - 1] = 1.0;
        matrix_A[i][i] = -3.0;
        matrix_A[i][i + 1] = 1.0;
    }
    matrixUfromA(matrix_A, matrix_U, x_len);

    ofstream ml_lu_error;
    ml_lu_error.open("ml_lu_error.txt");

    while(start_t < max_t) {
        LU_decomposition(matrix_L, matrix_U, indexes, x_len);
        vector_y = L_matrixEquation(matrix_L, vector_b, indexes, x_len);
        vector_u = U_matrixEquation(matrix_U, vector_y, indexes, x_len);

        for(int i = 0; i < x_len; i++) {
            error = fabs(vector_u[i] - exact_value(start_x + i * h, start_t));
            if(error > max_error)
                max_error = error;
        }
        ml_lu_error << start_t << " " << max_error << endl;
        max_error = 0.0;

        for (int i = 0; i < x_len; i++) {
            vector_b[i] = vector_u[i] * (-1.0);
        }

        vector_b[0] = 0.0;
        vector_b[x_len - 1] = 0.0;
        vector_u[0] = 0.0;
        vector_u[x_len - 1] = 0.0;

        indexesInput(indexes, x_len);

        matrixZeros(matrix_A, x_len, x_len);
        matrixZeros(matrix_L, x_len, x_len);
        matrix_A[0][0] = 1.0;
        matrix_A[x_len - 1][x_len - 1] = 1.0;
        for(int i = 1; i < x_len - 1; i++) {
            matrix_A[i][i - 1] = 1.0;
            matrix_A[i][i] = -3.0;
            matrix_A[i][i + 1] = 1.0;
        }
        matrixUfromA(matrix_A, matrix_U, x_len);

        start_t += delta_t;
    }

    destroyMatrix(matrix_A, x_len, x_len);
    destroyMatrix(matrix_L, x_len, x_len);
    destroyMatrix(matrix_U, x_len, x_len);
    destroyVector(vector_b);
    destroyVector(vector_u);
    destroyVector(vector_y);
    delete[] indexes;
}

void MCN_LU_error(double start_x, double max_x, double h, double start_t, double max_t, double delta_t) {
    double max_error = 0.0;
    double error;
    int x_len = int((max_x - start_x) / h) + 1;

    double **matrix_A = createDoubleMatrix(x_len, x_len);
    double **matrix_L = createDoubleMatrix(x_len, x_len);
    double **matrix_U = createDoubleMatrix(x_len, x_len);
    matrixZeros(matrix_A, x_len, x_len);
    matrixZeros(matrix_L, x_len, x_len);

    double *vector_b = new double[x_len];
    double *vector_u = new double[x_len];
    double *vector_y = new double[x_len];
    int *indexes = new int[x_len];
    indexesInput(indexes, x_len);

    for (int i = 0; i < x_len; i++) {
        vector_b[i] = initial_condition(start_x + h * i) * (-1.0);
    }

    vector_b[0] = 0.0;
    vector_b[x_len - 1] = 0.0;
    vector_u[0] = 0.0;
    vector_u[x_len - 1] = 0.0;

    matrix_A[0][0] = 1.0;
    matrix_A[x_len - 1][x_len - 1] = 1.0;
    for(int i = 1; i < x_len - 1; i++) {
        matrix_A[i][i - 1] = 0.5;
        matrix_A[i][i] = -2.0;
        matrix_A[i][i + 1] = 0.5;
    }
    matrixUfromA(matrix_A, matrix_U, x_len);

    ofstream mcn_lu_error;
    mcn_lu_error.open("mcn_lu_error.txt");

    while(start_t < max_t) {
        LU_decomposition(matrix_L, matrix_U, indexes, x_len);
        vector_y = L_matrixEquation(matrix_L, vector_b, indexes, x_len);
        vector_u = U_matrixEquation(matrix_U, vector_y, indexes, x_len);

        for(int i = 0; i < x_len; i++) {
            error = fabs(vector_u[i] - exact_value(start_x + i * h, start_t));
            if(error > max_error)
                max_error = error;
        }
        mcn_lu_error << start_t << " " << max_error << endl;
        max_error = 0.0;

        for (int i = 0; i < x_len; i++) {
            vector_b[i] = (0.5 * vector_u[i - 1] + 0.5 * vector_u[i + 1]) * (-1.0);
        }

        vector_b[0] = 0.0;
        vector_b[x_len - 1] = 0.0;
        vector_u[0] = 0.0;
        vector_u[x_len - 1] = 0.0;

        indexesInput(indexes, x_len);

        matrixZeros(matrix_A, x_len, x_len);
        matrixZeros(matrix_L, x_len, x_len);
        matrix_A[0][0] = 1.0;
        matrix_A[x_len - 1][x_len - 1] = 1.0;
        for(int i = 1; i < x_len - 1; i++) {
            matrix_A[i][i - 1] = 0.5;
            matrix_A[i][i] = -2.0;
            matrix_A[i][i + 1] = 0.5;
        }
        matrixUfromA(matrix_A, matrix_U, x_len);

        start_t += delta_t;
    }

    destroyMatrix(matrix_A, x_len, x_len);
    destroyMatrix(matrix_L, x_len, x_len);
    destroyMatrix(matrix_U, x_len, x_len);
    destroyVector(vector_b);
    destroyVector(vector_u);
    destroyVector(vector_y);
    delete[] indexes;
}



double ML_Thomass_maxt(double start_x, double max_x, double h, double start_t, double max_t, double delta_t) {
    double max_error = 0.0;
    double error;
    int x_len = int((max_x - start_x) / h) + 1;
    double **matrix_A = new double *[3];

    for (int i = 0; i < 3; i++) {
        matrix_A[i] = new double[x_len];
    }
    double *vector_b = new double[x_len];
    double *vector_u = new double[x_len];

    for (int i = 0; i < x_len; i++) {
        vector_b[i] = initial_condition(start_x + h * i) * (-1.0);
    }

    vector_b[0] = 0.0;
    vector_b[x_len - 1] = 0.0;
    vector_u[0] = 0.0;
    vector_u[x_len - 1] = 0.0;

    matrix_A[0][0] = 0.0;
    matrix_A[1][0] = 1.0;
    matrix_A[2][0] = 0.0;
    matrix_A[0][x_len - 1] = 0.0;
    matrix_A[1][x_len - 1] = 1.0;
    matrix_A[2][x_len - 1] = 0.0;
    for (int i = 1; i < x_len - 1; i++) {
        matrix_A[0][i] = 1.0;
        matrix_A[1][i] = -3.0;
        matrix_A[2][i] = 1.0;
    }

    while(start_t < max_t) {
        Thomass(matrix_A[0], matrix_A[1], matrix_A[2], vector_b, vector_u, x_len);

        if(start_t + delta_t >= max_t) {
            for(int i = 0; i < x_len; i++) {
                error = fabs(vector_u[i] - exact_value(start_x + i * h, start_t));
                if(error > max_error)
                    max_error = error;
            }
        }

        for (int i = 0; i < x_len; i++) {
            vector_b[i] = vector_u[i] * (-1.0);
        }

        vector_b[0] = 0.0;
        vector_b[x_len - 1] = 0.0;
        vector_u[0] = 0.0;
        vector_u[x_len - 1] = 0.0;
        start_t += delta_t;
    }

    destroyMatrix(matrix_A, 3, x_len);
    destroyVector(vector_u);
    destroyVector(vector_b);

    return max_error;
}

double MCN_Thomass_maxt(double start_x, double max_x, double h, double start_t, double max_t, double delta_t) {
    double max_error = 0.0;
    double error;
    int x_len = int((max_x - start_x) / h) + 1;
    double **matrix_A = new double *[3];

    for (int i = 0; i < 3; i++) {
        matrix_A[i] = new double[x_len];
    }
    double *vector_b = new double[x_len];
    double *vector_u = new double[x_len];

    for (int i = 0; i < x_len; i++) {
        vector_b[i] = initial_condition(start_x + h * i) * (-1.0);
    }

    vector_b[0] = 0.0;
    vector_b[x_len - 1] = 0.0;
    vector_u[0] = 0.0;
    vector_u[x_len - 1] = 0.0;

    matrix_A[0][0] = 0.0;
    matrix_A[1][0] = 1.0;
    matrix_A[2][0] = 0.0;
    matrix_A[0][x_len - 1] = 0.0;
    matrix_A[1][x_len - 1] = 1.0;
    matrix_A[2][x_len - 1] = 0.0;
    for (int i = 1; i < x_len - 1; i++) {
        matrix_A[0][i] = 0.5;
        matrix_A[1][i] = -2.0;
        matrix_A[2][i] = 0.5;
    }

    while(start_t < max_t) {
        Thomass(matrix_A[0], matrix_A[1], matrix_A[2], vector_b, vector_u, x_len);

        if(start_t + delta_t >= max_t) {
            for(int i = 0; i < x_len; i++) {
                error = fabs(vector_u[i] - exact_value(start_x + i * h, start_t));
                if(error > max_error)
                    max_error = error;
            }
        }

        for (int i = 1; i < x_len - 1; i++) {
            vector_b[i] = (0.5 * vector_u[i - 1] + 0.5 * vector_u[i + 1]) * (-1.0);
        }

        vector_b[0] = 0.0;
        vector_b[x_len - 1] = 0.0;
        vector_u[0] = 0.0;
        vector_u[x_len - 1] = 0.0;
        start_t += delta_t;
    }

    destroyMatrix(matrix_A, 3, x_len);
    destroyVector(vector_u);
    destroyVector(vector_b);

    return max_error;
}

double ML_LU_maxt(double start_x, double max_x, double h, double start_t, double max_t, double delta_t) {
    double max_error = 0.0;
    double error;
    int x_len = int((max_x - start_x) / h) + 1;

    double **matrix_A = createDoubleMatrix(x_len, x_len);
    double **matrix_L = createDoubleMatrix(x_len, x_len);
    double **matrix_U = createDoubleMatrix(x_len, x_len);
    matrixZeros(matrix_A, x_len, x_len);
    matrixZeros(matrix_L, x_len, x_len);

    double *vector_b = new double[x_len];
    double *vector_u = new double[x_len];
    double *vector_y = new double[x_len];
    int *indexes = new int[x_len];
    indexesInput(indexes, x_len);

    for (int i = 0; i < x_len; i++) {
        vector_b[i] = initial_condition(start_x + h * i) * (-1.0);
    }

    vector_b[0] = 0.0;
    vector_b[x_len - 1] = 0.0;
    vector_u[0] = 0.0;
    vector_u[x_len - 1] = 0.0;

    matrix_A[0][0] = 1.0;
    matrix_A[x_len - 1][x_len - 1] = 1.0;
    for(int i = 1; i < x_len - 1; i++) {
        matrix_A[i][i - 1] = 1.0;
        matrix_A[i][i] = -3.0;
        matrix_A[i][i + 1] = 1.0;
    }
    matrixUfromA(matrix_A, matrix_U, x_len);

    while(start_t < max_t) {
        LU_decomposition(matrix_L, matrix_U, indexes, x_len);
        vector_y = L_matrixEquation(matrix_L, vector_b, indexes, x_len);
        vector_u = U_matrixEquation(matrix_U, vector_y, indexes, x_len);

        if(start_t + delta_t >= max_t) {
            for(int i = 0; i < x_len; i++) {
                error = fabs(vector_u[i] - exact_value(start_x + i * h, start_t));
                if(error > max_error)
                    max_error = error;
            }
        }

        for (int i = 0; i < x_len; i++) {
            vector_b[i] = vector_u[i] * (-1.0);
        }

        vector_b[0] = 0.0;
        vector_b[x_len - 1] = 0.0;
        vector_u[0] = 0.0;
        vector_u[x_len - 1] = 0.0;

        indexesInput(indexes, x_len);

        matrixZeros(matrix_A, x_len, x_len);
        matrixZeros(matrix_L, x_len, x_len);
        matrix_A[0][0] = 1.0;
        matrix_A[x_len - 1][x_len - 1] = 1.0;
        for(int i = 1; i < x_len - 1; i++) {
            matrix_A[i][i - 1] = 1.0;
            matrix_A[i][i] = -3.0;
            matrix_A[i][i + 1] = 1.0;
        }
        matrixUfromA(matrix_A, matrix_U, x_len);

        start_t += delta_t;
    }

    destroyMatrix(matrix_A, x_len, x_len);
    destroyMatrix(matrix_L, x_len, x_len);
    destroyMatrix(matrix_U, x_len, x_len);
    destroyVector(vector_b);
    destroyVector(vector_u);
    destroyVector(vector_y);
    delete[] indexes;

    return max_error;
}

double MCN_LU_maxt(double start_x, double max_x, double h, double start_t, double max_t, double delta_t) {
    double max_error = 0.0;
    double error;
    int x_len = int((max_x - start_x) / h) + 1;

    double **matrix_A = createDoubleMatrix(x_len, x_len);
    double **matrix_L = createDoubleMatrix(x_len, x_len);
    double **matrix_U = createDoubleMatrix(x_len, x_len);
    matrixZeros(matrix_A, x_len, x_len);
    matrixZeros(matrix_L, x_len, x_len);

    double *vector_b = new double[x_len];
    double *vector_u = new double[x_len];
    double *vector_y = new double[x_len];
    int *indexes = new int[x_len];
    indexesInput(indexes, x_len);

    for (int i = 0; i < x_len; i++) {
        vector_b[i] = initial_condition(start_x + h * i) * (-1.0);
    }

    vector_b[0] = 0.0;
    vector_b[x_len - 1] = 0.0;
    vector_u[0] = 0.0;
    vector_u[x_len - 1] = 0.0;

    matrix_A[0][0] = 1.0;
    matrix_A[x_len - 1][x_len - 1] = 1.0;
    for(int i = 1; i < x_len - 1; i++) {
        matrix_A[i][i - 1] = 0.5;
        matrix_A[i][i] = -2.0;
        matrix_A[i][i + 1] = 0.5;
    }
    matrixUfromA(matrix_A, matrix_U, x_len);

    while(start_t < max_t) {
        LU_decomposition(matrix_L, matrix_U, indexes, x_len);
        vector_y = L_matrixEquation(matrix_L, vector_b, indexes, x_len);
        vector_u = U_matrixEquation(matrix_U, vector_y, indexes, x_len);

        if(start_t + delta_t >= max_t) {
            for(int i = 0; i < x_len; i++) {
                error = fabs(vector_u[i] - exact_value(start_x + i * h, start_t));
                if(error > max_error)
                    max_error = error;
            }
        }

        for (int i = 0; i < x_len; i++) {
            vector_b[i] = (0.5 * vector_u[i - 1] + 0.5 * vector_u[i + 1]) * (-1.0);
        }

        vector_b[0] = 0.0;
        vector_b[x_len - 1] = 0.0;
        vector_u[0] = 0.0;
        vector_u[x_len - 1] = 0.0;

        indexesInput(indexes, x_len);

        matrixZeros(matrix_A, x_len, x_len);
        matrixZeros(matrix_L, x_len, x_len);
        matrix_A[0][0] = 1.0;
        matrix_A[x_len - 1][x_len - 1] = 1.0;
        for(int i = 1; i < x_len - 1; i++) {
            matrix_A[i][i - 1] = 0.5;
            matrix_A[i][i] = -2.0;
            matrix_A[i][i + 1] = 0.5;
        }
        matrixUfromA(matrix_A, matrix_U, x_len);

        start_t += delta_t;
    }

    destroyMatrix(matrix_A, x_len, x_len);
    destroyMatrix(matrix_L, x_len, x_len);
    destroyMatrix(matrix_U, x_len, x_len);
    destroyVector(vector_b);
    destroyVector(vector_u);
    destroyVector(vector_y);
    delete[] indexes;

    return max_error;
}



void LU_decomposition(double **matrixL, double **matrixU, int *indexes, int n) {
    double factor;

    for(int i = 0; i < n - 1; i++) {
        for(int j = i + 1; j < n; j++) {
            if(matrixU[indexes[i]][i] == 0) {
                if(!choice(matrixU, indexes, i, n)) {
                    break;
                }
            }
            factor = matrixU[indexes[j]][i] / matrixU[indexes[i]][i];
            matrixL[indexes[j]][i] = factor;
            for(int k = i; k < n; k++) {
                matrixU[indexes[j]][k] -= matrixU[indexes[i]][k] * factor;
            }
        }
    }

    for(int i = 0; i < n; i++)
        matrixL[indexes[i]][i] = 1;
}

double *L_matrixEquation(double **matrixL, double *vectorb, int *indexes, int n) {
    double *vectory = new double[n];
    double sum;

    for(int i = 0; i < n; i++) {
        sum = 0;
        for(int j = 0; j < i; j++)
            sum += vectory[j] * matrixL[indexes[i]][j];
        vectory[i] = (vectorb[indexes[i]] - sum) / matrixL[indexes[i]][i];
    }

    return vectory;
}

double *U_matrixEquation(double **matrixU, double *vectory, int *indexes, int n) {
    double *vectorx = new double[n];
    double sum;

    for(int i = n - 1; i >= 0; i--) {
        sum = 0;
        for(int j = i + 1; j < n; j++)
            sum += vectorx[j] * matrixU[indexes[i]][j];
        vectorx[i] = (vectory[i] - sum) / matrixU[indexes[i]][i];
    }

    return vectorx;
}

void matrixUfromA(double **matrixA, double **matrixU, int n) {
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            matrixU[i][j] = matrixA[i][j];
}

bool choice(double **matrixA, int *indexes, int x, int n) {
    double max = 0;
    int indexes_max = x, temp;
    bool find = false;

    for (int i = x; i < n; ++i) {
        if (fabs(matrixA[indexes[i]][x]) > max) {
            max = fabs(matrixA[indexes[i]][x]);
            indexes_max = i;
            find = true;
        }
    }

    if (find) {
        temp = indexes[indexes_max];
        indexes[indexes_max] = indexes[x];
        indexes[x] = temp;
    }

    return find;
}

void indexesInput(int *vector, int columns) {
    for(int i = 0; i < columns; i++)
        vector[i] = i;
}

void matrixZeros(double **matrix, int rows, int columns) {
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < columns; j++)
            matrix[i][j] = 0.0;
}

double **createDoubleMatrix(int rows, int columns) {
    double **matrix = new double* [rows];
    for(int i = 0; i < rows; i++)
        matrix[i] = new double[columns];

    return matrix;
}

void destroyMatrix(double **matrix, int rows, int columns) {
    for(int i = 0; i < rows; i++)
        delete [] matrix[i];
    delete [] matrix;
}

double *createDoubleVector(int columns) {
    double *vector = new double[columns];
    return vector;
}

void destroyVector(double *vector) {
    delete [] vector;
}

