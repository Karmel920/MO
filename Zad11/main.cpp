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

void Thomass(double *u, double *d, double *l, double *b, double *x, int n);
void get_exact(double start_x, double max_x, double h, double start_t, double max_t, double delta_t);
double exact_value(double x, double t);
void ML_Thomass(double start_x, double max_x, double h, double start_t, double max_t, double delta_t);
void ML_LU();
void MCN_Thomass();
void MCN_LU();


int main() {
    double start_x = -A, start_t = 0.0, max_x = A, max_t = T_MAX, h = 0.0625, delta_t = h * h;

    get_exact(start_x, max_x, h, start_t, max_t, delta_t);

//    ML_Thomas(start_x, max_x, delta_x, start_t, max_t, delta_t);
//    ML_LU(start_x, max_x, delta_x, start_t, max_t, delta_t);
//    MCN_Thomass();
//    MCN_LU();

//    get_errors_t_max();
//    plot();
//    get_errors_t();

    return 0;
}



void get_exact(double start_x, double max_x, double h, double start_t, double max_t, double delta_t) {
    ofstream exact;
    exact.open("exact.txt");
    double iter_x = start_x;
    while(start_t <= max_t) {
        while(iter_x <= max_x) {
            exact << start_t << " " << iter_x << " " << exact_value(iter_x, start_t) << endl;
            iter_x += h;
        }
        iter_x = start_x;
        start_t += delta_t;
    }
}

double exact_value(double x, double t) {
    return (1.0 / (2.0 * sqrt(M_PIl * D * (TAU + t)))) * exp(-(x*x) / (4 * D * (TAU + t)));
}

double initial_condition(double x) {
    return (1.0 / (2.0 * sqrt(M_PIl * D * TAU)) * exp(-(x*x) / (4 * D * TAU));
}

void Thomass(double *u, double *d, double *l, double *b, double *x, int n) {
    for(int i = 1; i < n; i++) {
        d[i] = d[i] - ((l[i] * u[i - 1]) / d[i - 1]);
        b[i] = b[i] - ((l[i] * b[i - 1]) / d[i - 1]);
    }

    x[n - 1] = b[n - 1] / d[n - 1];
    for(int i = n - 2; i >= 0; i--)
        x[i] = (b[i] - (u[i] * x[i + 1])) / d[i];
}

void ML_Thomass(double start_x, double max_x, double h, double start_t, double max_t, double delta_t) {
    ofstream ML;
    ML.open("ML_plot_Thomas.txt");
    double lambda = delta_t / (h * h);
    int t_len = int((max_t - start_t) / delta_t) + 1;
    int x_len = int((max_x - start_x) / h) + 1;

    double **matrix_A = new double *[3];

    for (int i = 0; i < 3; i++) {
        matrix_A[i] = new double[x_len];
    }
    double *vector_b = new double[x_len];
    double *vector_u = new double[x_len];

    for (int i = 0; i < x_len; i++) {
        vector_b[i] = initial_condition(start_x + h * i);
    }

    vector_b[0] = 0.0;
    vector_b[x_len - 1] = 0.0;

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

    double *wektor_eta = new double[x_len];
    double *rozwiazanie = new double[x_len];
    thomas_eta(macierz_A, wektor_eta);
    ofstream thomas_05;
    thomas_05.open("thomas_05.txt");
    ofstream thomas_1;
    thomas_1.open("thomas_1.txt");
    ofstream thomas_15;
    thomas_15.open("thomas_15.txt");
    ofstream thomas_blad;
    thomas_blad.open("thomas_blad.txt");
    double _t = 0.0;
    for (int j = 0; j < przedzialy_t - 1; j++) {

        double *wektor_r = new double[wiersze];
        thomas_r(wektor_b, wektor_r, macierz_A[0], wektor_eta);
        rozwiaz(wektor_r, rozwiazanie, macierz_A[2], wektor_eta);

        for (int i = 0; i < przedzialy_x; i++) {
            wektor_b[i] = -rozwiazanie[i];
            if (j == 50)
                thomas_05 << -a + h * i << " " << wektor_b[i] << " " << endl;
            if (j == 100)
                thomas_1 << -a + h * i << " " << wektor_b[i] << " " << endl;
            if (j == 150)
                thomas_15 << -a + h * i << " " << wektor_b[i] << " " << endl;
        }

        wektor_b[0] = 0.;
        wektor_b[wiersze - 1] = 0.;
        delete[] wektor_r;

        if (j == przedzialy_t - 1) {
            wektor_rozw = wektor_b;
        }

        blad_t = blad_bezwzgledny_t(wektor_b, dt, h);
        thomas_blad << _t << " " << blad_t << endl;
        _t += dt;
    }
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

