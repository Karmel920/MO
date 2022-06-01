//
// Created by volie on 13.04.2022.
//

#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>

using namespace std;

const int n_max = 100;
const double EPS = 1e-13;

double **createDoubleMatrix(int rows, int columns);
double *createDoubleVector(int columns);
void destroyMatrix(double **matrix, int rows);
void destroyVector(double *vector);
void matrixInput(double **matrix, int rows, int columns);
void vectorInput(double *vector, int columns);
double norm(double *x, int n);
double *residuum(double **matrix, double *x, double *b, int n);
void printVector(int k, double *x, int n, double est, double res);
void Jacobi(double **matrix, double *b, double *x, int n);
void Gauss_Seidel(double **matrix, double *b, double *x, int n);
void SOR(double **matrix, double *b, double *x, int n, double omega);


int main() {
    int n = 4;
    double omega = 0.5;
    double **matrixA = createDoubleMatrix(n, n);
    double *vectorb = createDoubleVector(n);
    double *vectorx = createDoubleVector(n);

    matrixInput(matrixA, n, n);
    vectorInput(vectorb, n);

    vectorx[0] = 2.0; vectorx[1] = 2.0; vectorx[2] = 2.0; vectorx[3] = 2.0;
    cout << "Metoda Jacobiego: " << endl;
    Jacobi(matrixA, vectorb, vectorx, n);

    vectorx[0] = 2.0; vectorx[1] = 2.0; vectorx[2] = 2.0; vectorx[3] = 2.0;
    cout << "\nMetoda Gaussa-Seidela:" << endl;
    Gauss_Seidel(matrixA, vectorb, vectorx, n);

    vectorx[0] = 2.0; vectorx[1] = 2.0; vectorx[2] = 2.0; vectorx[3] = 2.0;
    cout << "\nMetoda SOR:" << endl;
    SOR(matrixA, vectorb, vectorx, n, omega);

    destroyMatrix(matrixA, n);
    destroyVector(vectorx);
    destroyVector(vectorb);

    return 0;
}


double **createDoubleMatrix(int rows, int columns) {
    double **matrix = new double* [rows];
    for(int i = 0; i < rows; i++)
        matrix[i] = new double[columns];

    return matrix;
}

void destroyMatrix(double **matrix, int rows) {
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


void matrixInput(double **matrix, int rows, int columns) {
    ifstream myfile;
    myfile.open ("matrix.txt");
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < columns; j++)
            myfile >> matrix[i][j];
}

void vectorInput(double *vector, int columns) {
    ifstream myfile;
    myfile.open ("vector.txt");
    for(int i = 0; i < columns; i++)
        myfile >> vector[i];
}

double norm(double *x, int n) {
    double max = fabs(x[0]);
    for(int i = 1; i < n; i++){
        if(fabs(x[i]) > max)
            max = fabs(x[i]);
    }
    return max;
}

double *residuum(double **matrix, double *x, double *b, int n) {
    double sum = 0;
    double *c = createDoubleVector(n);
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            sum += matrix[i][j]*x[j];
        }
        c[i] = sum - b[i];
        sum = 0;
    }
    return c;
}

void printVector(int k, double *x, int n, double est, double res) {
    cout.width(8);
    cout << left << k << "\t";
    for(int i = 0; i < n; i++) {
        cout.width(8);
        cout << left << x[i] << "\t";
    }
    cout.width(8);
    cout << left << est << "\t" << left << res << "\t" << endl;
}

void Jacobi(double **matrix, double *b, double *x, int n) {
    double *u = createDoubleVector(n);
    double sum = 0, est = 0, res = 0;
    int k = 0;

    cout.width(8);
    cout << "k \t\t x1 \t\t x2 \t\tx3 \t\t x4 \t\tEST \t\t RES" << endl;
    printVector(k, x, n, est, res);

    for(k = 1; k <= n_max; k++) {
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                if(j != i) {
                    sum += matrix[i][j] * x[j];
                }
            }
            u[i] = (b[i] - sum) / matrix[i][i];
            sum = 0;
        }

        est = fabs(norm(u, n) - norm(x, n));
        res = fabs(norm(residuum(matrix, u, b, n), n));
        printVector(k, u, n, est, res);
        if(est < EPS && res < EPS) break;

        for(int i = 0; i < n; i++) {
            x[i] = u[i];
        }
    }
    destroyVector(u);
}

void Gauss_Seidel(double **matrix, double *b, double *x, int n) {
    double *temp = createDoubleVector(n);
    double sum = 0.0, est = 0, res = 0;
    int k = 0;

    cout.width(8);
    cout << "k \t\t x1 \t\t x2 \t\tx3 \t\t x4 \t\tEST \t\t RES" << endl;
    printVector(k, x, n, est, res);

    for(k = 1; k <= n_max; k++) {
        for(int i = 0; i < n; i++) {
            temp[i] = x[i];
            for(int j = 0; j < n; j++) {
                if(j != i) {
                    sum += matrix[i][j] * x[j];
                }
            }
            x[i] = (b[i] - sum) / matrix[i][i];
            sum = 0;
        }

        est = fabs(norm(x, n) - norm(temp, n));
        res = fabs(norm(residuum(matrix, x, b, n), n));
        printVector(k, x, n, est, res);
        if(est < EPS && res < EPS) break;
    }
    destroyVector(temp);
}

void SOR(double **matrix, double *b, double *x, int n, double omega) {
    double *temp = createDoubleVector(n);
    double sum = 0.0, est = 0, res = 0;
    int k = 0;

    cout.width(8);
    cout << "k \t\t x1 \t\t x2 \t\tx3 \t\t x4 \t\tEST \t\t RES" << endl;
    printVector(k, x, n, est, res);

    for(k = 1; k <= n_max; k++) {
        for(int i = 0; i < n; i++) {
            temp[i] = x[i];
            for(int j = 0; j < n; j++) {
                if(j == i)
                    sum += (1.0 - (1.0/omega)) * matrix[i][i] * x[j];
                else
                    sum += matrix[i][j] * x[j];
            }
            x[i] = (b[i] - sum) / (matrix[i][i] / omega);                  //rozwiazanie ukladu z dolna macierza
            sum = 0;
        }

        est = fabs(norm(x, n) - norm(temp, n));
        res = fabs(norm(residuum(matrix, x, b, n), n));
        printVector(k, x, n, est, res);
        if(est < EPS && res < EPS) break;
    }
    destroyVector(temp);
}



