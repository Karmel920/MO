//
// Created by volie on 30.03.2022.
//

#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

double **createDoubleMatrix(int rows, int columns);
double *createDoubleVector(int columns);
int *createIntVector(int columns);
void destroyMatrix(double **matrix, int rows, int columns);
void destroyVector(double *vector);
void printMatrix(double **matrix, int rows, int columns);
void printMatrixIndexes(double **matrix, int* indexes, int n);
void printVector(double *vector, int columns);
void matrixInput(double **matrix, int rows, int columns);
void matrixZeros(double **matrix, int rows, int columns);
void vectorInput(double *vector, int columns);
void matrixUfromA(double **matrixA, double **matrixU, int n);
void indexesInput(int *vector, int columns);
void LU_decomposition(double **matrixL, double **matrixU, int *indexes, int n);
double *L_matrixEquation(double **matrixL, double *vectorb, int *indexes, int n);
double *U_matrixEquation(double **matrixU, double *vectory, int *indexes, int n);
bool choice(double **matrixA, int *indexes, int x, int n);

int main() {
    int n = 4;
    double **matrixA = createDoubleMatrix(n, n);
    double **matrixL = createDoubleMatrix(n, n);
    double **matrixU = createDoubleMatrix(n, n);
    double *vectorb = createDoubleVector(n);
    double *vectory = createDoubleVector(n);
    double *vectorx = createDoubleVector(n);
    int *indexes = createIntVector(n);

    matrixInput(matrixA, n, n);
    vectorInput(vectorb, n);
    matrixZeros(matrixL, n, n);
    indexesInput(indexes, n);
    matrixUfromA(matrixA, matrixU, n);

    cout << "\n\tUklad Ax = b\n" << endl;
    cout << "\tMacierz A:" <<endl;
    printMatrix(matrixA, n, n);
    cout << "\n\tWektor b:" <<endl;
    printVector(vectorb, n);

    LU_decomposition(matrixL, matrixU, indexes, n);
    cout << "\n\tMacierz L:" <<endl;
    printMatrixIndexes(matrixL, indexes, n);
    cout << "\n\tMacierz U:" <<endl;
    printMatrixIndexes(matrixU, indexes, n);

    vectory = L_matrixEquation(matrixL, vectorb, indexes, n);
    vectorx = U_matrixEquation(matrixU, vectory, indexes, n);
    cout << "\n\tWektor y:" <<endl;
    printVector(vectory, n);
    cout << "\n\tWektor x:" <<endl;
    printVector(vectorx, n);

    destroyMatrix(matrixA, n, n);
    destroyMatrix(matrixL, n, n);
    destroyMatrix(matrixU, n, n);
    destroyVector(vectorb);
    destroyVector(vectory);
    destroyVector(vectorx);

    return 0;
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

int *createIntVector(int columns) {
    int *vector = new int[columns];
    return vector;
}

void destroyVector(double *vector) {
    delete [] vector;
}

void printMatrix(double **matrix, int rows, int columns) {
    for(int i = 0; i < rows; i++) {
        for(int j = 0; j < columns; j++) {
            cout << matrix[i][j] << "\t";
        }
        cout << endl;
    }
}

void printVector(double *vector, int columns) {
    for (int i = 0; i < columns; i++)
        cout << vector[i] << "\t";
    cout << endl;
}

void matrixInput(double **matrix, int rows, int columns) {
    ifstream myfile;
    myfile.open ("matrix.txt");
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < columns; j++)
            myfile >> matrix[i][j];
}

void matrixZeros(double **matrix, int rows, int columns) {
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < columns; j++)
            matrix[i][j] = 0;
}

void vectorInput(double *vector, int columns) {
    ifstream myfile;
    myfile.open ("vector.txt");
    for(int i = 0; i < columns; i++)
        myfile >> vector[i];
}

void indexesInput(int *vector, int columns) {
    for(int i = 0; i < columns; i++)
        vector[i] = i;
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

void printMatrixIndexes(double **matrix, int *indexes, int n) {
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            cout << matrix[indexes[i]][j] << "\t";
        }
        cout << endl;
    }
}




