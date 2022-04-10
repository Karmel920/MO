//
// Created by volie on 08.04.2022.
//

#include <iostream>
#include <fstream>

using namespace std;

double *createDoubleVector(int columns);
void printVector(double *vector, int columns);
void vectorsInput(double *u, double *d, double *l, double*b, int columns);
void destroyVector(double *vector);
void Thomass(double *et, double *r, double *u, double *d, double *l, double *b, int n);
void vectorX(double *et, double *r, double *u, double *x, int n);

int main(){
    int n = 6;
    double *u = createDoubleVector(n - 1);
    double *d = createDoubleVector(n);
    double *l = createDoubleVector(n - 1);
    double *x = createDoubleVector(n);
    double *b = createDoubleVector(n);
    double *et = createDoubleVector(n);
    double *r = createDoubleVector(n);

    u[0] = 1.0 / 2.0; u[1] = 1.0 / 4.0; u[2] = 1.0 / 6.0; u[3] = 1.0 / 8.0; u[4] = 1.0 / 10.0;
    d[0] = 10.0; d[1] = 20.0; d[2] = 30.0; d[3] = 30.0; d[4] = 20.0; d[5] = 10.0;
    l[0] = 1.0 / 3.0; l[1] = 1.0 / 5.0; l[2] = 1.0 / 7.0; l[3] = 1.0 / 9.0; l[4] = 1.0 / 11.0;
    b[0] = 31.0; b[1] = 165.0 / 4.0; b[2] = 917.0 / 30.0; b[3] = 851.0 / 28.0; b[4] = 3637.0 / 90.0;
    b[5] = 332.0 / 11.0;

//    for(int i = 0; i < n; i++) {
//        x[i] = 1;
//        r[i] = 1;
//        et[i] = 1;
//    }

    cout << endl << "\n\tWektor u: " << endl;
    printVector(u, n-1);
    cout << endl << "\n\tWektor d: " << endl;
    printVector(d, n);
    cout << endl << "\n\tWektor l: " << endl;
    printVector(l, n-1);
    cout << endl << "\n\tWektor b: " << endl;
    printVector(b, n);
    cout << endl;

    Thomass(et, r, u, d, l, b, n);

    cout << endl << "\n\tWektor eta: " << endl;
    printVector(et, n);
    cout << endl << "\n\tWektor r: " << endl;
    printVector(r, n);
    cout << endl;

    vectorX(et, r, u, x, n);

    cout << endl << "\n\tWektor x: " << endl;
    printVector(x, n);
    cout << endl;

    destroyVector(u);
    destroyVector(d);
    destroyVector(l);
    destroyVector(b);
//    destroyVector(x);
//    destroyVector(et);
//    destroyVector(r);

    return 0;
}



double *createDoubleVector(int columns) {
    double *vector = new double[columns];
    return vector;
}

void vectorsInput(double *u, double *d, double *l, double*b, int columns) {
    ifstream myfile;
    myfile.open ("vectors.txt");
    for(int i = 0; i < columns - 1; i++)
        myfile >> u[i];
    for(int i = 0; i < columns; i++)
        myfile >> d[i];
    for(int i = 0; i < columns - 1; i++)
        myfile >> l[i];
    for(int i = 0; i < columns; i++)
        myfile >> b[i];
}

void printVector(double *vector, int columns) {
    for (int i = 0; i < columns; i++)
        cout << vector[i] << " ";
    cout << endl;
}

void destroyVector(double *vector) {
    delete [] vector;
}

void Thomass(double *et, double *r, double *u, double *d, double *l, double *b, int n) {
    et[0] = d[0];
    r[0] = b[0];
    for(int i = 1; i <= n; i++) {
        et[i] = d[i] - ((l[i - 1] * u[i - 1]) / et[i - 1]);
        r[i] = b[i] - ((l[i - 1] * r[i - 1]) / et[i - 1]);
    }
}

void vectorX(double* et, double *r, double *u, double *x, int n) {
    x[n] = r[n] / et[n];
    for(int i = n - 1; i >= 0; i--)
        x[i] = (r[i] - (u[i] * x[i + 1])) / et[i];
}