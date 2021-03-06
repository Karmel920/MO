#include <iostream>
#include <math.h>
#include <algorithm>

using namespace std;

//funkcje z ukladu
double funkcja1(double x, double y, double z) {
    return x*x + y*y + z*z - 2;
}

double funkcja2(double x, double y) {
    return x*x + y*y - 1;
}

double funkcja3(double x, double y) {
    return x*x - y;
}


//pochodne czastkowe pierwszej funkcji
double f1x(double x) {
    return 2*x;
}

double f1y(double y) {
    return 2*y;
}

double f1z(double z) {
    return 2*z;
}


//pochodne czastkowe drugiej funkcji
double f2x(double x) {
    return 2*x;
}

double f2y(double y) {
    return 2*y;
}

double f2z() {
    return 0;
}


//pochodne czastkowe trzeciej funkcji
double f3x(double x) {
    return 2*x;
}

double f3y() {
    return -1;
}

double f3z() {
    return 0;
}


//jakobian(wyznacznik macierzy jakobiego)
double detJ(double x, double y, double z) {
    return -4*x*z * (1 + 2*y);
}


int main() {
    const double TOLX = 1e-12;
    const double TOLF = 1e-12;
    int n_max = 50;
    double x = 3.0, y = 4.0, z = 5.0;
    double delta_x, delta_y, delta_z, x1, y1, z1;
    double est_x, est_y, est_z, est;
    double fn1, fn2, fn3, fn;

    cout << "   n \t  xn \t\t   yn\t\tzn\t EST\t\t\t |fn|" << endl;
//    |fn2|	 |fn3|

    for(int i = 1; i < n_max; i++) {
        delta_x = (funkcja2(x, y) * f3y() * f1z(z) - f1z(z) * f2y(y) * funkcja3(x, y)) / detJ(x, y, z);
        delta_y = (f2x(x) * funkcja3(x, y) * f1z(z) - f1z(z) * funkcja2(x, y) * f3x(x)) / detJ(x, y, z);
        delta_z = (f1x(x) * f2y(y) * funkcja3(x, y) + f2x(x) * f3y() * funkcja1(x, y, z) +
                   f3x(x) * f1y(y) * funkcja2(x, y) -
                   (funkcja1(x, y, z) * f2y(y) * f3x(x) + funkcja2(x, y) * f3y() * f1x(x) +
                    funkcja3(x, y) * f1y(y) * f2x(x))) / detJ(x, y, z);
        x1 = x - delta_x;
        y1 = y - delta_y;
        z1 = z - delta_z;

        est_x = fabs(x - x1);
        est_y = fabs(y - y1);
        est_z = fabs(z - z1);

        fn1 = fabs(funkcja1(x1, y1, z1));
        fn2 = fabs(funkcja2(x, y));
        fn3 = fabs(funkcja3(x, y));

        est = max(max(est_x, est_y),est_z);
        fn = max(max(fn1, fn2), fn3);

        cout.width(4);
        cout << i << "\t" << x1 << " \t" << y1 << " \t" << z1 << "    \t"<< est << "      \t" << fabs(fn)  << endl;

        if(est < TOLX && fn < TOLF) break;

        x = x1;
        y = y1;
        z = z1;
    }

    return 0;
}