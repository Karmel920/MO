//
// Created by volie on 16.03.2022.
//

#include <iostream>
#include <math.h>

using namespace std;

const int n_max = 50;
const double TOLX = 1e-12;
const double TOLF = 1e-12;


double pierwsza_funkcja(double x) {
    return pow(sin(x/4.0), 2.0) - x;
}

double pochodna_pierwsza_funkcja(double x) {
    return 0.5 * sin(x/4.0) * cos(x/4.0) - 1.0;
}

double pierwsza_funckja_picard(double x) {
    return pow(sin(x/4.0), 2.0);
}

double pierwsza_funkcja_picard_pochodna(double x) {
    return 0.5 * sin(x/4.0) * cos(x/4.0);
}


double druga_funkcja(double x) {
    return tan(2.0*x) - x - 1;
}

double pochodna_druga_funkcja(double x) {
    return 2.0 / pow(cos(2.0*x), 2.0) - 1.0;
}

double druga_funkcja_picard(double x) {
    return tan(2.0*x) - 1;
}
double druga_funkcja_picard_pochodna(double x) {
    return 2.0 / pow(cos(2.0*x), 2.0);
}


void metoda_Picarda(double(*funkcja) (double), double(*funkcja_picard)(double), double(*funkcja_pochodna)(double), double xn) {
    cout << "Metoda Picarda: " << endl;

    if(fabs(funkcja_pochodna(xn)) >= 1) {
        cout << "Fukcja rozbiezna" << endl;
    }
    else {
        double xn1, fxn, EST;
        cout << "   n \t\t   xn \t\t    EST \t \t  |f(xn)|" << endl;
        for (int i = 0; i < n_max; i++) {
            fxn = fabs(funkcja(xn));
            xn1 = funkcja_picard(xn);
            EST = fabs(xn1 - xn);
            xn = xn1;

            cout.width(4);
            cout << i << "\t\t" << xn << "\t\t" << EST << " \t\t" << fxn <<endl;

            if (EST < TOLX && fxn < TOLF) break;
        }
    }
}

void metoda_bisekcji(double(*funkcja)(double), double a, double b) {
    cout << "Metoda Bisekcji: " << endl;

    if (funkcja(a) * funkcja(b) > 0) {
        cout << "Zly przedzial." << endl;
    }
    else {
        double EST, xn, fxn;
        cout << "   n \t\t    xn \t\t       EST \t \t |f(xn)|" << endl;
        for (int i = 0; i < n_max; i++) {
            xn = (a + b) / 2;
            EST = fabs(b - a) / 2;
            fxn = funkcja(xn);
            if (funkcja(a) * fxn < 0)
                b = xn;
            else
                a = xn;

            cout.width(4);
            cout << i << "\t\t" << xn << "\t\t" << EST << " \t\t" << fabs(fxn) << endl;

            if (fabs(fxn) < TOLF) break;
//            (EST < TOLX && fabs(fxn) < TOLF) ||
        }
    }
}

void metoda_Newtona(double(*funkcja)(double), double(*funkcja_pochodna)(double), double xn) {
    cout << "Metoda Newtona: " << endl;

    double fxn, xn1, EST;
    cout << "   n \t\t    xn \t\t       EST \t \t |f(xn)|" << endl;
    for (int i = 0; i < n_max; i++) {
        xn1 = xn;
        fxn = funkcja(xn1);
        xn = xn1 - fxn / funkcja_pochodna(xn1);
        EST = fabs(xn1 - xn);

        cout.width(4);
        cout << i << "\t\t" << xn << "\t\t" << EST << " \t\t" << fabs(fxn) << endl;

        if (EST < TOLX && fabs(fxn) < TOLF) break;
    }
}

void metoda_siecznych(double(*funkcja)(double), double xn0, double xn1) {
    cout << "Metoda Siecznych: " << endl;

    double xn2, EST, fxn2;
    cout << "   n \t\t    xn \t\t       EST \t \t |f(xn)|" << endl;
    for (int i = 0; i < n_max; i++) {
        xn2 = xn1 - funkcja(xn1) / ((funkcja(xn1) - funkcja(xn0)) / (xn1 - xn0));
        xn0 = xn1;
        xn1 = xn2;
        fxn2 = fabs(funkcja(xn1));
        EST = fabs(xn1 - xn0);

        cout.width(4);
        cout << i << "\t\t" << xn2 << "\t\t" << EST << " \t\t" << fabs(fxn2) << endl;

        if (EST < TOLX && fxn2 < TOLF) break;
    }
}


int main() {
    cout << "Pierwsza funkcja: sin^2(x/4) - x = 0" << endl;
    metoda_Picarda(pierwsza_funkcja, pierwsza_funckja_picard, pierwsza_funkcja_picard_pochodna, 0.5);
    metoda_bisekcji(pierwsza_funkcja, -4.0, 6.2);
    metoda_Newtona(pierwsza_funkcja, pochodna_pierwsza_funkcja, 0.6);
    metoda_siecznych(pierwsza_funkcja, 5.0, 10.0);

    cout << endl << "Druga funkcja: tan(2x) - x - 1 = 0" << endl;
    metoda_Picarda(druga_funkcja, druga_funkcja_picard, druga_funkcja_picard_pochodna, 5.0);
    metoda_bisekcji(druga_funkcja, -0.1, 0.7);
    metoda_Newtona(druga_funkcja, pochodna_druga_funkcja, 5.0);
    metoda_siecznych(druga_funkcja, 10.0, 5.0);

    return 0;
}