//
// Created by volie on 16.03.2022.
//

#include <iostream>
#include <math.h>

using namespace std;

double pierwsza_funkcja(double x) {
    return pow(sin(x/4.0), 2.0) - x;
}

double pochodna_pierwsza_funkcja(double x) {
    return 1.0/2.0 * sin(x/4.0) * cos(x/4.0) - 1.0;
}

double druga_funkcja(double x) {
    return tan(2.0*x) - x - 1;
}

double pochodna_druga_funkcja(double x) {
    return 2.0 / pow(cos(2.0*x), 2.0) - 1.0;
}

void metoda_Picarda() {
    cout << "Metoda Picarda: " << endl;
}

void metoda_bisekcji() {
    cout << "Metoda Bisekcji: " << endl;
}

void metoda_Newtona() {
    cout << "Metoda Newtona: " << endl;
}

void metoda_siecznych() {
    cout << "Metoda Siecznych: " << endl;
}


int main() {
    cout << "Pierwsza funkcja: sin^2(x/4) - x = 0" << endl;
    metoda_Picarda();
    metoda_bisekcji();
    metoda_Newtona();
    metoda_siecznych();

    cout << "Druga funkcja: tan(2x) - x - 1 = 0" << endl;
    metoda_Picarda();
    metoda_bisekcji();
    metoda_Newtona();
    metoda_siecznych();

    return 0;
}