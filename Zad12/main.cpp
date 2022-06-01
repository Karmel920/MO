//
// Created by volie on 31.05.2022.
//

#include <iostream>
#include <fstream>
#include <cmath>

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#ifndef M_PIl
#define M_PIl (3.14159265358979323846264338327950288)
#endif

using namespace std;

double exact_value(double x);

void equidistant_nodes(double *x);

void czebyszew_nodes(double *x);

void newton(string name, double *x);

void exact_values_file(string name);

const double a = -1.0;
const double b = 1.0;
const double h = 0.01;
const int N = 13;


int main() {
    double equidistant_x[N];
    double czebyszew_x[N + 1];

    equidistant_nodes(equidistant_x);
    newton("equidistant.txt", equidistant_x);

    czebyszew_nodes(czebyszew_x);
    newton("czebyszew.txt", czebyszew_x);

    exact_values_file("exact.txt");

    system("gnuplot -p -e \"set ylabel 'y'; set xlabel 'x'; plot 'exact.txt' with lines, 'czebyszew.txt' , 'equidistant.txt'\"");
    return 0;
}

double exact_value(double x) {
    return 1.0 / (1.0 + 10.0 * pow(x, 6.0));
}

void equidistant_nodes(double *x) {
    double H = ((b - a) / (N - 1.0));
    for (int i = 0; i < N; i++) {
        x[i] = a + i * H;
    }
}

void czebyszew_nodes(double *x) {
    for (int i = 0; i <= N; i++) {
        x[i] = (b + a) / 2.0 + ((b - a) / 2.0) * cos(((2.0 * i + 1.0) / (2.0 * N + 2.0)) * M_PI);
    }
}

void newton(string name, double *x) {
    double c[N][N];
    double result;
    fstream file;
    file.open(name, fstream::out);
    file << scientific;

    //tworzenie bazy newtona
    for (int i = 0; i < N; ++i) {
        c[i][0] = exact_value(x[i]);
    }
    for (int i = 1; i < N; ++i) {
        for (int j = 0; j < N - i; ++j) {
            c[j][i] = (c[j + 1][i - 1] - c[j][i - 1]) / (x[i + j] - x[j]);
        }
    }

    //algrytm Hornera
    for (double xi = a; xi <= b; xi += h) {
        result = c[0][N - 1];
        for (int j = N - 1; j > 0; --j) {
            result *= (xi - x[j - 1]);
            result += c[0][j - 1];
        }
        file << xi << "\t" << result << endl;
    }
    file.close();
}

void exact_values_file(string name) {
    fstream file;
    file.open(name, fstream::out);
    file << scientific;
    for (double x = a; x <= b; x += h) {
        file << x << "\t" << exact_value(x) << endl;
    }
    file.close();
}