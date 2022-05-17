#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

const double alfa = 0.0, beta = 1.0, gamma = -1.0;
const double phi = 0.0, psi = 1.0, theta = 0.0;
const double a = 0, b = 1;

const double p = 1.0, q = 0.0, r = -4.0;

fstream plik_bledy, plik_konwencjonalna, plik_numerow;

double s(double x) {
    return -x;
}

double rozw_analityczne(double x) {
    return (exp(2.0-2.0*x)-4.0*exp(4.0-2.0*x)+4.0*exp(2*x)-exp(2.0+2.0*x)-x+x*exp(4.0)) / (4.0-4.0*exp(4.0));
}

void Thomas(double* L, double* D, double*U, double* b, double* x, int n) {
    double *temp_b = new double[n];
    double *temp_d = new double[n];

    temp_b[0] = b[0];
    temp_d[0] = D[0];

    for (int i = 1; i < n; i++) {
        temp_d[i] = D[i] - L[i - 1] * (U[i - 1] / temp_d[i - 1]);
        temp_b[i] = b[i] - L[i - 1] * temp_b[i - 1] / temp_d[i - 1];
    }

    x[n - 1] = temp_b[n - 1] / temp_d[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        x[i] = (temp_b[i] - U[i] * x[i + 1]) / temp_d[i];
    }
    delete[] temp_d;
    delete[] temp_b;
}

double oblicz_blad(double* bledy, int n) {
    double max = fabs(bledy[0]);
    for (int i = 0; i < n; i++) {
        if (fabs(bledy[i]) > max) {
            max = bledy[i];
        }
    }
    return max;
}

double dyskretyzacja_numerowa(double h, int n) {
    double *L = new double[n];
    double *D = new double[n];
    double *U = new double[n];
    double *b = new double[n];
    double *x = new double[n];
    double *bledy = new double[n];

    double x1 = a;
    double x2 = a;
    D[0] = beta - alfa / h;
    U[0] = alfa / h;
    b[0] = -gamma;

    for (int i = 1; i < n - 1; i++) {
        L[i - 1] = p / (h * h) + r / 12.0;
        D[i] = (-2.0 * p) / (h * h) + r * (10.0 / 12.0);
        U[i] = p / (h * h) + r / 12.0;
        b[i] = -(s(x1 + i * h - h) + 10.0 * s(x1 + i * h) + s(x1 + i * h + h)) / 12.0;
    }

    L[n - 2] = -phi / h;
    D[n - 1] = -phi / h + psi;
    b[n - 1] = -theta;
    Thomas(L, D, U, b, x, n);

    if (n == 19900) {
        for (int i = 0; i < n; i++) {
            plik_numerow << x2 << "\t" << rozw_analityczne(x2) << "\t" << x[i] << endl;
            x2 += h;
        }
    }

    for (int i = 0; i < n; i++) {
        bledy[i] = fabs(x[i] - rozw_analityczne(x1));
        x1 += h;
    }

    double max_blad = oblicz_blad(bledy, n);

    delete[] L;
    delete[] U;
    delete[] D;
    delete[] x;
    delete[] b;
    delete[] bledy;

    return max_blad;
}

double dyskretyzacja_konwencjonalna(double h, int n) {
    double *L = new double[n];
    double *D = new double[n];
    double *U = new double[n];
    double *b = new double[n];
    double *x = new double[n];
    double *bledy = new double[n];

    double x1 = a;
    double x2 = a;
    D[0] = beta - alfa / h;
    U[0] = alfa / h;
    b[0] = -gamma;

    for (int i = 1; i < n - 1; i++) {
        L[i - 1] = p / (h * h) - q / (2.0 * h);
        D[i] = r + (-2.0 * p) / (h * h);
        U[i] = p / (h * h) + q / (2.0 * h);
        b[i] = -s((x1 + i * h));
    }

    L[n - 2] = -phi / h;
    D[n - 1] = -phi / h + psi;
    b[n - 1] = -theta;

    Thomas(L, D, U, b, x, n);

    if (n == 19900) {
        for (int i = 0; i < n; i++) {
            plik_konwencjonalna << x2 << "\t" << rozw_analityczne(x2) << "\t" << x[i] << endl;
            x2 += h;
        }
    }

    for (int i = 0; i < n; i++) {
        bledy[i] = fabs(x[i] - rozw_analityczne(x1));
        x1 += h;
    }

    double max_blad = oblicz_blad(bledy, n);
    delete[] L;
    delete[] U;
    delete[] D;
    delete[] x;
    delete[] b;
    delete[] bledy;

    return max_blad;
}

int main() {
    double h;
    double blad_konwencjonalna, blad_numerow;
    plik_bledy.open("bledy.txt", ios::out);
    plik_numerow.open("numerow.txt", ios::out);
    plik_konwencjonalna.open("konwencjonalna.txt", ios::out);

    if (!plik_bledy.good() || !plik_numerow.good() || !plik_konwencjonalna.good()) {
        cout << "Problem z otwarciem plikow!" << endl;
        exit(0);
    }

    plik_numerow << "x " << "analityczne " << "numerow" << endl;
    plik_konwencjonalna << "x " << "analityczne " << "konwencjonalna" << endl;
    plik_bledy << "log10(h) " << "KonwencjonalnaError " << "NumerowError" << endl;

    plik_bledy << scientific;
    plik_numerow << scientific;
    plik_konwencjonalna << scientific;

    vector<double> blad_konw;
    vector<double> blad_num;
    vector<double> logh;

    for (int i = 10; i < 20000; i += 10) {
        h = (b - a) / (i - 1);
        blad_konwencjonalna = dyskretyzacja_konwencjonalna(h, i);
        blad_numerow = dyskretyzacja_numerowa(h, i);
        plik_bledy << log10(h) << "\t" << log10(blad_konwencjonalna) << "\t" << log10(blad_numerow) << endl;
        logh.push_back(log10(h));
        blad_num.push_back(log10(blad_numerow));
        blad_konw.push_back(log10(blad_konwencjonalna));
    }

    double rzad_konwencjonalna, rzad_numerow;
    rzad_konwencjonalna = (blad_konw[10] - blad_konw[9]) / (logh[10] - logh[9]);
    rzad_numerow = (blad_num[10] - blad_num[9]) / (logh[10] - logh[9]);
    cout << "Rzad dokladnosci konwencjonalna: " << rzad_konwencjonalna << endl;
    cout << "Rzad dokladnosci Numerow: " << rzad_numerow << endl;

    system("gnuplot -p -e \" plot for [col=2:3] 'konwencjonalna.txt' using 1:col title columnheader w l\"");
    system("gnuplot -p -e \" plot for [col=2:3] 'numerow.txt' using 1:col title columnheader w l\"");
    system("gnuplot -p -e \"set ylabel 'log10(maksymalny blad bezwzgledny)'; set xlabel 'log10(h)';  plot for [col=2:3] 'bledy.txt' using 1:col title columnheader w l\"");

    plik_konwencjonalna.close();
    plik_numerow.close();
    plik_bledy.close();

    return 0;
}