#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

const int N = 500;

double rozw_analityczne(double t) {
    return 1.0 - exp(-10.0 * (t + atan(t)));
}

double BME(double h, double t_max) {
    double y = 0.0;
    for (double t = 0.0; t < t_max; t += h) {
        y += h * (-((10.0 * t * t + 20.0) / (t * t + 1.0) * (y - 1.0)));
    }
    return y;
}

double PME(double h, double t_max) {
    double y = 0.0;
    double ulamek;
    for (double t = 0.0; t < t_max; t += h) {
        ulamek = (10.0 * (t + h) * (t + h) + 20.0) / ((t + h) * (t + h) + 1.0);
        y = (y + h * ulamek) / (1.0 + h * ulamek);
    }
    return y;
}

double PMT(double h, double t_max) {
    double y = 0.0;
    double ulamek1, ulamek2;
    for (double t = 0.0; t < t_max; t += h) {
        ulamek1 = (10.0 * t * t + 20.0) / (t * t + 1.0);
        ulamek2 = (10.0 * (t + h) * (t + h) + 20.0) / ((t + h) * (t + h) + 1.0);
        y = (y - (h * (ulamek1 * (y - 1.0) - ulamek2)) * 0.5) / (1.0 + ((h * ulamek2) * 0.5));
    }
    return y;
}

double oblicz_blad(double(*f)(double, double), double h) {
    double t = 0.0;
    double dokladna, przyblizona;
    double error = fabs(rozw_analityczne(t) - f(h, t));

    for (int i = 0; i < N; i++) {
        dokladna = rozw_analityczne(t);
        przyblizona = f(h, t);
        dokladna = fabs(dokladna - przyblizona);
        if (dokladna > error) {
            error = dokladna;
        }
        t += h;
    }
    return error;
}

int main() {

    double h_stabilne = 0.025;
    double h_niestabilne = 0.25;
    fstream bme_stabilne, bme_niestabilne, pme, pmt, bledy;

    bledy.open("bledy.txt", ios::out);
    bme_stabilne.open("bezposredni_stabilne.txt", ios::out);
    bme_niestabilne.open("bezposredni_niestabilne.txt", ios::out);
    pme.open("posredni.txt", ios::out);
    pmt.open("trapezow.txt", ios::out);

    bledy << scientific;
    bme_stabilne << scientific;
    bme_niestabilne << scientific;
    pme << scientific;
    pmt << scientific;
    bme_stabilne << "t " << "analityczne " << "BME" << endl;
    bme_niestabilne << "t " << "analityczne " << "BME" << endl;
    pme << "t " << "analityczne " << "PME" << endl;
    pmt << "t " << "analityczne " << "PMT" << endl;
    bledy << "log10(dt) " << "log10(bladBME) " << "log10(bladPME) " << "log10(bladPMT)" << endl;

    if (!bme_stabilne.good() || !pme.good() || !pmt.good() || !bme_niestabilne.good() || !bledy.good()) {
        cout << "Problem z otwarciem plikow" << endl;
    }

    for (double t = 0.0; t < 5.0; t += h_stabilne) {
        bme_stabilne << t << " " << rozw_analityczne(t) << " " << BME(h_stabilne, t) << endl;
        pme << t << " " << rozw_analityczne(t) << " " << PME(h_stabilne, t) << endl;
        pmt << t << " " << rozw_analityczne(t) << " " << PMT(h_stabilne, t) << endl;
    }

    for (double t = 0.0; t < 5.0; t += h_niestabilne) {
        bme_niestabilne << t << " " << rozw_analityczne(t) << " " << BME(h_niestabilne, t) << endl;
    }

    vector<double> log_dt;
    vector<double> bme_blad;
    vector<double> pme_blad;
    vector<double> pmt_blad;

    for (double i = 0.1; i > 1e-20; i /= 1.2) {
        bledy << log10(i) << " " << log10(oblicz_blad(BME, i)) << " " << log10(oblicz_blad(PME, i)) << " "
              << log10(oblicz_blad(PMT, i)) << endl;
        log_dt.push_back(log10(i));
        bme_blad.push_back(log10(oblicz_blad(BME, i)));
        pme_blad.push_back(log10(oblicz_blad(PME, i)));
        pmt_blad.push_back(log10(oblicz_blad(PMT, i)));
    }

    double rzad_bme = (bme_blad[31] - bme_blad[30]) / (log_dt[31] - log_dt[30]);
    double rzad_pme = (pme_blad[31] - pme_blad[30]) / (log_dt[31] - log_dt[30]);
    double rzad_pmt = (pmt_blad[31] - pmt_blad[30]) / (log_dt[31] - log_dt[30]);
    system("gnuplot -p -e \"set ylabel 'y'; set xlabel 't'; plot 'bezposredni_stabilne.txt' u 1:2 title columnheader w l, '' u 1:3 title columnheader\"");
    system("gnuplot -p -e \"set ylabel 'y'; set xlabel 't'; plot 'bezposredni_niestabilne.txt' u 1:2 title columnheader w l, '' u 1:3 title columnheader\"");
    system("gnuplot -p -e \"set ylabel 'y'; set xlabel 't'; plot 'posredni.txt' u 1:2 title columnheader w l, '' u 1:3 title columnheader\"");
    system("gnuplot -p -e \"set ylabel 'y'; set xlabel 't'; plot 'trapezow.txt' u 1:2 title columnheader w l, '' u 1:3 title columnheader\"");
    system("gnuplot -p -e \"set ylabel 'log10(blad)'; set xlabel 'log10(dt)';plot for [col=2:4] 'bledy.txt' using 1:col title columnheader w l\"");

    cout << "Rzad dokladnosci BME: " << rzad_bme << endl;
    cout << "Rzad dokladnosci PME: " << rzad_pme << endl;
    cout << "Rzad dokladnosci PMT: " << rzad_pmt << endl;

    return 0;
}