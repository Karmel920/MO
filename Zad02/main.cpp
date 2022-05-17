//
// Created by volie on 13.03.2022.
//

#include <iostream>
#include <math.h>
#include <fstream>
#include <limits>
#include <vector>

using namespace std;

double taylor(double x) {
    double w = 1.0;
    double suma= 1.0;

    for (int i = 1; i < 40; i++) {
        w *= -x / (i + 1.0);
        suma += w;
    }
    return suma;
}

//double taylor(double x)
//{
//    double suma = 1.0;
//    double poprzedni = 1.0;
//    for(int i=2; i<40; i++)
//    {
//        poprzedni*=x/(double)i;
//        poprzedni = -poprzedni;
//        suma+=poprzedni;
//    }
//    return suma;
//}

double wartosc_przyblizona(double x) {
    double f = (1 - exp(-x)) / x;
    return f;
}

double blad_wzgledny(double w, double wp) {
    return fabs((wp - w) / w);
}

double log_blad_wzgledny(double blad) {
    return log10(blad);
}

int main() {
    vector <double> log_10, x, fx, log_10_blad, log_10_git;
    double e1, e2, e3;
    ifstream input("dane.txt");
    while(input >> e1 >> e2 >> e3) {
        log_10.push_back(e1);
        x.push_back(e2);
        fx.push_back(e3);
    }

    double wartosc_przyb, blad_wzg, log_blad_wzg, log_git, blad_wzg_tay;

    for(int i = 0; i < x.size(); i++) {
        wartosc_przyb = wartosc_przyblizona(x[i]);
        blad_wzg = blad_wzgledny(fx[i], wartosc_przyb);

        if(blad_wzg > numeric_limits<double>::epsilon()) {
            wartosc_przyb = taylor(x[i]);
            blad_wzg_tay = blad_wzgledny(fx[i], wartosc_przyb);
        }

        if(blad_wzg > blad_wzg_tay)
            blad_wzg = blad_wzg_tay;

        log_blad_wzg = log_blad_wzgledny(blad_wzg);
        log_10_blad.push_back(log_blad_wzg);

        log_git = log_blad_wzgledny(fx[i]);
        log_10_git.push_back(log_git);
    }

    fstream output;
    output.open("log10(blad)_od_log10(x)taylor.txt", ios::out);
    for(int i = 0; i < log_10.size(); i++)
        if(output.is_open()) {
            output << scientific << log_10[i] << "\t\t\t";
            output << scientific << log_10_blad[i] << endl;
        }
    output.close();

//    fstream output2;
//    output2.open("lot.txt", ios::out);
//    for(int i = 0; i < log_10.size(); i++)
//        if(output2.is_open()) {
//            output2 << scientific << log_10[i] << "\t\t\t";
//            output2 << scientific << log_10_git[i] << endl;
//        }
//    output2.close();

//    fstream output2;
//    output2.open("log10(bladwzgl).txt", ios::out);
//    for(int i = 0; i < log_10_blad.size(); i++)
//        if(output2.is_open())
//            output2 << scientific << log_10_blad[i] << endl;
//    output2.close();

//    for(int i = 0; i < log_10.size(); i++)
//        cout << scientific << log_10[i] << " " << scientific <<  x[i] << " " << scientific << fx[i] << endl;

    log_10.clear();
    x.clear();
    fx.clear();
    log_10_blad.clear();
    log_10_git.clear();

    return 0;
}
