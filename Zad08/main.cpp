//
// Created by volie on 13.04.2022.
//

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#ifndef M_PIl
#define M_PIl (3.14159265358979323846264338327950288)
#endif

using namespace std;

template <typename type> type derivative_function(type x);
template <typename type> type forward_difference(type x, type h);
template <typename type> type backward_difference(type x, type h);
template <typename type> type central_difference(type x, type h);
template <typename type> type three_point_forward_difference(type x, type h);
template <typename type> type three_point_backward_difference(type x, type h);
template <typename type> type compute_epsilon();
template <typename type> void write_data();


int main() {

    write_data<double>();
    write_data<float>();

    system("gnuplot -p -e \"set ylabel 'blad bezwzgledny(double)'; set xlabel 'log10(h)'; plot for [col=2:8] 'double.txt' using 1:col title columnheader w l\"");
    system("gnuplot -p -e \"set ylabel 'blad bezwzgledny(float)'; set xlabel 'log10(h)'; plot for [col=2:8] 'float.txt' using 1:col title columnheader w l\"");

    return 0;
}


template <typename type> type derivative_function(type x) {
    return cos(x);
}

template <typename type> type forward_difference(type x, type h) {
    return (sin(x + h) - sin(x)) / h;
}

template <typename type> type backward_difference(type x, type h) {
    return (sin(x) - sin(x - h)) / h;
}

template <typename type> type central_difference(type x, type h) {
    return (sin(x + h) - sin(x - h)) / (2.0 * h);
}

template <typename type> type three_point_forward_difference(type x, type h) {
    return ((-3.0 / 2.0) * sin(x) + 2 * sin(x + h) - (1.0 / 2.0) * sin(x + 2.0 * h)) / h;
}

template<typename type> type three_point_backward_difference(type x, type h) {
    return ((1.0 / 2.0) * sin(x - 2.0 * h) - 2 * sin(x - h) + (3.0 / 2.0) * sin(x)) / h;
}

template<typename type>
type compute_epsilon() {
    type eps = 1;
    type eps2 = eps;

    while (1 + eps2 > 1) {
        eps = eps2;
        eps2 /= 2;
    }

    return eps;
}

template<typename type>
void write_data() {
    type a = 0.0, mid = M_PI / 4, b = M_PI / 2, h = 0.1;
    type epsilon = compute_epsilon<type>();
    type data[7];
    type blad[2][8];
    type rzad[7];
    ofstream file;
    string name = typeid(type).name();
    name += ".txt";
    int count = 0;

    file.open(name, ios::out);
    if(!file.good())
        return;
    file << scientific << setprecision(6);

    file << "log10(h) " << "forwardDifferenceA " << "forwardDifference3A " <<
          "forwardDifferenceMid " << "centralDifferenceMid " << "backwardDifferenceMid " <<
          "backwardDifferenceB " << "backwardDifference3b" << endl;
    while(h > epsilon) {
        data[0] = log10(fabs(derivative_function(a) - forward_difference(a, h)));
        data[1] = log10(fabs(derivative_function(a) - three_point_forward_difference(a, h)));
        data[2] = log10(fabs(derivative_function(mid) - forward_difference(mid, h)));
        data[3] = log10(fabs(derivative_function(mid) - central_difference(mid, h)));
        data[4] = log10(fabs(derivative_function(mid) - backward_difference(mid, h)));
        data[5] = log10(fabs(derivative_function(b) - backward_difference(b, h)));
        data[6] = log10(fabs(derivative_function(b) - three_point_backward_difference(b, h)));
        file << log10(h) << " ";
        if(count == 0) {
            for(int j = 0; j < 7; j++) {
                blad[0][j] = data[j];
            }
            blad[0][7] = log10(h);
        }
        if(count == 50) {
            for(int j = 0; j < 7; j++) {
                blad[1][j] = data[j];
            }
            blad[1][7] = log10(h*1.1);
        }
        for(int i = 0; i < 7; i++) {
            file << data[i] << " ";
        }
        file << endl;
        h /= 1.1;
        count++;
    }

    for(int i = 0; i < 7; i++) {
        rzad[i] = (blad[1][i] - blad[0][i]) / (blad[1][7] - blad[0][7]);
    }

    cout << typeid(type).name() << endl;
    cout << "rzad forwardDifferenceA: " << rzad[0] << "\t\trzad teoretyczny: "<< 1 << endl;
    cout << "rzad forwardDifference3A: " << rzad[1] << "\t\trzad teoretyczny: "<< 2 << endl;
    cout << "rzad forwardDifferenceMid: " << rzad[2] << "\t\trzad teoretyczny: "<< 1 << endl;
    cout << "rzad centralDifferenceMid: " << rzad[3] << "\t\trzad teoretyczny: "<< 2 << endl;
    cout << "rzad backwardDifferenceMid: " << rzad[4] << "\t\trzad teoretyczny: "<< 1 << endl;
    cout << "rzad backwardDifferenceB: " << rzad[5] << "\t\trzad teoretyczny: "<< 1 << endl;
    cout << "rzad backwardDifference3b: " << rzad[6] << "\t\trzad teoretyczny: "<< 2 << endl << endl;

}
