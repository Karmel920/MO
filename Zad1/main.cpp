#include <iostream>

using namespace std;

int main(){
    /*
     * espylon = 2ni = 2^(-t)
     * x- <---> x <---> x+
     * log10((fp(x)-f(x)/f(x)) -|- log10(x)
     */

    int t = -1;
    double epsylon;
    double stala = 1.0;

    epsylon = 1.0;

    double a;
    a = 1.0 + epsylon;

    while(a > stala) {
        epsylon /= 2.0;
        t += 1;
        a = 1.0f + epsylon;
    }

    cout << "liczba bitow: " << t << endl;
    cout << "epsilon: " << epsylon * 2;

    return 0;
}


