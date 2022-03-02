#include <iostream>

using namespace std;

int main(){
    /*
     * espylon = 2ni = 2^(-t)
     * x- <---> x <---> x+
     * log10((fp(x)-f(x)/f(x)) -|- log10(x)
     */
    int t = 0;
    double epsylon;
    double stala = 1.0;

    epsylon = 1.0;

    while(1 + epsylon > stala) {
        epsylon /= 2.0;
        t += 1;
    }

    cout << t << endl;
    cout << epsylon;

    return 0;
}


