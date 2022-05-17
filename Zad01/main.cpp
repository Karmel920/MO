#include <iostream>

using namespace std;

int main(){
    /*
     * espylon = 2ni = 2^(-t)
     * x- <---> x <---> x+
     * log10((fp(x)-f(x)/f(x)) -|- log10(x)
     */

    int t = -1;
    float epsylon = 1.0f;
    float stala = 1.0f;

    while(1.0f + epsylon > stala) {
        epsylon /= 2.0f;
        t += 1;

    }

    cout << "liczba bitow: " << t << endl;
    cout << "epsilon: " << epsylon * 2;

    return 0;
}


