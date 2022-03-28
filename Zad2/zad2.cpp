// MO_Lab2.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <math.h>

using namespace std;

double szereg(int x) {
	double y = 1.0;
	double n = 1.0, m = 1.0;
	for (double i = 1; i <= 500; i++)
	{
		m *= i;
		y += ((pow((double)x, (double)i))  / (double)m);
	}
	return y;
}


double liczenie(double x) {
	 double z = 1, y;
	 for (double i = 0; i < x; i++) {
		 z *= szereg(1);
	 }
	 y = 1 / z;
	 return y;

}

double liczenie2(double x) {
	 double z=1, y;
	 for (double i = 0; i < x; i++) {
		 z *= szereg(1);
	 }
	
	 return z;
}

int main(int argc, char *argv[])
{
	cout << "x\tWartosc dokladna\tszereg Taylora\t\tBlad wzgledny\n" << endl;
	double blad_wzgledny,blad_bezwzgledny;
	double t=1;
	for (double i = -30.0; i <= 30.0; i+=1)
	{
		double z = exp(i);
		if (i < 0)
		{
			t = liczenie(fabs(i));
		}
		else if (i>0)
		{
			t = liczenie2(i);
		}
		else
			t = szereg(0);
		blad_wzgledny = (t - z) / z;

		blad_bezwzgledny = t - z;
		//cout << "\tBlad wzgledny : "<<blad_wzgledny << "\t\tBlad bezwzgledny :" << blad_bezwzgledny << endl;
		//double wart_konc = o + blad_wzgledny;

		double wart_konc = t + blad_bezwzgledny;
		if (i >= 0 && i < 14)
		{
			cout << i << "\t" << z << "\t\t\t" << wart_konc<< "\t\t\t" << blad_wzgledny  << endl;
		}
		else
			cout << i << "\t" << z << "\t\t" << wart_konc << "\t\t" << blad_wzgledny << endl;
		
		//printf("%d \t\t %.8lf \t\t %.8lf\n", i, z, wart_konc);
	}


	system("PAUSE");
	return 0;
}

