#include <iostream>
#include <cmath>

using namespace std;

double convert(double degree){
    double pi = 3.14159265359;
    return (degree*(pi/180));
}


int main()
{
    double  result3;
//    result = tgamma(1 + 1.5);
//    result2 = tgamma((1+1.5)/2);

    result3 =pow((tgamma(2.5)*sin(135)/tgamma(1.25)*1.5*pow(2,0.25)),0.66666667);
    cout<<"ANSWER: "<<result3<<endl;

    return 0;
}
