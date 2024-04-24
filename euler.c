#include <stdio.h>
#include <math.h>

void euler(double x0, double y0, double h, double xn) {
    double x = x0;
    double y = y0;
    double yt = (x + 1)*(x + 1) - 0.5*exp(x);
    double slope = y - (x*x) + 1;
    double err = fabs(yt - y);
    double errb = 0.5*h*(0.5*exp(2)-2)*(exp(x)-1);

    while (x < xn) {
        printf("x = %f, y = %f, yt = %f, error = %f, error bound = %f\n", x, y, yt, err, errb);
        
        x = x + h;
        y = y + h * slope;
        yt = (x + 1)*(x + 1) - 0.5*exp(x);
        slope = y - (x*x) + 1;
        err = fabs(yt - y);
        errb = 0.5*h*(0.5*exp(2)-2)*(exp(x)-1);
    }

    printf("x = %f, y = %f, yt = %f, error = %f, error bound = %f\n", x, y, yt, err, errb);
}

int main() {
    double x0 = 0;
    double y0 = 0.5;
    double h = 0.2;
    double xn = 2;

    euler(x0, y0, h, xn);

    return 0;
}
