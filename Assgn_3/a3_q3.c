#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

double sinc(double x) {
    if (x == 0.0){
        return 1.0;
    } else {
        return sin(x)/x;
    }
}

int main() {
    int N = 10001;
    double x0 = -5000;
    double dx = 2 * x0 / (N - 1);
    double data[2*N];
    
    gsl_fft_complex_wavetable *wt;
    gsl_fft_complex_workspace *ws;

    for (int i=0; i < N; i++) {
        REAL(data,i) = sinc(-x0 + i * dx);
        IMAG(data,i) = 0.0;
    }

    wt = gsl_fft_complex_wavetable_alloc(N);
    ws = gsl_fft_complex_workspace_alloc(N);
    
    gsl_fft_complex_forward(data,1,N,wt,ws);

    FILE *file = fopen("q3_fft.csv","w");
    if (file == NULL) {
        printf("Error Opening File!\n");
        return 1;
    }

    for (int i = 0; i < N; i++) {
        fprintf(file, "%g, %g\n", REAL(data,i), IMAG(data,i));
    }
    
    gsl_fft_complex_wavetable_free (wt);
    gsl_fft_complex_workspace_free (ws);
    fclose(file);

    printf("Fourier transform data has been written to 'q3_fft.csv'.\n");

    return 0;
}

