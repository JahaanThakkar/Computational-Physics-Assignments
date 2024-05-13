#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

double gaussian(double x) {
    return exp(-x*x);
}

int main() {
    int N = 10001;
    double x0 = -5000;
    fftw_complex signal[N], FT[N];
    fftw_plan dft;
    double dx = 2 * x0 / (N - 1);

    for (int i=0; i < N; i++) {
        signal[i] = gaussian(-x0 + i * dx) + I * 0.0;
    }

    dft = fftw_plan_dft_1d(N, signal, FT, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(dft);

    FILE *file = fopen("q4_fft.csv","w");
    if (file == NULL) {
        printf("Error Opening File!\n");
        return 1;
    }

    for (int i = 0; i < N; i++) {
        fprintf(file, "%g, %g\n", creal(FT[i]), cimag(FT[i]));
    }

    fclose(file);

    fftw_destroy_plan(dft);
    fftw_cleanup();

    printf("Fourier transform data has been written to 'q4_fft.csv'.\n");

    return 0;
}

