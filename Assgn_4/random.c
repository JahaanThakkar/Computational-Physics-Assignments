#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 10000
#define MEAN 0.5

double generate_exponential(double mean) {
    double u = (double)rand() / RAND_MAX;
    return -mean * log(1 - u);
}

int main() {
    double numbers[N];
    for (int i = 0; i < N; i++) {
        numbers[i] = generate_exponential(MEAN);
    }

    // Write the numbers to a file
    FILE *file = fopen("exponential_numbers.txt", "w");
    if (file == NULL) {
        fprintf(stderr, "Error opening file.\n");
        return 1;
    }
    for (int i = 0; i < N; i++) {
        fprintf(file, "%f\n", numbers[i]);
    }
    fclose(file);

    return 0;
}

