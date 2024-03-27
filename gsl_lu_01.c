#include <stdio.h>
#include <gsl/gsl_linalg.h>

void lu_decomposition(double A[], int size, double L[], double U[]) {
    gsl_matrix_view gsl_A = gsl_matrix_view_array(A, size, size);
    gsl_permutation *p = gsl_permutation_alloc(size);
    int signum;

    gsl_linalg_LU_decomp(&gsl_A.matrix, p, &signum);

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (i > j)
                L[i * size + j] = gsl_matrix_get(&gsl_A.matrix, i, j);
            else if (i == j)
                L[i * size + j] = 1.0;
            else
                L[i * size + j] = 0.0;

            if (i <= j)
                U[i * size + j] = gsl_matrix_get(&gsl_A.matrix, i, j);
            else
                U[i * size + j] = 0.0;
        }
    }

    gsl_permutation_free(p);
}

int main() {
    double A1[] = {3, -1, 1, 3, 6, 2, 3, 3, 7};
    double A2[] = {10, -1, 0, -1, 10, -2, 0, -2, 10};
    double A3[] = {10, 5, 0, 0, 5, 10, -4, 0, 0, -4, 8, -1, 0, 0, -3, 5};
    double A4[] = {4, 1, 1, 0, 1, -1, -3, 1, 1, 0, 2, 1, 5, -1, -1, -1, -1 -1, 4, 0, 0, 2, -1, 1, 4};

    // Allocate space for L and U matrices
    double L1[3 * 3];
    double U1[3 * 3];
    double L2[3 * 3];
    double U2[3 * 3];
    double L3[4 * 4];
    double U3[4 * 4];
    double L4[5 * 5];
    double U4[5 * 5];

    // Perform LU decomposition
    lu_decomposition(A1, 3, L1, U1);
    lu_decomposition(A2, 3, L2, U2);
    lu_decomposition(A3, 4, L3, U3);
    lu_decomposition(A4, 5, L4, U4);

    // Print the results
    printf("Matrix L1:\n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%8.7f ", L1[i * 3 + j]);
        }
        printf("\n");
    }

    printf("Matrix U1:\n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%8.7f ", U1[i * 3 + j]);
        }
        printf("\n");
    }
    printf("Matrix L2:\n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%8.7f ", L2[i * 3 + j]);
        }
        printf("\n");
    }

    printf("Matrix U2:\n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%8.7f ", U2[i * 3 + j]);
        }
        printf("\n");
    }
    printf("Matrix L3:\n");
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            printf("%8.7f ", L3[i * 4 + j]);
        }
        printf("\n");
    }

    printf("Matrix U3:\n");
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            printf("%8.7f ", U3[i * 4 + j]);
        }
        printf("\n");
    }
    printf("Matrix L4:\n");
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            printf("%8.7f ", L4[i * 5 + j]);
        }
        printf("\n");
    }

    printf("Matrix U4:\n");
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            printf("%8.7f ", U4[i * 5 + j]);
        }
        printf("\n");
    }

    // Verification of LU decompsition


    // For A1
    gsl_matrix_view L1M = gsl_matrix_view_array(L1, 3, 3);
    gsl_matrix_view U1M = gsl_matrix_view_array(U1, 3, 3);

    // Defining matrix R1
    gsl_matrix *R1 = gsl_matrix_alloc(3, 3);

    // Product LU = R1
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &L1M.matrix, &U1M.matrix, 0.0, R1);

    // Printing R1
    printf("Result of matrix multiplication:\n");
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            printf("%8.3f ", gsl_matrix_get(R1, i, j));
        }
        printf("\n");
    }

    // Empty R1
    gsl_matrix_free(R1);


    // For A2
    gsl_matrix_view L2M = gsl_matrix_view_array(L2, 3, 3);
    gsl_matrix_view U2M = gsl_matrix_view_array(U2, 3, 3);

    // Defining matrix R2
    gsl_matrix *R2 = gsl_matrix_alloc(3, 3);

    // Product L1U1 = R2
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &L2M.matrix, &U2M.matrix, 0.0, R2);

    // Printing R2
    printf("Result of matrix multiplication:\n");
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            printf("%8.3f ", gsl_matrix_get(R2, i, j));
        }
        printf("\n");
    }

    // Empty R2
    gsl_matrix_free(R2);


    // For A3
    gsl_matrix_view L3M = gsl_matrix_view_array(L3, 4, 4);
    gsl_matrix_view U3M = gsl_matrix_view_array(U3, 4, 4);

    // Defining matrix R3
    gsl_matrix *R3 = gsl_matrix_alloc(4, 4);

    // Product LU = R3
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &L3M.matrix, &U3M.matrix, 0.0, R3);

    // Printing R3
    printf("Result of matrix multiplication:\n");
    for (size_t i = 0; i < 4; ++i) {
        for (size_t j = 0; j < 4; ++j) {
            printf("%8.3f ", gsl_matrix_get(R3, i, j));
        }
        printf("\n");
    }

    // Empty R
    gsl_matrix_free(R3);


    // For A4
    gsl_matrix_view L4M = gsl_matrix_view_array(L4, 5, 5);
    gsl_matrix_view U4M = gsl_matrix_view_array(U4, 5, 5);

    // Defining matrix R4
    gsl_matrix *R4 = gsl_matrix_alloc(5, 5);

    // Product LU = R4
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &L4M.matrix, &U4M.matrix, 0.0, R4);

    // Printing R4
    printf("Result of matrix multiplication:\n");
    for (size_t i = 0; i < 5; ++i) {
        for (size_t j = 0; j < 5; ++j) {
            printf("%8.3f ", gsl_matrix_get(R4, i, j));
        }
        printf("\n");
    }

    // Empty R
    gsl_matrix_free(R4);

    return 0;
}