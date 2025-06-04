// sekwencyjne_sito.c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int is_prime(int n) {
    if (n < 2) return 0;
    for (int i = 2; i <= sqrt(n); i++)
        if (n % i == 0)
            return 0;
    return 1;
}

int main() {
    FILE* fin = fopen("../input.txt", "r");
    FILE* fout = fopen("output_sekwencyjny.txt", "w");
    if (!fin || !fout) {
        printf("Nie udało się otworzyć pliku.\n");
        return 1;
    }

    clock_t start = clock();

    int num;
    while (fscanf(fin, "%d,", &num) == 1) {
        if (is_prime(num))
            fprintf(fout, "%d,", num);
    }

    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Czas wykonania (sekwencyjny): %.4f sekund\n", elapsed);

    fclose(fin);
    fclose(fout);
    return 0;
}
