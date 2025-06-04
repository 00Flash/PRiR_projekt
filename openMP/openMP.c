#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

int is_prime(int n) {
    if (n < 2) return 0;
    for (int i = 2; i <= sqrt(n); i++)
        if (n % i == 0)
            return 0;
    return 1;
}

int main() {
    int num_threads = 8; // Ustaw liczbê w¹tków
    omp_set_num_threads(num_threads);

    // Otwórz pliki
    FILE* fin = fopen("../input.txt", "r");
    FILE* fout = fopen("output_openmp.txt", "w");
    if (!fin || !fout) {
        printf("Nie uda³o siê otworzyæ pliku.\n");
        return 1;
    }

    // Wczytaj liczby
    int max_numbers = 100000000; // Dostosuj do potrzeb
    int* numbers = (int*)malloc(max_numbers * sizeof(int));
    int count = 0;
    while (fscanf(fin, "%d,", &numbers[count]) == 1)
        count++;
    fclose(fin);

    printf("Wczytano %d liczb\n", count);

    double start_time = omp_get_wtime();

    // Przygotuj do obliczeñ
    int* primes = (int*)malloc(count * sizeof(int));

    int real_num_threads;
#pragma omp parallel
    {
#pragma omp single
        real_num_threads = omp_get_num_threads();
    }

    int* local_counts = (int*)malloc(real_num_threads * sizeof(int));
    int* offsets = (int*)malloc(real_num_threads * sizeof(int));

#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int* local_primes = (int*)malloc(count * sizeof(int)); // lokalna tablica w¹tku
        int local_count = 0;

        int j;
#pragma omp for
        for (j = 0; j < count; j++) {
            if (is_prime(numbers[j])) {
                local_primes[local_count++] = numbers[j];
            }
        }

        local_counts[tid] = local_count;

#pragma omp barrier

        // Oblicz offsety
#pragma omp single
        {
            offsets[0] = 0;
            for (int i = 1; i < real_num_threads; i++)
                offsets[i] = offsets[i - 1] + local_counts[i - 1];
        }

        // Skopiuj lokalne wyniki do globalnej tablicy
        for (int j = 0; j < local_count; j++)
            primes[offsets[tid] + j] = local_primes[j];

        free(local_primes);
    }

    // Oblicz ca³kowit¹ liczbê liczb pierwszych
    int prime_count = 0;
    for (int i = 0; i < real_num_threads; i++)
        prime_count += local_counts[i];

    // Bezpieczne sprawdzenie
    if (prime_count > count) {
        printf("B³¹d: prime_count (%d) > count (%d)!\n", prime_count, count);
        exit(1);
    }

    // Zapisz wyniki
    for (int i = 0; i < prime_count; i++)
        fprintf(fout, "%d,", primes[i]);

    double end_time = omp_get_wtime();
    printf("Znaleziono %d liczb pierwszych.\n", prime_count);
    printf("Czas wykonania (OpenMP): %.4f sekund\n", end_time - start_time);

    // Zwolnij pamiêæ
    fclose(fout);
    free(primes);
    free(numbers);
    free(local_counts);
    free(offsets);

    return 0;
}
