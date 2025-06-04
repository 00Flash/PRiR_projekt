// mpi_sito.c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

int is_prime(int n) {
    if (n < 2) return 0;
    for (int i = 2; i <= sqrt(n); i++)
        if (n % i == 0)
            return 0;
    return 1;
}

int main(int argc, char** argv) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int* numbers = NULL;
    int count = 0;

    // Proces 0 wczytuje liczby z pliku
    if (rank == 0) {
        FILE* fin = fopen("../input.txt", "r");
        if (!fin) {
            printf("Nie uda³o siê otworzyæ pliku.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        int tmp;
        while (fscanf(fin, "%d,", &tmp) == 1)
            count++;
        rewind(fin);

        numbers = malloc(count * sizeof(int));
        for (int i = 0; i < count; i++)
            fscanf(fin, "%d,", &numbers[i]);
        fclose(fin);
    }

    // Rozsy³amy liczbê elementów do wszystkich procesów
    MPI_Bcast(&count, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Obliczanie sendcounts i displs
    int chunk = count / size;
    int rem = count % size;

    int* sendcounts = NULL;
    int* displs = NULL;
    if (rank == 0) {
        sendcounts = malloc(size * sizeof(int));
        displs = malloc(size * sizeof(int));
        for (int i = 0; i < size; i++)
            sendcounts[i] = chunk + (i < rem ? 1 : 0);
        displs[0] = 0;
        for (int i = 1; i < size; i++)
            displs[i] = displs[i - 1] + sendcounts[i - 1];
    }

    int local_count = chunk + (rank < rem ? 1 : 0);
    int* local_numbers = malloc(local_count * sizeof(int));

    // Rozsy³anie danych do procesów
    MPI_Scatterv(numbers, sendcounts, displs, MPI_INT,
        local_numbers, local_count, MPI_INT,
        0, MPI_COMM_WORLD);

    double start_time = MPI_Wtime();

    // Sprawdzanie liczb pierwszych i zapisywanie do pliku
    char filename[64];
    sprintf(filename, "output_mpi_rank%d.txt", rank);
    FILE* fout = fopen(filename, "w");

    for (int i = 0; i < local_count; i++)
        if (is_prime(local_numbers[i]))
            fprintf(fout, "%d,", local_numbers[i]);

    double end_time = MPI_Wtime();
    printf("Rank %d: czas wykonania (MPI): %.4f sekund\n", rank, end_time - start_time);

    fclose(fout);

    // Zwolnienie pamiêci
    if (rank == 0) {
        free(numbers);
        free(sendcounts);
        free(displs);
    }
    free(local_numbers);

    MPI_Finalize();
    return 0;
}
