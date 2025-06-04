#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
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

    int num_threads = 8;  // liczba w¹tków
    omp_set_num_threads(num_threads);

    // Proces 0 wczytuje liczby
    int* numbers = NULL;
    int count = 0;

    if (rank == 0) {
        FILE* fin = fopen("../input.txt", "r");
        if (!fin) {
            printf("Nie uda³o siê otworzyæ pliku.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        numbers = (int*)malloc(100000000 * sizeof(int));
        while (fscanf(fin, "%d,", &numbers[count]) == 1)
            count++;
        fclose(fin);

        printf("Wczytano %d liczb.\n", count);
    }

    // Rozeslij liczbe elementow
    MPI_Bcast(&count, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Podzia³ pracy
    int chunk_size = count / size;
    int remainder = count % size;
    int my_count = (rank < remainder) ? chunk_size + 1 : chunk_size;

    int* my_numbers = (int*)malloc(my_count * sizeof(int));
    int* sendcounts = NULL;
    int* displs = NULL;

    if (rank == 0) {
        sendcounts = (int*)malloc(size * sizeof(int));
        displs = (int*)malloc(size * sizeof(int));
        int offset = 0;
        for (int i = 0; i < size; i++) {
            sendcounts[i] = (i < remainder) ? chunk_size + 1 : chunk_size;
            displs[i] = offset;
            offset += sendcounts[i];
        }
    }

    // Rozeslij liczby
    MPI_Scatterv(numbers, sendcounts, displs, MPI_INT,
        my_numbers, my_count, MPI_INT,
        0, MPI_COMM_WORLD);

    double start_time = MPI_Wtime();

    // Lokalna tablica liczb pierwszych
    int* local_primes = (int*)malloc(my_count * sizeof(int));
    int local_count = 0;

#pragma omp parallel
    {
        int* thread_primes = (int*)malloc(my_count * sizeof(int));
        int thread_count = 0;

        int i;
#pragma omp for
        for (i = 0; i < my_count; i++) {
            if (is_prime(my_numbers[i]))
                thread_primes[thread_count++] = my_numbers[i];
        }

#pragma omp critical
        {
            for (int j = 0; j < thread_count; j++)
                local_primes[local_count++] = thread_primes[j];
        }

        free(thread_primes);
    }

    // Zbierz wyniki
    int* recv_counts = NULL;
    int* recv_displs = NULL;
    int* all_primes = NULL;

    if (rank == 0) {
        recv_counts = (int*)malloc(size * sizeof(int));
    }

    MPI_Gather(&local_count, 1, MPI_INT,
        recv_counts, 1, MPI_INT,
        0, MPI_COMM_WORLD);

    if (rank == 0) {
        int total_primes = 0;
        recv_displs = (int*)malloc(size * sizeof(int));
        recv_displs[0] = 0;
        for (int i = 0; i < size; i++) {
            total_primes += recv_counts[i];
            if (i > 0)
                recv_displs[i] = recv_displs[i - 1] + recv_counts[i - 1];
        }
        all_primes = (int*)malloc(total_primes * sizeof(int));
    }

    MPI_Gatherv(local_primes, local_count, MPI_INT,
        all_primes, recv_counts, recv_displs, MPI_INT,
        0, MPI_COMM_WORLD);

    double end_time = MPI_Wtime();

    if (rank == 0) {
        FILE* fout = fopen("output_mpi_openmp.txt", "w");
        for (int i = 0; i < recv_displs[size - 1] + recv_counts[size - 1]; i++)
            fprintf(fout, "%d,", all_primes[i]);
        fclose(fout);

        printf("Czas wykonania (MPI+OpenMP, %d procesy, %d watki/proces): %.4f sekund\n",
            size, num_threads, end_time - start_time);

        free(numbers);
        free(sendcounts);
        free(displs);
        free(recv_counts);
        free(recv_displs);
        free(all_primes);
    }

    free(my_numbers);
    free(local_primes);

    MPI_Finalize();
    return 0;
}
