#include <iostream>
#include <mpi.h>

int main(int argc, char* argv[]) 
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int value = rank + 1;
    int* result = new int[size];

    // Gather the values from all the processes to the root process
    MPI_Gather(&value, 1, MPI_INT, result, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Brodcast the gathered values to all the processes
    MPI_Bcast(result, size, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "Gatherd and broadcasted values across all processes:\n";
        for (int i = 0; i < size; i++) {
            std::cout << "Process " << i << ": " << result[i] << std::endl;
        }
    }

    delete[] result;
    MPI_Finalize();
    return 0;
}