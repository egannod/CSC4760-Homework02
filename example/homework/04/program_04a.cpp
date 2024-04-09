#include <iostream>
#include <mpi.h>

int main(int argc, char* argv[]) 
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int value = rank + 1;
    int result;

    // Broadcast the value to all processes
    MPI_Bcast(&value, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Reduce the values from all the processes to the root process
    MPI_Reduce(&value, &result, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "The sum of all values is: " << result << std::endl;
    }

    MPI_Finalize();
    return 0;
}