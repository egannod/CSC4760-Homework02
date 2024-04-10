#include <iostream>
#include <mpi.h>

int main(int argc, char** argv)
{
    int rank, size;
    int color1, color2;
    MPI_Comm new_comm1, new_comm2;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Choose P and Q such that P*Q = size
    int P = 2;
    int Q = 2;

    // Check if the world size is compatible
    if (size != P*Q) {
        std::cout << "The world size is not compatible with P and Q" << std::endl;
        MPI_Finalize();
        return 1;
    }

    // First split: color based on rank divided by Q
    color1 = rank / Q;
    MPI_Comm_split(MPI_COMM_WORLD, color1, rank, &new_comm1);

    // Second split: color based on rank modulo Q
    color2 = rank % Q;
    MPI_Comm_split(new_comm1, color2, rank, &new_comm2);

    // Finalize MPI
    MPI_Finalize();

    return 0;
}