#include <iostream>
#include <mpi.h>

int main(int argc, char* argv[]) 
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int send_data = rank + 1;
    int* recv_data = new int[size];

    // Alltoall communication to exchange data between all processes
    MPI_Alltoall(&send_data, 1, MPI_INT, recv_data, 1, MPI_INT, MPI_COMM_WORLD);

    std::cout << "Process " << rank << " received data: ";
    for (int i = 0; i < size; i++) {
        std::cout << recv_data[i] << " ";
    }
    std::cout << std::endl;

    delete[] recv_data;
    MPI_Finalize();

    return 0;
}