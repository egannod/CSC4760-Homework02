#include <iostream>
#include "mpi.h"

int main(int argc, char* argv[]) 
{
    int rank, size;
    int ring_num;
    int N;

    int stop = -1;
    if (argc < 2) {
        std::cout << "Please input the times around the ring" << std::endl;
        return 1;
    }
    else {
        N = atoi(argv[1]);
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    while(1)
    {
        if (rank == 0) {
            MPI_Send(&ring_num, 1, MPI_INT, 2, 0, MPI_COMM_WORLD);
            if (ring_num/(size-1) == N) {
                MPI_Send(&stop, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
                break;
            }
            MPI_Recv(&ring_num, 1, MPI_INT, size-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            ring_num++;
        }
        else {
            MPI_Recv(&ring_num, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (ring_num == -1) {
                MPI_Send(&stop, 1, MPI_INT, (rank+1)%size, 0, MPI_COMM_WORLD);
                break;
            }
            ring_num++;
            MPI_Send(&ring_num, 1, MPI_INT, (rank+1)%size, 0, MPI_COMM_WORLD);
        }
    }
    if (rank == 0) {
        std::cout << "The ring has been around " << N << " times" << std::endl;
        std::cout << "end number: " << ring_num << std::endl;
    }
    MPI_Finalize();

    return 0;
}

