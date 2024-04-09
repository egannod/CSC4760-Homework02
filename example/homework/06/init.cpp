#include <Kokkos_Core.hpp>
#include <mpi.h>
#include <cstdio>
#include <iostream>

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  Kokkos::initialize(argc, argv);
  {
  // Make View and init values
  Kokkos::View<int*> view_0("view_0", 10);
  Kokkos::parallel_for("initialize_view_0", Kokkos::RangePolicy<>(0, 10), KOKKOS_LAMBDA(const int i) {
    view_0(i) = i*i;
  });
  Kokkos::fence();

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int send_rank = 0;
  int recv_rank = 1;
 
  // Send View values with MPI functions
  if (rank == send){
    MPI_Send(view_0.data(), view_0.size(), MPI_INT, recv_rank, 0, MPI_COMM_WORLD);
  } else if (rank == recv_rank){
    Kokkos::View<int*> recv_view("view_0", 10);
    MPI_Recv(recv_view.data(), recv_view.size(), MPI_INT, send_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Print received values
    Kokkos::parallel_for("print_recieved_view_0", Kokkos::RangePolicy<>(0, recv_view.size()), KOKKOS_LAMBDA(const int i) {
      std::cout << "recv_view[" << i << "] = " << recv_view[i] << std::endl;
    });
    Kokkos::fence();
  }

  // Output
  printf("\nhello world\n");
  }
  Kokkos::finalize();
  MPI_Finalize();

  return 0;
}
