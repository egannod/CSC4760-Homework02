#include <Kokkos_Core.hpp>
#include <cstdio>
#include <iostream>

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);
  {
  // Make View of length n > 1000 and add values
  const int n = 1000;

  Kokkos::View<int*> view_0("view_0", n);
  Kokkos::parallel_for("initialize_view_0", n, KOKKOS_LAMBDA(const int i) {
    view_0(i) = i*i;
  });
  
  // create two additional views of same size and datatype
  Kokkos::View<int*> view_1("view_1", n);
  Kokkos::View<int*> view_2("view_2", n);
 
  // deep_copy
  Kokkos::Timer timer;
  Kokkos::deep_copy(view_1, view_0);
  Kokkos::fence();
  double time1 = timer.seconds();
  timer.reset();
  // user copy
  Kokkos::parallel_for("copy_view_0_to_view_2", n, KOKKOS_LAMBDA(const int i) {
    view_2(i) = view_0(i);
  });
  Kokkos::fence();
  double time2 = timer.seconds();
  
  // Output times 
  std::cout << "Deep copy time: " << time1 << std::endl;
  std::cout << "Parallel_for copy time: " << time2 << std::endl;
  }
  Kokkos::finalize();
}
