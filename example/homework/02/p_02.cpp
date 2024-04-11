#include <iostream>
#include <string>
#include <cmath>

#include <mpi.h>

class Domain {
    public:
        // Constructor
        Domain(int xStart, int xlen, int yStart, int ylen, const std::string &name) 
        : xStart(xStart), xlen(xlen), yStart(yStart), ylen(ylen), domain(new char[(xlen+1)*(ylen+1)]), name(name) {}
        // Destructor
        ~Domain() { delete[] domain; }
        char &operator()(int i, int j) { return domain[i*(ylen+1) + j]; }
        char operator()(int i, int j) const { return domain[i*(ylen+1) + j]; }

        int get rows() const { return xlen;}
        int get cols() const { return ylen;}
        const std::string &get_name() const { return name; }

        char *rawprt() { return domain; }
    private:
        int xStart, yStart;
        int xlen, ylen;
        char *domain;

        std::string name;
};

void zero_domain(Domain &d) {
    for (int i = 0; i < d.rows(); i++) {
        for (int j = 0; j < d.cols(); j++) {
            d(i, j) = 0;
        }
    }
}

void print_domain(Domain &d, int rank);
void update_domain(Domain &new_domain, Domai &old_domain, int size, int myrank, MPI_Comm comm);
void parallel_code(int M, int N, int iterations, int size, int myrank, MPI_Comm comm);

int main(int argc, char **argv)
{
    int M, N;
    int iterations;

    if(argc < 4){
        std::cout << "Usage: " << argv[0] << " <M> <N> <iterations>" << std::endl;
        exit(0);
    }
    MPI_Init(&argc, &argv);

    int size, myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int array[3];
    if (myrank == 0){
        M = atoi(argv[1]);
        N = atoi(argv[2]);
        iterations = atoi(argv[3]);
        array[0] = M;
        array[1] = N;
        array[2] = iterations;
    }
    MPI_Bcast(array, 3, MPI_INT, 0, MPI_COMM_WORLD);
    if(myrank != 0){
        M = array[0];
        N = array[1];
        iterations = array[2];
    }

    parallel_code(M, N, iterations, size, myrank, MPI_COMM_WORLD);

    MPI_Finalize();
    return 0;
}

void parallel_code(int M, int N, int iterations, int size, int myrank, MPI_Comm comm)
{
    int numSubM = sqrt(size);
    int numSubN = sqrt(size);
    int subM = M / numSubM;
    int subN = N / numSubN;
    int subMStart = (myrank / numSubM) * subM;
    int subNStart = (myrank / numSubN) * subN;

    Domain old_domain(subMStart, subM, subNStart, subN, "old");
    Domain new_domain(subMStart, subM, subNStart, subN, "new");

    zero_domain(old_domain);
    zero_domain(new_domain);

    Domain *old, *new;
    old = &old_domain;
    new = &new_domain;

    for(int i=0; i < iterations; i++){
        update_domain(*new, *old, size, myrank, comm);
        Domain *temp = old;
        old = new;
        new = temp;
    }
}

void print_domain(Domain &d, int rank)
{
    std::cout << "Rank: " << rank << " Domain: " << d.get_name() << std::endl;
    for (int i = 0; i < d.rows(); i++) {
        for (int j = 0; j < d.cols(); j++) {
            std::cout << (d(i,j) ? "*" : ".");
        }
        std::cout << std::endl;
    }
}

inline char update_the_cell(char cell, int neighbor_count)
{
    char newcell;
    if(cell==0) // Dead cell
    {
        if(neighbor_count==3)
            newcell = 1;
        else
            newcell = 0;
    }
    else // Live cell
    {
        if(neighbor_count==2 || neighbor_count==3)
            newcell = 1;
        else
            newcell = 0;
    }
}

void update_domain(Domain &new_domain, Domain &old_domain, int size, int myrank, MPI_Comm comm)
{
    MPI_Request request[4];

    int m = new_domain.rows();
    int n = new_domain.cols();

    char *top_row = new char[n];
    char *right_col = new char[m];
    char *bottom_row = new char[n];
    char *left_col = new char[m];

    char *top_halo = new char[n];
    char *right_halo = new char[m];
    char *bottom_halo = new char[n];
    char *left_halo = new char[m];

    const int top_row_index = 0;
    const int bottom_row_index = m-1;
    const int left_col_index = 0;
    const int right_col_index = n-1;

    const int TOP_HALO = 0, BOTTOM_HALO = 1, LEFT_HALO = 2, RIGHT_HALO = 3;

    // post receives for the halos from the neighbors
    MPI_Irecv(top_halo, n, MPI_CHAR, (myrank-1), BOTTOM_HALO, comm, &request[0]);
    MPI_Irecv(right_halo, m, MPI_CHAR, (myrank+numSubM), 3, comm, &request[3]);
    MPI_Irecv(bottom_halo, n, MPI_CHAR, (myrank+1), TOP_HALO, comm, &request[1]);
    MPI_Irecv(left_halo, m, MPI_CHAR, (myrank-numSubM), 2, comm, &request[2]);

    // gather and send the halos to the neighbors
    for(int i=0; i<m; i++)
    {
        top_row[i] = old_domain(top_row_index, i);
    }
    MPI_Isend(top_row, n, MPI_CHAR, (myrank-1+size)%size, BOTTOM_HALO, comm, &request[0]);
    for(int i=0; i<m; i++)
    {
        bottom_row[i] = old_domain(bottom_row_index, i);
    }
    MPI_Isend(bottom_row, n, MPI_CHAR, (myrank+1)%size, TOP_HALO, comm, &request[1]);
    for(int i=0; i<n; i++)
    {
        left_col[i] = old_domain(i, left_col_index);
    }
    MPI_Isend(left_col, m, MPI_CHAR, (myrank-numSubM+size)%size, RIGHT_HALO, comm, &request[2]);
    for(int i=0; i<n; i++)
    {
        right_col[i] = old_domain(i, right_col_index);
    }
    MPI_Isend(right_col, m, MPI_CHAR, (myrank+numSubM)%size, LEFT_HALO, comm, &request[3]);

    // complete all 8 transfers
    MPI_Waitall(8, request, MPI_STATUSES_IGNORE);

    // update the interior of the domain
    for(int i=1; i<m-1; i++)
    {
        for(int j=1; j<n-1; j++)
        {
            int neighbor_count = 0;
            neighbor_count += top_row[j-1];
            neighbor_count += top_row[j];
            neighbor_count += top_row[j+1];
            neighbor_count += old_domain(i, j-1);
            neighbor_count += old_domain(i, j+1);
            neighbor_count += bottom_row[j-1];
            neighbor_count += bottom_row[j];
            neighbor_count += bottom_row[j+1];

            new_domain(i, j) = update_the_cell(old_domain(i, j), neighbor_count);
        }
    }
}