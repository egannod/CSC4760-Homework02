// Parallelized Game of Life 2D Decomposition
// I could not complete fully. Did not complete sending and receiving the halo cells and updating based on that data.
#include <iostream>
#include <vector>
#include <unistd.h>
#include <mpi.h>

void parallel_code(int X, int Y, int iterations, int size, int rank, MPI_Comm comm);

void print_board(std::vector<std::vector<int>>& board)
{
    for (int x = 0; x < board.size(); x++)
    {
        std::cout << " ";
        for (int y = 0; y < board[0].size(); y++)
        {
            if (board[x][y] == 0) 
                std::cout << ".";
            else
                std::cout << "*";
        }
        std::cout << std::endl;
    }
}

bool isAlive(std::vector<std::vector<int>>& board, const int x, const int y, std::vector<int> top_row, std::vector<int> bottom_row, std::vector<int> left_column, std::vector<int> right_column)
{
    int X = board.size();
    int Y = board[0].size();
    int alive = 0;
    // testing left
    if (board[(x-1+X)%X][y] == 1) alive++;
    // testing right
    if (board[(x+1)%X][y] == 1) alive++;
    // testing up
    if (board[x][(y-1+Y)%Y] == 1) alive++;
    // testing down
    if (board[x][(y+1)%Y] == 1) alive++;
    // testing up-left
    if (board[(x-1+X)%X][(y-1+Y)%Y] == 1) alive++;
    // testing up-right
    if (board[(x+1)%X][(y-1+Y)%Y] == 1) alive++;
    // testing down-left
    if (board[(x-1+X)%X][(y+1)%Y] == 1) alive++;
    // testing down-right
    if (board[(x+1)%X][(y+1)%Y] == 1) alive++;

    if (board[x][y] == 1)
    {
        if (alive < 2) return false;
        if (alive > 3) return false;
        return true;
    }
    else
    {
        if (alive == 3) return true;
        return false;
    }
}

std::vector<std::vector<int>> random_board(int X, int Y)
{
    std::vector<std::vector<int>> board(X, std::vector<int>(Y, 0));
    for (int x = 0; x < X; x++)
    {
        for (int y = 0; y < Y; y++)
        {
            board[x][y] = rand() % 2;
        }
    }
    return board;
}

int main(int argc, char **argv)
{
    int X, Y;
    int iterations;

    if(argc < 4){
        std::cout << "Usage: " << argv[0] << " <X> <Y> <iterations>" << std::endl;
        exit(0);
    }
    MPI_Init(&argc, &argv);
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int array[3];
    if (rank == 0){
        X = atoi(argv[1]);
        Y = atoi(argv[2]);
        iterations = atoi(argv[3]);
        array[0] = X;
        array[1] = Y;
        array[2] = iterations;
    }
    MPI_Bcast(array, 3, MPI_INT, 0, MPI_COMM_WORLD);
    if(rank != 0){
        X = array[0];
        Y = array[1];
        iterations = array[2];
    }

    parallel_code(X, Y, iterations, size, rank, MPI_COMM_WORLD);

    MPI_Finalize();
    return 0;
}

void parallel_code(int X, int Y, int iterations, int size, int rank, MPI_Comm comm)
{
    int rect_width = X / size;
    int rect_height = Y / size;

    int start_x = rank % size * rect_width;
    int start_y = rank / size * rect_height;

    // initialize board 
    std::vector<std::vector<int>> global_board(X, std::vector<int>(Y));
    if (rank == 0) {
        global_board = random_board(X, Y);
        MPI_Bcast(&global_board, X*Y, MPI_Type_vector, 0, MPI_COMM_WORLD);
    }

    // Game Loop
    for (int i = 0; i < iterations; i++)
    {
        print_board(board);
        std::vector<std::vector<int>> new_board(X, std::vector<int>(Y, 0));
        std::vector<int> top_row(X);
        std::vector<int> bottom_row(X);
        std::vector<int> left_column(Y);
        std::vector<int> right_column(Y);
        for (int i = 0; i < X; i++){
            top_row[i] = global_board[i][start_y];
            bottom_row[i] = global_board[i][start_y + rect_height - 1];
        }
        for (int i = 0; i < Y; i++){
            left_column[i] = global_board[start_x][i];
            right_column[i] = global_board[start_x + rect_width - 1][i];
        }
        if (size == 8){
            MPI_Send(&top_row, X, MPI_INT, (rank + 4) % size, 0, MPI_COMM_WORLD);
            MPI_Send(&bottom_row, X, MPI_INT, (rank + 4) % size, 1, MPI_COMM_WORLD);
            if (rank == 0 || rank == 4) {
                MPI_Send(&left_column, Y, MPI_INT, rank + 3, 2, MPI_COMM_WORLD);
            } else {
                MPI_Send(&left_column, Y, MPI_INT, rank - 1, 2, MPI_COMM_WORLD);
            }
            if (rank == 3 || rank == 7) {
                MPI_Send(&right_column, Y, MPI_INT, rank - 3, 3, MPI_COMM_WORLD);
            } else {
                MPI_Send(&right_column, Y, MPI_INT, rank + 1, 3, MPI_COMM_WORLD);
            }
            MPI_Recv(&top_row, X, MPI_INT, (rank + 4) % size, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&bottom_row, X, MPI_INT, (rank + 4) % size, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (rank == 0 || rank == 4) {
                MPI_Recv(&left_column, Y, MPI_INT, rank + 3, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            } else {
                MPI_Recv(&left_column, Y, MPI_INT, rank - 1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if (rank == 3 || rank == 7) {
                MPI_Recv(&right_column, Y, MPI_INT, rank - 3, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            } else {
                MPI_Recv(&right_column, Y, MPI_INT, rank + 1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }

        for (int x = start_x; x < rect_width; x++)
        {
            for (int y = start_y; y < rect_height; y++)
            {
                new_board[x][y] = isAlive(board, x, y) ? 1 : 0;
            }
        }
        if (rank == 0){

        }
        board = new_board;
        usleep(100000);
    }


}