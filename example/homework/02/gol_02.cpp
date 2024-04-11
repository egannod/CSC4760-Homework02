// Sequential Game of Life for Sanity's Sake
#include <iostream>
#include <vector>
#include <unistd.h>


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

bool isAlive(std::vector<std::vector<int>>& board, const int x, const int y)
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
    X = atoi(argv[1]);
    Y = atoi(argv[2]);
    iterations = atoi(argv[3]);

    std::cout << "X: " << X << " Y: " << Y << " iterations: " << iterations << std::endl;


    // initialize board 
    std::vector<std::vector<int>> board(X, std::vector<int>(Y));
    board = random_board(X, Y);

    // Game Loop
    for (int i = 0; i < iterations; i++)
    {
        print_board(board);
        std::vector<std::vector<int>> new_board(X, std::vector<int>(Y, 0));
        for (int x = 0; x < X; x++)
        {
            for (int y = 0; y < Y; y++)
            {
                new_board[x][y] = isAlive(board, x, y) ? 1 : 0;
            }
        }
        board = new_board;
        usleep(100000);
    }

    return 0;
}
