#include <vector>
#include <igraph.h>
#include <igraph_matrix.h>
#include <iostream>
using namespace std;

void igraph_matrix_view(igraph_matrix_t *A, igraph_real_t *data, int nrows, int ncols)
{
    A->ncol = nrows;
    A->nrow = ncols;
    A->data.stor_begin = data;
    A->data.stor_end = data+ncols*nrows;
    A->data.end = A->data.stor_end;
}

int main(int argc, char *argv[])
{
    vector<double> x(16);
    for (int i=0; i<16; i++)
        x[i]=i;

    igraph_matrix_t A;

    igraph_matrix_view(&A, x.data(),4,4);

    for (int i=0; i<4; i++)
    {
        for (int j=0; j<4; j++)
        {
            cout << MATRIX(A,i,j) << " ";
        }
        cout << endl;
    }
    return 0;
}
