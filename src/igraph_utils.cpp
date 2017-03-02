#include "igraph_utils.h"

void igraph_matrix_view(igraph_matrix_t *A, igraph_real_t *data, int nrows, int ncols)
{
    A->ncol = nrows;
    A->nrow = ncols;
    A->data.stor_begin = data;
    A->data.stor_end = data+ncols*nrows;
    A->data.end = A->data.stor_end;
}

std::vector<double> read_adj_matrix(const std::string &filename)
{
    std::vector<double> adj_mat;
    using std::ifstream;
    ifstream file;
    file.open(filename.c_str(),std::ios::in);
    std::string line;

    while(!std::getline(file, line, '\n').eof())
    {
        std::istringstream reader(line);
        std::vector<double> lineData;
        while(!reader.eof())
        {
            double val;
            reader >> val;
            if(reader.fail())
                break;
            adj_mat.push_back(val);
        }
    }

    return adj_mat;
}
