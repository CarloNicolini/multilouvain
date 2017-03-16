#include "igraph_utils.h"

void igraph_matrix_view(igraph_matrix_t *A, igraph_real_t *data, int nrows, int ncols)
{
    A->ncol = nrows;
    A->nrow = ncols;
    A->data.stor_begin = data;
    A->data.stor_end = data+ncols*nrows;
    A->data.end = A->data.stor_end;
}


//Eigen::MatrixXd read_adj_matrix(const std::string &filename)
//{
//    Eigen::MatrixXd adj_mat;

//    std::ifstream file;
//    file.open(filename.c_str(),std::ios::in);

//    std::vector<double> data;
//    std::string line;

//    unsigned int nRows=1;
//    unsigned int nCols=1;
//    while(!std::getline(file, line, '\n').eof())
//    {
//        std::istringstream reader(line);
//        std::vector<double> lineData;
//        nCols=0;
//        while(!reader.eof())
//        {
//            double val;
//            reader >> val;
//            if(reader.fail())
//                break;
//            lineData.push_back(val);
//            nCols++;
//        }
//        data.push_back(lineData);
//        nRows++;
//    }
//    nRows--;
//    // Deep copy of data array to the adjacency matrix
//    adj_mat.setZero(nRows,nCols);
//    for (unsigned int i=0; i<nRows; ++i)
//        for (unsigned int j=0; j<nCols; ++j)
//            adj_mat.coeffRef(i,j)=data[i][j];

//    return adj_mat;
//}

