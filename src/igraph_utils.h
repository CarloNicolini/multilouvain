#ifndef _IGRAPH_ADDITIONAL_UTILS_
#define _IGRAPH_ADDITIONAL_UTILS_

#include <igraph.h>
#include <igraph_error.h>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <Eigen/Core>

using namespace Eigen;

#define IGRAPH_TRY(call){\
    int __result = call;\
    std::stringstream ss; ss << __result; \
    if (__result != 0)\
{\
    throw std::runtime_error("Igraph Error " + ss.str());\
    } \
    }

#define IGRAPHPP_TRY_NEW(variable, type)   \
    try {                                \
    variable = new type;             \
    } catch (const std::bad_alloc&) { \
    IGRAPH_ERROR("std::bad_alloc thrown in C++ code", IGRAPH_ENOMEM); \
    }

void igraph_matrix_view(igraph_matrix_t *A, igraph_real_t *data, int nrows, int ncols);

template<typename M>
M read_adj_matrix(const std::string & path, const char sep=' ')
{
    std::ifstream indata;
    indata.open(path);
    std::string line;
    std::vector<double> values;
    uint rows = 0;
    while (std::getline(indata, line))
    {
        std::stringstream lineStream(line);
        std::string cell;
        while (std::getline(lineStream, cell, sep))
        {
            values.push_back(std::stod(cell));
        }
        ++rows;
    }
    return Eigen::Map<const Eigen::Matrix<typename M::Scalar, M::RowsAtCompileTime, M::ColsAtCompileTime, RowMajor>>(values.data(), rows, values.size()/rows);
}


// Vector pretty printer
template<typename T>
std::ostream & operator<<(std::ostream & os, std::vector<T> vec)
{
    os<<"[ ";
    std::copy(vec.begin(), vec.end(), std::ostream_iterator<T>(os, " "));
    os<<"]";
    return os;
}


namespace Eigen{
template<class Matrix>
void write_binary(const char* filename, const Matrix& matrix){
    std::ofstream out(filename,std::ios::out | std::ios::binary | std::ios::trunc);
    typename Matrix::Index rows=matrix.rows(), cols=matrix.cols();
    out.write((char*) (&rows), sizeof(typename Matrix::Index));
    out.write((char*) (&cols), sizeof(typename Matrix::Index));
    out.write((char*) matrix.data(), rows*cols*sizeof(typename Matrix::Scalar) );
    out.close();
}
template<class Matrix>
void read_binary(const char* filename, Matrix& matrix){
    std::ifstream in(filename,std::ios::in | std::ios::binary);
    typename Matrix::Index rows=0, cols=0;
    in.read((char*) (&rows),sizeof(typename Matrix::Index));
    in.read((char*) (&cols),sizeof(typename Matrix::Index));
    matrix.resize(rows, cols);
    in.read( (char *) matrix.data() , rows*cols*sizeof(typename Matrix::Scalar) );
    in.close();
}
}

#endif
