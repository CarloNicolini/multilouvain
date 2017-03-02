#ifndef _IGRAPH_ADDITIONAL_UTILS_
#define _IGRAPH_ADDITIONAL_UTILS_
#include <igraph.h>
#include <igraph_error.h>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include <vector>
#include <Eigen/Core>

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

std::vector<double> read_adj_matrix(const std::string &filename);


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
} // Eigen::

#endif
