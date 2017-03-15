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

// Vector pretty printer
template<typename T>
std::ostream & operator<<(std::ostream & os, std::vector<T> vec)
{
    os<<"[ ";
    std::copy(vec.begin(), vec.end(), std::ostream_iterator<T>(os, " "));
    os<<"]";
    return os;
}

#endif
