#ifndef _IGRAPH_ADDITIONAL_UTILS_
#define _IGRAPH_ADDITIONAL_UTILS_
#include <igraph.h>
#include <igraph_error.h>
#include <stdexcept>
#include <sstream>

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

#endif
