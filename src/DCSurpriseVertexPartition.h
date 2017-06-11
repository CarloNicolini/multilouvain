#ifndef DCSURPRISEVERTEXPARTITION_H_
#define DCSURPRISEVERTEXPARTITION_H_

#include "MutableVertexPartition.h"

class DCSurpriseVertexPartition: public MutableVertexPartition
{
public:
    DCSurpriseVertexPartition(Graph* graph, vector<size_t> membership);
    DCSurpriseVertexPartition(Graph* graph, DCSurpriseVertexPartition* partition);
    DCSurpriseVertexPartition(Graph* graph);
    virtual ~DCSurpriseVertexPartition();
    virtual DCSurpriseVertexPartition* create(Graph* graph);

    virtual double diff_move(size_t v, size_t new_comm);
    virtual double quality();
protected:
private:
};

#endif // DCSURPRISEVERTEXPARTITION_H_
