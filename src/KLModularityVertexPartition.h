#ifndef KLMODULARITYVERTEXPARTITION_H
#define KLMODULARITYVERTEXPARTITION_H

#include "MutableVertexPartition.h"

class KLModularityVertexPartition: public MutableVertexPartition
{
public:
    KLModularityVertexPartition(Graph* graph, vector<size_t> membership);
    KLModularityVertexPartition(Graph* graph, KLModularityVertexPartition* partition);
    KLModularityVertexPartition(Graph* graph);
    virtual ~KLModularityVertexPartition();
    virtual KLModularityVertexPartition* create(Graph* graph);

    virtual double diff_move(size_t v, size_t new_comm);
    virtual double quality();
protected:
private:
};

#endif // KLMODULARITYVERTEXPARTITION_H
