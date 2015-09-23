#ifndef MODULARITYVERTEXPARTITION_H
#define MODULARITYVERTEXPARTITION_H

#include "MutableVertexPartition.h"

class ModularityVertexPartition : public MutableVertexPartition
{
public:
    ModularityVertexPartition(Graph* graph,
                              vector<size_t> membership);
    ModularityVertexPartition(Graph* graph);
    virtual ~ModularityVertexPartition();
    virtual ModularityVertexPartition* create(Graph* graph);

    virtual double diff_move(size_t v, size_t new_comm);
    virtual double quality();

protected:
private:
};

#endif // MODULARITYVERTEXPARTITION_H
