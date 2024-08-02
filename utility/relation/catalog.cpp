#include "catalog.h"
#include <vector>

void catalog::initialize_catalog_info() {
     for (uint32_t u = 0; u < num_sets_; ++u) {
        for (uint32_t v = u + 1; v < num_sets_; ++v) {
            if (query_graph_->checkEdgeExistence(u, v)) {
                RelationMetaInfo meta_info;
                if (query_graph_->getCoreValue(u) > 1 && query_graph_->getCoreValue(v) > 1) {
                    meta_info.type = EdgeType::CoreEdge;
                }
                else if (query_graph_->getVertexDegree(u) == 1 || query_graph_->getVertexDegree(v) == 1) {
                    meta_info.type = EdgeType::LeafEdge;
                }
                else {
                    meta_info.type = EdgeType::TreeEdge;
                }
                catalog_info_.insert(std::make_pair(std::make_pair(u, v), meta_info));
            }
        }
    }
}
