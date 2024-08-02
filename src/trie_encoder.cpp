#include "trie_encoder.h"

TrieRelation::TrieRelation(catalog* storage, Vertex src, Vertex dst){
    Vertex src_tmp = std::min(src, dst);
    Vertex dst_tmp = std::max(src, dst);
    src_ = src;
    dst_ = dst;
    edge_relation& target_edge_relation = storage->edge_relations_[src_tmp][dst_tmp];
    edge* edges = target_edge_relation.edges_;
    uint32_t edge_size = target_edge_relation.size_;
    uint32_t src_idx = 0, dst_idx = 1;
    size_ = edge_size;
    if(src > dst){
        std::sort(edges, edges + edge_size, [](edge& l, edge& r) -> bool {
            if (l.vertices_[1] == r.vertices_[1])
                return l.vertices_[0] < r.vertices_[0];
            return l.vertices_[1] < r.vertices_[1];
        });
        swap(src_idx, dst_idx);
    }

    // start encoding
    uint32_t start_offset=0;
    children_ = new Vertex [edge_size];
    Vertex current_vertex = edges[0].vertices_[src_idx];
    children_[0] = edges[0].vertices_[dst_idx];
    for(uint32_t i=1; i<edge_size; ++i){
        Vertex candidate_src = edges[i].vertices_[src_idx];
        Vertex candidate_dst = edges[i].vertices_[dst_idx];

        if(candidate_src != current_vertex){
            offset_map[current_vertex] = {start_offset+children_, i-start_offset};
            start_offset = i;
            current_vertex = candidate_src;
        }
        children_[i] = candidate_dst;
    }
    if(start_offset<edge_size){
        Vertex candidate_src = edges[start_offset].vertices_[src_idx];
        offset_map[candidate_src] = {start_offset+children_, edge_size-start_offset};
    }
}

void TrieRelation::remove_candidate(Vertex src_u, Vertex src_v, Graph* data_graph){
    if(src_u == src_){
        offset_map.erase(src_v);
    }else{
        uint32_t nbrs_count;
        const Vertex* nbrs = data_graph->getVertexNeighbors(src_v, nbrs_count);
        for(int i=0;i<nbrs_count;++i){
            Vertex v_n = nbrs[i];
            auto itf = offset_map.find(v_n);
            if(itf != offset_map.end()){
                Vertex* start = itf->second.first;
                const uint32_t size = itf->second.second;
                uint32_t valid_size = 0;
                for(int j=0;j<size;++j){
                    if(start[j] != src_v){
                        start[valid_size] = start[j];
                        valid_size ++;
                    }
                }
                itf->second.second = valid_size;
            }
        }
    }
}

TrieRelation::~TrieRelation(){
    delete children_;
}

TrieEncoder::TrieEncoder(catalog* storage, vector<Vertex>& order, vector<Vertex>& order_index, Graph* query_graph, Graph* data_graph, vector<vector<uint32_t>>& successor_depth){
    order_ = order;
    storage_ = storage;
    candidate_edges.resize(order.size());
    for(int i=0;i<order.size();++i){
        candidate_edges[i] = vector<TrieRelation*>(order.size(), NULL);
    }
    order_index_ = order_index;
    query_graph_ = query_graph;
    data_graph_ = data_graph;
    candidates_set_depth.resize(order.size()+1);
    first_successor_depth_.resize(order.size()+1);
    // start encoding
    for(int i=1;i<order_.size();++i){
        Vertex src = order_[i];

        uint32_t nbrs_count;
        const Vertex* nbrs = query_graph->getVertexNeighbors(src, nbrs_count);

        for(int j=0;j<nbrs_count;++j){
            Vertex dst = nbrs[j];
            uint32_t dst_depth = order_index[dst];
            if(order_index[src] < order_index[dst]){
                candidate_edges[i][dst_depth] = new TrieRelation(storage, src, dst);
            }
        }
        if(successor_depth[i].size() > 0){
            first_successor_depth_[i] = successor_depth[i][0];
        }else{
            first_successor_depth_[i] = 0;
        }
    }
}

void TrieEncoder::get_candidates(uint32_t u_depth, vector<Vertex>& result){
    for(int i=0;i<order_.size();++i){
        if(candidate_edges[u_depth][i] != NULL){
            for(auto& p : candidate_edges[u_depth][i]->offset_map){
                result.push_back(p.first);
            }
            sort(result.begin(), result.end());
            return;
        }
    }
}

unordered_map<Vertex, pair<Vertex*, uint32_t>>* TrieEncoder::get_candidates_map(uint32_t u_depth){
    int suc_depth = first_successor_depth_[u_depth];
    if(suc_depth != 0){
        unordered_map<Vertex, pair<Vertex*, uint32_t>>& offset_map = candidate_edges[u_depth][suc_depth]->offset_map;
        return &offset_map;
    }else{
        if(candidates_set_depth[u_depth].size() == 0){
            // construct
            uint32_t count;
            const Vertex* nbrs = query_graph_->getVertexNeighbors(order_[u_depth], count);
            int pred_depth = order_index_[nbrs[0]];
            Vertex* children = candidate_edges[pred_depth][u_depth]->children_;
            for(uint32_t i=0;i<candidate_edges[pred_depth][u_depth]->size_;++i){
                candidates_set_depth[u_depth].insert({children[i], {NULL, 0}});
            }
        }
        return &(candidates_set_depth[u_depth]);
    }
}

bool TrieEncoder::validate_neighbor_containment(uint32_t cur_depth, Vertex root1, Vertex root2){
    unordered_map<Vertex, pair<Vertex*, uint32_t>>* candidate_set = get_candidates_map(cur_depth);
    uint32_t count;
    const Vertex* nbrs = data_graph_->getVertexNeighbors(root2, count);
    Vertex v_n=0;
    for(int i=0;i<count;++i){
        v_n = nbrs[i];
        if(candidate_set->find(v_n) != candidate_set->end()){
            if(v_n != root1 && data_graph_->checkEdgeExistence(v_n, root1) == false){
                return false;
            }
        }
    }
    return true;
}

Vertex* TrieEncoder::get_edge_candidate(uint32_t src_depth, uint32_t dst_depth, Vertex src_v, uint32_t& count){
    auto& p = candidate_edges[src_depth][dst_depth]->offset_map[src_v];
    count = p.second;
    // if(count == 0){
    //     uint32_t nbrs_count;
    //     cout<<"remove:"<<src_depth<<":"<<src_v<<endl;
    //     const Vertex* nbrs = query_graph_->getVertexNeighbors(order_[src_depth], nbrs_count);
    //     for(int i=0;i<count;++i){
    //         Vertex u_n = nbrs[i];
    //         uint32_t depth_n = order_index_[u_n];
    //         if(depth_n>src_depth){
    //             candidate_edges[src_depth][depth_n]->remove_candidate(src_depth, src_v, data_graph_);
    //         }else{
    //             candidate_edges[depth_n][src_depth]->remove_candidate(depth_n, src_v, data_graph_);
    //         }
    //     }
    // }
    return p.first;
}

// void TrieEncoder::get_joint_candidate(uint32_t suc_depth, uint32_t dst_depth, Vertex src_v, vector<Vertex>& result){
//     uint32_t count;
//     const Vertex* v_neighbors = data_graph_->getVertexNeighbors(src_v, count);
//     if(suc_depth != 0){
//         unordered_map<Vertex, pair<Vertex*, uint32_t>>& offset_map = candidate_edges[dst_depth][suc_depth]->offset_map;
//         for(int i=0;i<count;++i){
//             Vertex v_n = v_neighbors[i];
//             if(offset_map.find(v_n) != offset_map.end()){
//                 result.push_back(v_n);
//             }
//         }
//     }else{
//         if(candidates_set_depth[dst_depth].size() == 0){
//             // construct
//             uint32_t count;
//             const Vertex* nbrs = query_graph_->getVertexNeighbors(order_[dst_depth], count);
//             int pred_depth = order_index_[nbrs[0]];
//             Vertex* children = candidate_edges[pred_depth][dst_depth]->children_;
//             for(uint32_t i=0;i<candidate_edges[pred_depth][dst_depth]->size_;++i){
//                 candidates_set_depth[dst_depth].insert(children[i]);
//             }
//         }
//         for(int i=0;i<count;++i){
//             Vertex v_n = v_neighbors[i];
//             if(candidates_set_depth[dst_depth].find(v_n) != candidates_set_depth[dst_depth].end()){
//                 result.push_back(v_n);
//             }
//         }
//     }
// }

TrieEncoder::~TrieEncoder(){
    for(int i=0;i<candidate_edges.size();++i){
        for(int j=0;j<candidate_edges[i].size();++j){
            if(candidate_edges[i][j] != NULL){
                delete candidate_edges[i][j];
            }
        }
    }
}
