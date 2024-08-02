#include "containment_filter.h"

ContainmentFilter::ContainmentFilter(Graph* query_graph, Graph* data_graph, vector<vector<uint32_t>>& successors_depth, vector<vector<uint32_t>>& predecessor_depth, vector<Vertex>& search_order, vector<Vertex>& order_index, vector<bitset<MAX_QUERY_SIZE>>& ancestors_depth, TrieEncoder* encoder){
    query_graph_ = query_graph;
    data_graph_ = data_graph;
    order_ = search_order;
    encoder_ = encoder;

    query_vertex_count_ = query_graph_->getVerticesCount();
    
    history_.resize(query_vertex_count_+2);
    history_root_vertices_.resize(query_vertex_count_+2);
    empty_record_offset_.resize(query_vertex_count_+2);
    conflict_record_offset_.resize(query_vertex_count_+2);
    history_conflict_source_.resize(query_vertex_count_+2);

    successors_with_same_labels_.resize(query_vertex_count_+2);
    history_max_search_depth.resize(query_vertex_count_+2);
    ancestors_depth_ = ancestors_depth;
    
    failing_set_recorder_.resize(query_vertex_count_+2);
    successor_offset_map_.resize(query_vertex_count_+2);

    successors_depth_ = successors_depth;
    predecessor_depth_ = predecessor_depth;
    none_succesors_depth_.resize(query_vertex_count_+1);
    unmatched_query_depth_dist_.resize(query_vertex_count_+1);
    bitset<MAX_QUERY_SIZE> touch_query_depth;
    for(uint32_t d=1; d <= query_vertex_count_; ++d){
        none_succesors_depth_[d].clear();
        touch_query_depth.set(d);
        successors_with_same_labels_[d] = 0;
        for(uint32_t d_n = d+1; d_n<=query_vertex_count_; d_n++){
            bool is_successors = false;
            for(auto suc_depth : successors_depth_[d]){
                if(d_n == suc_depth){
                    is_successors = true;
                    break;
                }
            }
            if(is_successors == false){
                none_succesors_depth_[d].push_back(d_n);
                if(touch_query_depth.test(d_n) == true){
                    unmatched_query_depth_dist_[d].second.push_back(d_n);
                }
            }else{
                unmatched_query_depth_dist_[d].first.push_back(d_n);
                if(query_graph_->getVertexLabel(order_[d]) == query_graph_->getVertexLabel(order_[d_n])){
                    successors_with_same_labels_[d] ++;
                }
            }
        }
        successor_offset_map_[d].resize(query_vertex_count_+1);
        int offset = 0;
        for(auto suc_depth : successors_depth_[d]){
            touch_query_depth.set(suc_depth);
            successor_offset_map_[d][suc_depth] = offset;
            offset ++;
        }
    }

    mask_.resize(query_vertex_count_+1);
    for(int x=0;x<=query_vertex_count_;++x){
        mask_[x].reset();
        for(int y=x;y<=query_vertex_count_;++y){
            mask_[x].set(y);
        }
    }

#if DEBUG==1
    same_label_depth_.resize(query_vertex_count_+1);
    for(uint32_t d=1; d <= query_vertex_count_; ++d){
        for(uint32_t d_n = d+1; d_n<=query_vertex_count_; d_n++){
            if(query_graph->getVertexLabel(order_[d]) == query_graph->getVertexLabel(order_[d_n])){
                if(query_graph->checkEdgeExistence(order_[d], order_[d_n])==true){
                    same_label_depth_[d].first.push_back(d_n);
                }else{
                    same_label_depth_[d].second.push_back(d_n);
                }
            }
        }
    }
    preprocessed_depth_intersection.resize(query_vertex_count_+1);
    for(uint32_t d=2; d<=query_vertex_count_;++d){
        preprocessed_depth_intersection[d].resize(query_vertex_count_+1);
        for(uint32_t d_c = d-1; d_c>0; --d_c){
            for(uint32_t d_n=1;d_n<=query_vertex_count_;++d_n){
                if(query_graph->checkEdgeExistence(order_[d_n], order_[d_c]) == false && query_graph->checkEdgeExistence(order_[d_n], order_[d]) == true){
                    if(d_n < d_c){
                        preprocessed_depth_intersection[d][d_c].first.push_back(d_n);
                    }else if(d_n > d_c){
                        preprocessed_depth_intersection[d][d_c].second.push_back(d_n);
                    }
                }
            }
        }
    }
#endif
}

bool ContainmentFilter::validate_containment_for_homo(uint32_t depth, vector<vector<pair<Vertex*, uint32_t>>>& candidates_stack, Vertex root, CandidatesHistoryStack& recorder){
    for(auto offset : recorder.validate_completion_offset_){
        bool contained = true;
        
        Vertex history_root = recorder.history_root_vertices_[offset];
        int suc_offset = 0;
        vector<pair<Vertex*, uint32_t>>& history = recorder.history_candidates_[offset];

        for(auto suc_depth : successors_depth_[depth]){
            pair<Vertex*, uint32_t>& p = *(candidates_stack[suc_depth].rbegin());
            if(validate_vector_containment_strict(history[suc_offset].first, history[suc_offset].second, p.first, p.second) == false){ // validate_vector_containment(history[depth][offset][suc_offset], *(candidates_stack[suc_depth].rbegin())) == false
                contained = false;
                break;
            }
            suc_offset ++;
        }

        if(contained == true){
            return true;
        }
    }
    return false;
}

bool ContainmentFilter::validate_exclusion_containment(uint32_t depth, vector<vector<pair<Vertex*, uint32_t>>>& candidates_stack, Vertex root, CandidatesHistoryStack& recorder){
    // for(int offset=0;offset<recorder.history_root_vertices_.size();++offset){ // int i=0;i<failing_set_recorder_[depth].size();++i
    //     if(recorder.history_conflict_sources_[offset].any() == true || recorder.validation_marker_[offset] == 0){
    //         continue;
    //     }
    for(auto offset : recorder.validate_exclusion_offsets_){
        
        bool contained = true;
        
        Vertex history_root = recorder.history_root_vertices_[offset];
        int suc_offset = 0;
        vector<pair<Vertex*, uint32_t>>& history = recorder.history_candidates_[offset];

        for(auto suc_depth : successors_depth_[depth]){
            pair<Vertex*, uint32_t>& p = *(candidates_stack[suc_depth].rbegin());
            if(validate_vector_containment_strict(history[suc_offset].first, history[suc_offset].second, p.first, p.second) == false){ // validate_vector_containment(history[depth][offset][suc_offset], *(candidates_stack[suc_depth].rbegin())) == false
                contained = false;
                break;
            }
            suc_offset ++;
        }

        if(contained == true){
            return true;
        }
    }
    return false;
}


#if DEBUG==1
bool ContainmentFilter::validate_complete_containment(uint32_t depth, vector<vector<pair<Vertex*, uint32_t>>>& candidates_stack, Vertex root, CandidatesHistoryStack& recorder, vector<Vertex>& embedding_depth){
    for(auto offset : recorder.validate_completion_offset_){
        bool contained = true;
        Vertex history_root = recorder.history_root_vertices_[offset];
        vector<pair<Vertex*, uint32_t>>& history = recorder.history_candidates_[offset];
        // scan constraints
        int suc_offset = 0;
        for(auto d : unmatched_query_depth_dist_[depth].first){
            pair<Vertex*, uint32_t>& p = *(candidates_stack[d].rbegin());
            if(validate_vector_containment_relaxed(history[suc_offset].first, history[suc_offset].second, p.first, p.second, history_root, root) == false){
                contained = false;
                break;
            }
            suc_offset ++;
        }

#if COMPLETE_CONSTRAINT_OVERLAP == 0
        if(contained == true){
            suc_offset = 0;
            bitset<MAX_QUERY_SIZE> fs = recorder.history_failing_set_[offset];
            bitset<MAX_QUERY_SIZE> explored_depth = fs & mask_[depth];
            int max_depth = order_.size()-1;
            for(;max_depth>=depth;max_depth--){
                if(fs.test(max_depth) == true){
                    break;
                }
            }
            for(auto d : same_label_depth_[depth].first){
                if(d>max_depth){
                    continue;
                }
                pair<Vertex*, uint32_t>& p = *(candidates_stack[d].rbegin());
                if(validate_vector_element_existence(p.first, p.second, history_root) == true){
                    if(validate_vector_element_existence(history[suc_offset].first, history[suc_offset].second, root) == true){
                        for(auto d_v : preprocessed_depth_intersection[d][depth].second){
                            if(candidates_stack[d_v].empty()){
                                // candidates_containment_check.insert(d_v);
                                if((explored_depth & ancestors_depth_[d_v]).any() == false){
                                    continue;
                                }
                                if(encoder_->validate_neighbor_containment(d_v, root, history_root) == false){
                                    contained = false;
                                    break;
                                }
                            }else{
                                // containment_check.insert(d_v);
                                if((explored_depth & ancestors_depth_[d_v]).any() == false){
                                    continue;
                                }
                                pair<Vertex*, uint32_t>& p = *(candidates_stack[d_v].rbegin());
                                Vertex* cans = p.first;
                                Vertex v;
                                for(uint32_t x=0;x<p.second; ++x){
                                    v = cans[x];
                                    if(data_graph_->checkEdgeExistence(v, history_root) == true && v != root){
                                        if(data_graph_->checkEdgeExistence(v, root) == false){
                                            contained = false;
                                            break;
                                        }
                                    }
                                }
                                // if(contained == false){
                                //     break;
                                // }
                            }
                        }
                    }else{
                        contained = false;
                        break;
                    }
                }
                suc_offset++;
            }
            if(contained == false){
                continue;
            }
            for(auto d : same_label_depth_[depth].second){
                if(d>max_depth){
                    continue;
                }
                if(fs.test(d) == false){
                    continue;
                }
                pair<Vertex*, uint32_t>& p = *(candidates_stack[d].rbegin());
                if(candidates_stack[d].empty() == true){
                    for(auto d_v : preprocessed_depth_intersection[d][depth].first){
                        // connection_check.insert(embedding_depth[d_v]);
                        if(data_graph_->checkEdgeExistence(root, embedding_depth[d_v]) == false){
                            contained = false;
                            break;
                        }
                    }
                    for(auto d_v : preprocessed_depth_intersection[d][depth].second){
                        if(candidates_stack[d_v].empty()){
                            // candidates_containment_check.insert(d_v);
                            if((explored_depth & ancestors_depth_[d_v]).any() == false){
                                continue;
                            }
                            if(encoder_->validate_neighbor_containment(d_v, root, history_root) == false){
                                contained = false;
                                break;
                            }
                        }else{
                            // containment_check.insert(d_v);
                            if((explored_depth & ancestors_depth_[d_v]).any() == false){
                                continue;
                            }
                            pair<Vertex*, uint32_t>& p = *(candidates_stack[d_v].rbegin());
                            Vertex* cans = p.first;
                            Vertex v;
                            for(uint32_t x=0;x<p.second; ++x){
                                v = cans[x];
                                if(data_graph_->checkEdgeExistence(v, history_root) == true && v != root){
                                    if(data_graph_->checkEdgeExistence(v, root) == false){
                                        contained = false;
                                        break;
                                    }
                                }
                            }
                            // if(contained == false){
                            //     break;
                            // }
                        }
                    }
                }else if(validate_vector_element_existence(p.first, p.second, history_root) == true){
                    if(validate_vector_element_existence(p.first, p.second, root) == true){
                        for(auto d_v : preprocessed_depth_intersection[d][depth].second){
                            if(candidates_stack[d_v].empty()){
                                // candidates_containment_check.insert(d_v);
                                if((explored_depth & ancestors_depth_[d_v]).any() == false){
                                    continue;
                                }
                                if(encoder_->validate_neighbor_containment(d_v, root, history_root) == false){
                                    contained = false;
                                    break;
                                }
                            }else{
                                // containment_check.insert(d_v);
                                if((explored_depth & ancestors_depth_[d_v]).any() == false){
                                    continue;
                                }
                                pair<Vertex*, uint32_t>& p = *(candidates_stack[d_v].rbegin());
                                Vertex* cans = p.first;
                                Vertex v;
                                for(uint32_t x=0;x<p.second; ++x){
                                    v = cans[x];
                                    if(data_graph_->checkEdgeExistence(v, history_root) == true && v != root){
                                        if(data_graph_->checkEdgeExistence(v, root) == false){
                                            contained = false;
                                            break;
                                        }
                                    }
                                }
                                // if(contained == false){
                                //     break;
                                // }
                            }
                        }
                    }else{
                        contained = false;
                        break;
                    }
                }
            }
            if(contained == true){
                return true;
            }

        }


#else
        if(contained == true){
            // unordered_set<Label> neighbor_check_label;
            unordered_set<uint32_t> candidates_containment_check;
            unordered_set<Vertex> connection_check;
            unordered_set<int> containment_check;
            // suc_offset = 0;
            bitset<MAX_QUERY_SIZE> fs = recorder.history_failing_set_[offset];
            bitset<MAX_QUERY_SIZE> explored_depth = fs & mask_[depth];
            int max_depth = order_.size()-1;
            for(;max_depth>=depth;max_depth--){
                if(fs.test(max_depth) == true){
                    break;
                }
            }
            vector<uint32_t>& suc_offset_map = successor_offset_map_[depth];
            for(auto d : same_label_depth_[depth].first){
                if(d>max_depth){
                    continue;
                }
                pair<Vertex*, uint32_t>& p = *(candidates_stack[d].rbegin());
                if(validate_vector_element_existence(p.first, p.second, history_root) == true){
                    // validate_vector_element_existence(history[suc_offset].first, history[suc_offset].second, root)
                    if(validate_vector_element_existence(history[suc_offset_map[d]].first, history[suc_offset_map[d]].second, root) == true){
                        // for(auto d_v : preprocessed_depth_intersection[d][depth].first){
                        //     connection_check.insert(embedding_depth[d_v]);
                        // }
                        for(auto d_v : preprocessed_depth_intersection[d][depth].second){
                            // if(d_v>max_depth){
                            //     continue;
                            // }
                            if(candidates_stack[d_v].empty()){
                                // neighbor_check_label.insert(query_graph_->getVertexLabel(order_[d_v]));
                                candidates_containment_check.insert(d_v);
                            }else{
                                containment_check.insert(d_v);
                            }
                        }
                    }else{
                        contained = false;
                        break;
                    }
                }
                // suc_offset++;
            }
            if(contained == false){
                continue;
            }
            for(auto d : same_label_depth_[depth].second){
                if(d>max_depth){
                    continue;
                }
                if(fs.test(d) == false){
                    continue;
                }
                pair<Vertex*, uint32_t>& p = *(candidates_stack[d].rbegin());
                if(candidates_stack[d].empty() == true){
                    for(auto d_v : preprocessed_depth_intersection[d][depth].first){
                        connection_check.insert(embedding_depth[d_v]);
                    }
                    for(auto d_v : preprocessed_depth_intersection[d][depth].second){
                        // if(d_v>max_depth){
                        //     continue;
                        // }
                        if(candidates_stack[d_v].empty()){
                            // neighbor_check_label.insert(query_graph_->getVertexLabel(order_[d_v]));
                            candidates_containment_check.insert(d_v);
                        }else{
                            containment_check.insert(d_v);
                        }
                    }
                }else if(validate_vector_element_existence(p.first, p.second, history_root) == true){
                    if(validate_vector_element_existence(p.first, p.second, root) == true){
                        // for(auto d_v : preprocessed_depth_intersection[d][depth].first){
                        //     connection_check.insert(embedding_depth[d_v]);
                        // }
                        for(auto d_v : preprocessed_depth_intersection[d][depth].second){
                            // if(d_v>max_depth){
                            //     continue;
                            // }
                            if(candidates_stack[d_v].empty()){
                                // neighbor_check_label.insert(query_graph_->getVertexLabel(order_[d_v]));
                                candidates_containment_check.insert(d_v);
                            }else{
                                containment_check.insert(d_v);
                            }
                        }
                    }else{
                        contained = false;
                        break;
                    }
                }
            }
            
            // start validation
            if(contained == true){
                for(auto d : candidates_containment_check){
                    // if((explored_depth & ancestors_depth_[d]).any() == false){
                    //     continue;
                    // }

                    if(encoder_->validate_neighbor_containment(d, root, history_root) == false){
                        contained = false;
                        break;
                    }
                }
            }
            if(contained == true){
                for(auto& p : connection_check){
                    if(data_graph_->checkEdgeExistence(root, p) == false){
                        contained = false;
                        break;
                    }
                }
            }
            if(contained == true){
                for(auto& d : containment_check){
                    // if((explored_depth & ancestors_depth_[d]).any() == false){
                    //     continue;
                    // }
                    pair<Vertex*, uint32_t>& p = *(candidates_stack[d].rbegin());
                    Vertex* cans = p.first;
                    Vertex v;
                    for(uint32_t x=0;x<p.second; ++x){
                        v = cans[x];
#if HOMOMORPHISM == 1
                        if(data_graph_->checkEdgeExistence(v, history_root) == true){
                            if(v == root){
                                if(validate_vector_element_existence(p.first, p.second, history_root) == false){
                                    contained = false;
                                    break;
                                }
                                // for(auto pred_depth : predecessor_depth_[d]){
                                //     if(d > depth){
                                //         contained = false;
                                //         break;
                                //     }
                                // }

                            }
                            if(data_graph_->checkEdgeExistence(v, root) == false){
                                contained = false;
                                break;
                            }
                        }
#else
                        if(data_graph_->checkEdgeExistence(v, history_root) == true && v != root){
                            if(data_graph_->checkEdgeExistence(v, root) == false){
                                contained = false;
                                break;
                            }
                        }
#endif
                    }
                    if(contained == false){
                        break;
                    }
                }
            }
            if(contained == true){
                return true;
            }
        }
#endif
    }
    return false;
}
#endif


bool ContainmentFilter::validate_vector_containment_strict(Vertex* vec1, int32_t size1, Vertex* vec2, int32_t size2){
    if(size1<size2){
        return false;
    }
    int i=0,j=0;
    for(; i<size1 && j<size2;){
        if(vec1[i]==vec2[j]){
            ++i;++j;
        }else if(vec1[i]>vec2[j]){
            ++j;
            return false;
        }else{
            ++i;
        }
    }
    if(j<size2){
        return false;
    }
    return true;
}


bool ContainmentFilter::validate_vector_containment_relaxed(const Vertex* vec1, int32_t vec1_size, const Vertex* vec2, int32_t vec2_size, Vertex r1, Vertex r2){
    if(vec1_size<vec2_size-1){
        return false;
    }
    int i=0,j=0;
    for(; i<vec1_size && j<vec2_size;){
        if(vec1[i] == r2){
            ++i;
            if(i>=vec1_size){
                break;
            }
        }
        if(vec2[j] == r1){
            ++j;
            if(j>=vec2_size){
                break;
            }
        }
        if(vec1[i]==vec2[j]){
            ++i;++j;
        }else if(vec1[i]>vec2[j]){
            ++j;
            return false;
        }else{
            ++i;
        }
    }
    
    if(j<vec2_size){
        if(vec2[j] == r1){
            ++j;
            if(vec2_size == j){
                return true;
            }
        }
        return false;
    }
    return true;
}

bool ContainmentFilter::validate_vector_element_existence(Vertex* vec1, int32_t size1, Vertex v){
    for(int i=0;i<size1;++i){
        if(vec1[i] == v){
            return true;
        }else if(vec1[i] > v){
            return false;
        }
    }
    return false;
}

bool ContainmentFilter::validate_vector_containment_relaxed_full(Vertex* vec1, int32_t size1, Vertex* vec2, int32_t size2, Vertex r1, Vertex r2){
    if(size1<size2-1){
        return false;
    }
    int i=0,j=0;
    for(; i<size1 && j<size2;){
        if(vec1[i] == r2 || vec1[i] == r1){
            ++i;
            if(i>=size1){
                break;
            }
        }
        if(vec2[j] == r1 || vec2[j] == r2){
            ++j;
            if(j>=size2){
                break;
            }
        }
        if(vec1[i]==vec2[j]){
            ++i;++j;
        }else if(vec1[i]>vec2[j]){
            ++j;
            return false;
        }else{
            ++i;
        }
    }
    
    if(j<size2){
        if(vec2[j] == r1 || vec2[j] == r2){
            ++j;
            if(size2 == j){
                return true;
            }
        }
        return false;
    }
    return true;
}