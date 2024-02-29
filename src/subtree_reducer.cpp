#include "subtree_reducer.h"


// ======================================================================================
SubtreeReducer::SubtreeReducer(Graph* query_graph, Graph* data_graph, vector<vector<uint32_t>>& successors_depth, vector<vector<uint32_t>>& predecessor_depth, vector<Vertex>& search_order, vector<Vertex>& order_index, vector<bitset<MAX_QUERY_SIZE>>& ancestors_depth, TrieEncoder* encoder){
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
        for(auto suc_depth : successors_depth_[d]){
            touch_query_depth.set(suc_depth);
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

void SubtreeReducer::clear_record(uint32_t depth){
    for(auto& p : failing_set_recorder_[depth]){
        delete p.second;
    }
    failing_set_recorder_[depth].clear();
    history_root_vertices_[depth].clear();
    empty_record_offset_[depth].clear();
    conflict_record_offset_[depth].clear();
    history_[depth].clear();
    history_max_search_depth[depth].clear();
    history_conflict_source_[depth].clear();
}

bool SubtreeReducer::validate_standard_containment(uint32_t depth, vector<vector<vector<Vertex>>>& candidates_stack, Vertex root, FailingSet*& result_fs, bool full_descendent){
    result_fs = NULL;
    for(auto offset : empty_record_offset_[depth]){ // int i=0;i<failing_set_recorder_[depth].size();++i
        bool contained = true;
        
        Vertex history_root = history_root_vertices_[depth][offset];
        int suc_offset = 0;
        // if(full_descendent == false){
            for(auto suc_depth : successors_depth_[depth]){
                if(validate_vector_containment(history_[depth][offset][suc_offset], *(candidates_stack[suc_depth].rbegin())) == false){ // validate_vector_containment(history[depth][offset][suc_offset], *(candidates_stack[suc_depth].rbegin())) == false
                    contained = false;
                    break;
                }
                suc_offset ++;
            }
        // }else{
        //     for(auto suc_depth : successors_depth_[depth]){
        //         if(validate_vector_containment(history_[depth][offset][suc_depth-depth-1], *(candidates_stack[suc_depth].rbegin())) == false){ // validate_vector_containment(history[depth][offset][suc_depth-depth-1], *(candidates_stack[suc_depth].rbegin())) == false
        //             contained = false;
        //             break;
        //         }
        //         suc_offset ++;
        //     }
        // }
        
        if(contained == true){
            result_fs = failing_set_recorder_[depth][history_root];
            return true;
        }
    }
    return false;
}

bool SubtreeReducer::validate_standard_containment_test(uint32_t depth, vector<vector<vector<Vertex>>>& candidates_stack, Vertex root, FailingSet*& result_fs, bool full_descendent, Vertex& history_root_test){
    result_fs = NULL;
    // for(auto offset : empty_record_offset_[depth]){ // int i=0;i<failing_set_recorder_[depth].size();++i
    //     bool contained = true;
        
    //     Vertex history_root = history_root_vertices_[depth][offset];
    //     int suc_offset = 0;
    //     // if(full_descendent == false){
    //         for(auto suc_depth : successors_depth_[depth]){
    //             if(validate_vector_containment(history_[depth][offset][suc_offset], *(candidates_stack[suc_depth].rbegin())) == false){ // validate_vector_containment(history[depth][offset][suc_offset], *(candidates_stack[suc_depth].rbegin())) == false
    //                 contained = false;
    //                 break;
    //             }
    //             suc_offset ++;
    //         }
    //     if(contained == true){
    //         result_fs = failing_set_recorder_[depth][history_root];
    //         return false;
    //     }
    // }
    int contained_count = 0;
    bitset<MAX_QUERY_SIZE> conflict_depth;
    for(int offset=0;offset<failing_set_recorder_[depth].size();++offset){
        bool contained = true;
        Vertex history_root = history_root_vertices_[depth][offset];

        // scan constraints
        int suc_offset = 0;
        for(auto suc_depth : successors_depth_[depth]){
            if(validate_vector_containment_relaxed(history_[depth][offset][suc_offset], *(candidates_stack[suc_depth].rbegin()), history_root, root) == false){
            // if(validate_vector_containment(history_[depth][offset][suc_offset], *(candidates_stack[suc_depth].rbegin())) == false){
                contained = false;
                break;
            }
            suc_offset ++;
        }
        if(contained == true){
            contained_count ++;
            // conflict_depth |= history_conflict_source_[depth][offset];
            // if(contained_count > conflict_depth.count()+successors_with_same_labels_[depth] && contained_count > 1){
            //     history_root_test = history_root;
            //     return true;
            // }
            // result_fs = failing_set_recorder_[depth][history_root];
            if(contained_count > 1){
                return true;
            }
        }
    }
    return false;
}


#if DEBUG==1
bool SubtreeReducer::validate_extended_containment_debug(uint32_t depth, vector<vector<vector<Vertex>>>& candidates_stack, Vertex root, FailingSet*& result_fs, bool full_descendent, vector<Vertex>& embedding_depth){
    for(int offset=0;offset<failing_set_recorder_[depth].size();++offset){ // auto offset : conflict_record_offset[depth]
        bool contained = true;
        Vertex history_root = history_root_vertices_[depth][offset];
        uint32_t max_search_depth = history_max_search_depth[depth][offset];
        // scan constraints
        int suc_offset = 0;
        for(auto d : unmatched_query_depth_dist_[depth].first){
            if(validate_vector_containment_relaxed(history_[depth][offset][suc_offset], *(candidates_stack[d].rbegin()), history_root, root) == false){
                contained = false;
                break;
            }
            suc_offset ++;
        }

        if(contained == true){
            // unordered_set<Label> neighbor_check_label;
            unordered_set<uint32_t> candidates_containment_check;
            unordered_set<Vertex> connection_check;
            unordered_set<int> containment_check;
            suc_offset = 0;
            FailingSet* fs = failing_set_recorder_[depth][history_root];
            bitset<MAX_QUERY_SIZE> explored_depth = fs->failing_set_ & mask_[depth];
            int max_depth = order_.size()-1;
            for(;max_depth>=depth;max_depth--){
                if(fs->failing_set_.test(max_depth) == true){
                    break;
                }
            }

            for(auto d : same_label_depth_[depth].first){
                if(d>max_depth){
                    continue;
                }
                if(validate_vector_element_existence(*(candidates_stack[d].rbegin()), history_root) == true){
                    if(validate_vector_element_existence(history_[depth][offset][suc_offset], root) == true){
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
                suc_offset++;
            }
            if(contained == false){
                continue;
            }
            for(auto d : same_label_depth_[depth].second){
                if(d>max_depth){
                    continue;
                }
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
                }else if(validate_vector_element_existence(*(candidates_stack[d].rbegin()), history_root) == true){
                    if(validate_vector_element_existence(*(candidates_stack[d].rbegin()), root) == true){
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
                // if(neighbor_check_label.size() > 0){
                //     cout<<"require_checking:"<<depth<<":"<<root<<endl;
                // }
                // for(auto d : neighbor_check_label){
                for(auto d : candidates_containment_check){
                    if((explored_depth & ancestors_depth_[d]).any() == false){
                        continue;
                    }
                    vector<Vertex> v1_neighbor, v2_neighbor;

                    if(successors_depth_[d].size() > 0){
                        encoder_->get_joint_candidate(successors_depth_[d][0], d, history_root, v1_neighbor);
                        encoder_->get_joint_candidate(successors_depth_[d][0], d, root, v2_neighbor);
                    }else{
                        encoder_->get_joint_candidate(0, d, history_root, v1_neighbor);
                        encoder_->get_joint_candidate(0, d, root, v2_neighbor);
                    }
                    

                    if(validate_vector_containment_relaxed_full(v2_neighbor, v1_neighbor, history_root, root) == false){
                        contained = false;
                        break;
                    }
                    // uint32_t v1_neighbor_count;
                    // const Vertex* v1_neighbor = data_graph_->getVertexNeighbors(history_root, v1_neighbor_count);
                    // Vertex v1_n;
                    // for(int x=0;x<v1_neighbor_count;++x){
                    //     v1_n = v1_neighbor[x];
                    //     if(data_graph_->getVertexLabel(v1_n)==d && data_graph_->checkEdgeExistence(v1_n, root)== false && root != v1_n){
                    //         contained = false;
                    //         break;
                    //     }
                    // }
                }
                // if(contained==false){
                //     cout<<"failed due to checking"<<depth<<":"<<root<<endl;
                // }
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
                for(auto& p : containment_check){
                    if((explored_depth & ancestors_depth_[p]).any() == false){
                        continue;
                    }
                    // if(p > max_depth){
                    //     continue;
                    // }
                    vector<Vertex>& vec = *(candidates_stack[p].rbegin());
                    for(Vertex v : vec){
                        if(data_graph_->checkEdgeExistence(v, history_root) == true && v != root){
                            if(data_graph_->checkEdgeExistence(v, root) == false){
                                contained = false;
                                break;
                            }
                        }
                    }
                    if(contained == false){
                        break;
                    }
                }
            }
            if(contained == true){
                result_fs = failing_set_recorder_[depth][history_root];
                return true;
            }
        }
    }
    return false;
}
#endif

bool SubtreeReducer::validate_extended_containment(uint32_t depth, vector<vector<vector<Vertex>>>& candidates_stack, Vertex root, FailingSet*& result_fs){
    for(int offset=0;offset<failing_set_recorder_[depth].size();++offset){ // auto offset : conflict_record_offset[depth]
        bool contained = true;
        Vertex history_root = history_root_vertices_[depth][offset];

        for(auto d : unmatched_query_depth_dist_[depth].second){
            if(validate_vector_containment_with_conflict(*(candidates_stack[d].rbegin()), *(candidates_stack[d].rbegin()), history_root, root) == false){
                contained = false;
                break;
            }
        }
        int suc_offset = 0;
        for(auto d : unmatched_query_depth_dist_[depth].first){
            if(validate_vector_containment_with_conflict(history_[depth][offset][suc_offset], *(candidates_stack[d].rbegin()), history_root, root) == false){
                contained = false;
                break;
            }
            suc_offset ++;
        }
        
        FailingSet* fs = failing_set_recorder_[depth][history_root];
        bitset<MAX_QUERY_SIZE> explored_depth = fs->failing_set_ & mask_[depth];

        if(contained == true){
            for(auto d : none_succesors_depth_[depth]){
                if((explored_depth & ancestors_depth_[d]).any() == false){
                    continue;
                }
                for(auto v : *(candidates_stack[d].rbegin())){
                    if(data_graph_->checkEdgeExistence(v, history_root) == true && v != root){
                        if(data_graph_->checkEdgeExistence(v, root) == false){
                            contained = false;
                            break;
                        }
                    }
                }
            }
        }
        
        if(contained == true){
            result_fs = failing_set_recorder_[depth][history_root];
            return true;
        }
    }
    return false;
}

void SubtreeReducer::mark_candidates_to_be_validated(uint32_t depth, Vertex root, FailingSet* fs, vector<vector<vector<Vertex>>>& candidates_stack, bitset<MAX_QUERY_SIZE>& conflict_source, bool full_descendent, uint32_t max_search_depth){
    // is_conflicted[depth].push_back(conflicted);
    history_[depth].push_back({});
    vector<vector<Vertex>>& record = *(history_[depth].rbegin());
    
    for(auto suc_depth : successors_depth_[depth]){
        record.push_back(*(candidates_stack[suc_depth].rbegin()));
    }
    
    history_root_vertices_[depth].push_back(root);

    history_conflict_source_[depth].push_back(conflict_source);
    history_max_search_depth[depth].push_back(max_search_depth);

    if(conflict_source.any() == false){
        empty_record_offset_[depth].push_back(failing_set_recorder_[depth].size());
    }else{
        conflict_record_offset_[depth].push_back(failing_set_recorder_[depth].size());
    }
    failing_set_recorder_[depth].insert({root, fs});
}













bool SubtreeReducer::validate_vector_equivalent(vector<Vertex>& vec1, vector<Vertex>& vec2, Vertex r1, Vertex r2, bool& contain_r2){
    int v1_size = vec1.size();
    int v2_size = vec2.size();
    if(v1_size-v2_size>1 ||v2_size-v1_size>1){
        return false;
    }
    int i=0,j=0;
    for(; i<vec1.size() && j<vec2.size();){
        if(vec2[j]==r2 || vec2[j]==r1){
            ++j;
        }
        if(vec1[i]==r1){
            ++i;
        }else if(vec1[i]==r2){
            contain_r2 = true;
            ++i;
        }
        if(vec1[i]==vec2[j]){
            ++i;++j;
        }else if(vec1[i]>vec2[j]){
            return false;
        }else{
            return false;
        }
    }
    if(j<vec2.size()){
        if(vec2[j]==r1 || vec2[j]==r2){
            return true;
        }
        return false;
    }
    return true;
}



bool SubtreeReducer::validate_vector_containment_with_conflict(vector<Vertex>& vec1, vector<Vertex>& vec2, Vertex r1, Vertex r2){
    // if(vec1.size()<vec2.size()-1){
    //     return false;
    // }
    // int i=0,j=0;
    // bool v1_contained_r2 = false;
    // bool v2_contained_r1 = false;
    // for(; i<vec1.size() && j<vec2.size();){
    //     if(vec2[j]==r2){
    //         ++j;
    //         // prevent overflow
    //         if(j>=vec2.size()){
    //             break;
    //         }
    //     }
    //     if(vec2[j] == r1){
    //         v2_contained_r1 = true;
    //     }
    //     if(vec1[i] == r2){
    //         v1_contained_r2 = true;
    //     }
    //     if(vec1[i]==vec2[j]){
    //         ++i;++j;
    //     }else if(vec1[i]>vec2[j]){
    //         if(vec2[j] == r1){
    //             ++j;
    //             v2_contained_r1 = true;
    //         }else{
    //             return false;
    //         }
    //     }else{
    //         ++i;
    //     }
    // }
    // if(j<vec2.size()){
    //     if(j==vec2.size()-1){
    //         if(vec2[j] == r1){
    //             v2_contained_r1 = true;
    //         }else if(vec2[j] != r2){
    //             return false;
    //         }
    //         if(v2_contained_r1==true && v1_contained_r2==false){
    //             return false;
    //         }
    //         return true;
    //     }
    //     return false;
    // }
    // if(v2_contained_r1==true && v1_contained_r2==false){
    //     return false;
    // }
    // return true;


    if(vec1.size()<vec2.size()-1){
        return false;
    }
    int i=0,j=0;
    bool v1_contained_r2 = false;
    bool v2_contained_r1 = false;
    for(; i<vec1.size() && j<vec2.size();){
        if(vec2[j]==r1){
            ++j;
            v2_contained_r1 = true;
            continue;
        }else if(vec2[j] == r2){
            ++j;
            continue;
        }
        if(vec1[i] == r1){
            ++i;
            continue;
        }else if(vec1[i] == r2){
            ++i;
            v1_contained_r2 = true;
            continue;
        }
        if(i>=vec1.size() || j>=vec2.size()){
            break;
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
    if(j<vec2.size()){
        if(vec2.size() > j){
            return false;
        }else if(vec2.size() == j && (vec2[j] != r1 || vec2[j] != r2)){
            return false;
        }
        if(vec2[j] == r1){
            v2_contained_r1 == true;
        }
    }
    if(v1_contained_r2 == true){
        return true;
    }else{
        if(v2_contained_r1 == true){
            for(i;i<vec1.size();++i){
                if(vec1[i] == r2){
                    return true;
                }
            }
            return false;
        }
        return true;
    }
}

bool SubtreeReducer::validate_vector_containment_relaxed(vector<Vertex>& vec1, vector<Vertex>& vec2, Vertex r1, Vertex r2){
    if(vec1.size()<vec2.size()-1){
        return false;
    }
    int i=0,j=0;
    for(; i<vec1.size() && j<vec2.size();){
        if(vec1[i] == r2){
            ++i;
            if(i>=vec1.size()){
                break;
            }
        }
        if(vec2[j] == r1){
            ++j;
            if(j>=vec2.size()){
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
    
    if(j<vec2.size()){
        if(vec2[j] == r1){
            ++j;
            if(vec2.size() == j){
                return true;
            }
        }
        return false;
    }
    return true;
}

bool SubtreeReducer::validate_vector_containment_relaxed(const Vertex* vec1, uint32_t vec1_size, const Vertex* vec2, uint32_t vec2_size, Vertex r1, Vertex r2){
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

bool SubtreeReducer::validate_vector_containment_relaxed_full(vector<Vertex>& vec1, vector<Vertex>& vec2, Vertex r1, Vertex r2){
    if(vec1.size()<vec2.size()-1){
        return false;
    }
    int i=0,j=0;
    for(; i<vec1.size() && j<vec2.size();){
        if(vec1[i] == r2 || vec1[i] == r1){
            ++i;
            if(i>=vec1.size()){
                break;
            }
        }
        if(vec2[j] == r1 || vec2[j] == r2){
            ++j;
            if(j>=vec2.size()){
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
    
    if(j<vec2.size()){
        if(vec2[j] == r1 || vec2[j] == r2){
            ++j;
            if(vec2.size() == j){
                return true;
            }
        }
        return false;
    }
    return true;
}