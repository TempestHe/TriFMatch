#include "successor_equivalent_cache.h"
#include <set>

SuccessorEquivalentCache::SuccessorEquivalentCache(vector<Vertex>& order, vector<Vertex>& order_index, Graph* query_graph, Graph* data_graph, vector<vector<Vertex>>& successor_neighbors_in_depth, vector<bitset<MAX_QUERY_SIZE>>& ancestors_depth, uint64_t emb_limit){
    order_ = order;
    order_index_ = order_index;
    query_graph_ = query_graph;
    data_graph_ = data_graph;
    successor_neighbors_in_depth_ = successor_neighbors_in_depth;
    ancestors_depth_ = ancestors_depth;
    emb_limit_ = emb_limit;

    query_size_ = query_graph->getVerticesCount();

    // extract the equivalence zone
    equivalent_zones_.clear();
    start_index_ = vector<EquivalenceZone*>(order_.size()+1, NULL);
    end_index_ = vector<EquivalenceZone*>(order_.size()+1, NULL);
    vector<bool> equivalent_bit_map = vector<bool>(order_.size()+1, 0);
    result_count_recorder_ = vector<uint32_t>(order_.size()+1, 0);
    compressed_ = 0;
    for(int depth=(int)(order.size()*max_depth_ratio); depth>1; depth--){
        // cout<<endl<<"///";
        if(equivalent_bit_map[depth] == true){continue;}
        int start_depth = depth;
        EquivalenceZone* zone = new EquivalenceZone;
        zone->equivalent_depth.reset();
        zone->equivalent_depth.set(depth);
        for(int upper_depth = depth-1; upper_depth>0; upper_depth--){
            set<Vertex> end_set(successor_neighbors_in_depth[depth].begin(), successor_neighbors_in_depth[depth].end());
            set<Vertex> start_set(successor_neighbors_in_depth[upper_depth].begin(), successor_neighbors_in_depth[upper_depth].end());
            end_set.erase(upper_depth);
            start_set.erase(depth);
            if(start_set.size() == end_set.size() && query_graph->getVertexLabel(order_[depth]) == query_graph->getVertexLabel(order_[upper_depth])){
                bool equivalent = true;
                for(auto x : end_set){
                    if(start_set.find(x) == start_set.end()){
                        equivalent = false;
                        break;
                    }
                }
                if(equivalent == true){
                    start_depth = upper_depth;
                    equivalent_bit_map[upper_depth] = true;
                    zone->equivalent_depth.set(upper_depth);
                }
            }
        }
        if(start_depth < depth){
            zone->end_depth = depth;
            zone->start_depth = start_depth;
            if(zone->end_depth-zone->start_depth > max_zone_size){
                int max_end_depth = zone->start_depth;
                for(int i=zone->start_depth; i<zone->start_depth+max_zone_size; i++){
                    if(zone->equivalent_depth.test(i) == true){
                        max_end_depth = i;
                    }
                }
                if(max_end_depth != zone->start_depth){
                    zone->end_depth = max_end_depth;
                    start_index_[zone->start_depth] = zone;
                    end_index_[zone->end_depth] = zone;
                    zone->equivalent_depth_map.clear();
                    for(int x=zone->start_depth;x<=zone->end_depth;++x){
                        if(zone->equivalent_depth.test(x) == true){
                            zone->equivalent_depth_map.push_back(x);
                            compressed_++;
                            // cout<<order[x]<<":"<<x<<" ";
                        }
                    }
                    zone->equivalent_relative_offset.clear();
                    int offset = 0;
                    for(int x=zone->start_depth;x<=zone->end_depth;++x){
                        if(zone->equivalent_depth.test(x) == true){
                            zone->equivalent_relative_offset.push_back(offset);
                        }
                        ++offset;
                    }
                    
                    equivalent_zones_.push_back(zone);
                    
                }else{
                    delete zone;
                }
            }else{
                start_index_[start_depth] = zone;
                end_index_[depth] = zone;
                zone->equivalent_depth_map.clear();
                for(int x=zone->start_depth;x<=zone->end_depth;++x){
                    if(zone->equivalent_depth.test(x) == true){
                        zone->equivalent_depth_map.push_back(x);
                        compressed_++;
                        // cout<<order[x]<<":"<<x<<" ";
                    }
                }
                zone->equivalent_relative_offset.clear();
                int offset = 0;
                for(int x=zone->start_depth;x<=zone->end_depth;++x){
                    if(zone->equivalent_depth.test(x) == true){
                        zone->equivalent_relative_offset.push_back(offset);
                    }
                    ++offset;
                }
                
                equivalent_zones_.push_back(zone);
            }
        }else{
            delete zone;
        }
    }
    for(auto zone : equivalent_zones_){
        for(auto d : zone->equivalent_depth_map){
            zone->candidates_to_validate.push_back({});
            uint32_t nei_size;
            const Vertex* nei = query_graph_->getVertexNeighbors(order[d], nei_size);
            for(int x=0;x<nei_size;++x){
                Vertex n = nei[x];
                if(order_index[n] < zone->end_depth){
                    if(order_index[n] < d){
                        zone->candidates_to_validate.rbegin()->push_back(order_index[n]);
                    }else if(zone->equivalent_depth.test(order_index[n]) == false){
                        zone->candidates_to_validate.rbegin()->push_back(order_index[n]);
                    }
                }
            }
        }
    }
    // cout<<endl;
    // cout<<"compressed:"<<compressed_<<endl;
    // neq compression
    // compressed_neq_ = 0;
    // bitset<MAX_QUERY_SIZE> compressed_bitmap_neq;
    // for(int x=1;x<=query_graph_->getVerticesCount();++x){
    //     vector<uint32_t> group;
    //     Vertex ux = order_[x];
    //     group.push_back(ux);
    //     if(compressed_bitmap_neq.test(ux) == true){
    //         continue;
    //     }
    //     compressed_bitmap_neq.set(ux);
    //     for(int y=x+1;y<query_graph_->getVerticesCount();++y){
    //         Vertex uy = order_[y];
    //         if(query_graph_->getVertexLabel(ux) == query_graph_->getVertexLabel(uy) && query_graph_->getVertexDegree(ux) == query_graph_->getVertexDegree(uy)){
    //             bool shared = true;
    //             uint32_t nei_y_size;
    //             const Vertex* nei_uy = query_graph_->getVertexNeighbors(uy, nei_y_size);
    //             for(auto u : group){
    //                 shared = true;
    //                 uint32_t nei_u_size;
    //                 const Vertex* nei_u = query_graph_->getVertexNeighbors(u, nei_u_size);
    //                 if(nei_u_size != nei_y_size){
    //                     shared = false;
    //                     break;
    //                 }
    //                 for(int z=0;z<nei_u_size;++z){
    //                     bool find = false;
    //                     if(nei_u[z] == uy){
    //                         continue;
    //                     }
    //                     for(int t=0;t<nei_y_size;++t){
    //                         if(nei_u[z] == nei_uy[t]){
    //                             find = true;
    //                             break;
    //                         }    
    //                     }
    //                     if(find == false){
    //                         shared = false;
    //                         break;
    //                     }
    //                 }
    //                 if(shared == false){
    //                     break;
    //                 }
    //             }
    //             if(shared == false){
    //                 continue;
    //             }
    //             if(group.size() == 1){
    //                 group.push_back(uy);
    //                 compressed_bitmap_neq.set(uy);
    //             }else{
    //                 if(query_graph->checkEdgeExistence(group[0], group[1])){
    //                     bool connected = true;
    //                     for(auto u : group){
    //                         if(query_graph->checkEdgeExistence(u, uy) == false){
    //                             connected = false;
    //                             break;
    //                         }
    //                     }
    //                     if(connected){
    //                         group.push_back(uy);
    //                         compressed_bitmap_neq.set(uy);
    //                     }
    //                 }else{
    //                     bool connected = false;
    //                     for(auto u : group){
    //                         if(query_graph->checkEdgeExistence(u, uy) == true){
    //                             connected = true;
    //                             break;
    //                         }
    //                     }
    //                     if(connected==false){
    //                         group.push_back(uy);
    //                         compressed_bitmap_neq.set(uy);
    //                     }
    //                 }
    //             }
    //         }
    //     }
    //     if(group.size() > 1){
    //         compressed_neq_ += group.size();
    //         // cout<<"{";
    //         // for(auto v : group){
    //         //     cout<<v<<":"<<order_index[v]<<" ";
    //         // }
    //         // cout<<"}"<<endl;
    //     }
    // }
}

bool SuccessorEquivalentCache::validate_cache_count(uint32_t cur_depth, Vertex* embedding){
    EquivalenceZone* zone = end_index_[cur_depth];
    if(zone == NULL || zone->recorder.empty()){
        return false;
    }
    vector<Vertex> r(embedding+zone->start_depth, embedding+zone->end_depth+1);

    for(int x=0;x<zone->equivalent_relative_offset.size();++x){
        int ex_depth = zone->equivalent_relative_offset[x];
        for(int y=x+1;y<zone->equivalent_relative_offset.size();++y){
            int ey_depth = zone->equivalent_relative_offset[y];
            if(r[ex_depth] > r[ey_depth]){
                swap(r[ex_depth], r[ey_depth]);
            }
        }
    }
    // unordered_map<vector<Vertex>, MatchInfo, HashFunc_Vec>& recorder = zone->recorder_emb_count;
    auto itf = zone->recorder.find(r);
    if(itf != zone->recorder.end()){
        // result_count = itf->second.emb_count;


        // failing_set_depth = itf->second.empty_set_failure;
        // empty_failure_depth = itf->second.empty_set_failure;
        // for(auto v : itf->second.conflict_vertex){
        //     int d = visited_query_vertices[v];
        //     conflict_source_depth.set(d);
        //     failing_set_depth |= ancestors_depth_[d];
        // }
        // for(auto d : itf->second.conflict_depth_outer){
        //     conflict_source_depth.set(d);
        //     failing_set_depth |= ancestors_depth_[d];
        // }

        // // failing_set_depth = itf->second.empty_set_failure;
        // // for(auto v : itf->second.conflict_vertex){
        // //     int d = visited_query_vertices[v];
        // //     conflict_source_depth.set(d);
        // //     failing_set_depth |= ancestors_depth_[d];
        // // }
        return true;
    }
    return false;
}


// bool SuccessorEquivalentCache::validate_cache_result(uint32_t cur_depth, Vertex* embedding, bitset<MAX_QUERY_SIZE>& conflict_source_depth, bitset<MAX_QUERY_SIZE>& empty_failure_depth, bitset<MAX_QUERY_SIZE>& failing_set_depth, uint64_t& result_count, Vertex* visited_query_vertices, vector<vector<Vertex>>& results){
//     EquivalenceZone* zone = end_index_[cur_depth];
//     if(zone == NULL){
//         return false;
//     }
//     vector<Vertex> r(embedding+zone->start_depth, embedding+zone->end_depth+1);

//     for(int x=0;x<zone->equivalent_relative_offset.size();++x){
//         int ex_depth = zone->equivalent_relative_offset[x];
//         for(int y=x+1;y<zone->equivalent_relative_offset.size();++y){
//             int ey_depth = zone->equivalent_relative_offset[y];
//             if(r[ex_depth] > r[ey_depth]){
//                 swap(r[ex_depth], r[ey_depth]);
//             }
//         }
//     }
//     unordered_map<vector<Vertex>, vector<vector<Vertex>>, HashFunc_Vec> recorder_emb_result = zone->recorder_emb_result;
//     auto itf = recorder_emb_result.find(r);
//     if(itf != recorder_emb_result.end()){
//         for(auto res : itf->second){
//             memcpy(&(res[0]), embedding, sizeof(Vertex)*(zone->end_depth+1));
//             for(int i=1;i<res.size();++i){
//                 cout<<res[i]<<" ";
//             }
//             cout<<endl;
//             results.push_back(res);
//         }
//         return true;
//     }
//     return false;
// }

void SuccessorEquivalentCache::start_insert_cache(uint32_t cur_depth, uint64_t& result_count){
    EquivalenceZone* zone = end_index_[cur_depth];
    if(zone == NULL){
        return;
    }
    result_count_recorder_[cur_depth] = result_count;
}

bool SuccessorEquivalentCache::finish_insert_cache_count(uint32_t cur_depth, Vertex* embedding, uint64_t& result_count, Vertex* visited_query_vertices
#if PRINT_RESULT==1
    , vector<vector<Vertex>>& matches
#endif
    ){
    EquivalenceZone* zone = end_index_[cur_depth];
    if(zone == NULL){
        return false;
    }
    int equivalent_size = zone->equivalent_relative_offset.size();
    vector<Vertex> r(embedding+zone->start_depth, embedding+zone->end_depth+1);
    vector<Vertex> swappable_candidates;
    for(int x=0;x<equivalent_size;++x){
        int ex_depth = zone->equivalent_relative_offset[x];
        for(int y=x+1;y<equivalent_size;++y){
            int ey_depth = zone->equivalent_relative_offset[y];
            if(r[ex_depth] > r[ey_depth]){
                swap(r[ex_depth], r[ey_depth]);
            }
        }
    }
    // MatchInfo info;
    // info.emb_count = result_count-result_count_recorder_[cur_depth];
    // // info.empty_set_failure = empty_failure_depth;
    // // info.empty_set_failure = failing_set_depth;
    // // for(int x=1;x<=query_size_;++x){
    // //     if(conflict_source_depth.test(x) == true){
    // //         if(zone->equivalent_depth.test(x) == true){
    // //             info.conflict_vertex.push_back(embedding[x]);
    // //         }else{
    // //             info.conflict_depth_outer.push_back(x);
    // //         }
    // //     }
    // // }
    
    // unordered_map<vector<Vertex>, MatchInfo, HashFunc_Vec>& recorder = zone->recorder_emb_count;
    // if(recorder.size() > size_limit_){
    //     recorder.erase(recorder.begin());
    // }
    // recorder.insert({r, info});

    // start amplifying
    for(auto d : zone->equivalent_depth_map){
        swappable_candidates.push_back(embedding[d]);
    }
    uint32_t count = result_count-result_count_recorder_[cur_depth];
    sort(swappable_candidates.begin(), swappable_candidates.end());
    vector<vector<Vertex>> candidates;
    unordered_set<Vertex> matched;
    vector<Vertex> local_emb(embedding, embedding+query_size_+1);
    vector<Vertex> candidate_offset = vector<Vertex>(equivalent_size+1, 0);
    candidates.resize(equivalent_size+1);
    // initialize the first level candidates
    for(auto c : swappable_candidates){
        bool valid = true;
        for(auto d : zone->candidates_to_validate[0]){
            if(data_graph_->checkEdgeExistence(c, embedding[d]) == false){
                valid = false;
                break;
            }
        }
        if(valid == true){
            candidates[0].push_back(c);
        }
    }
    int local_depth = 0;
    uint32_t hit_count = 0;
    uint64_t extra_count = result_count;
    while(true){
        while(candidate_offset[local_depth] < candidates[local_depth].size()){
            if(local_depth == equivalent_size-1){
                // print results
#if PRINT_RESULT==1
                for(auto c : candidates[local_depth]){
                    local_emb[zone->equivalent_depth_map[local_depth]] = c;
                    bool same = true;
                    for(auto d : zone->equivalent_depth_map){
                        if(embedding[d] != local_emb[d]){
                            same = false;
                            break;
                        }
                    }
                    if(same == true){
                        continue;
                    }
                    for(int x=result_count_recorder_[cur_depth]; x<result_count_recorder_[cur_depth]+count; ++x){
                        memcpy(&(local_emb[zone->end_depth+1]), &(matches[x][zone->end_depth+1]), sizeof(Vertex)*(query_size_-zone->end_depth));
                        
                        for(int i=1;i<=query_size_;++i){
                            cout<<local_emb[i]<<" ";
                        }
                        cout<<endl;
                        matches.push_back(local_emb);
                    }
                    result_count += count;
                    if(result_count >= emb_limit_){
                        return true;
                    }
                }
#else
                // for(auto c : candidates[local_depth]){
                //     // local_emb[zone->equivalent_depth_map[local_depth]] = c;
                //     // bool same = true;
                //     // for(auto d : zone->equivalent_depth_map){
                //     //     if(embedding[d] != local_emb[d]){
                //     //         same = false;
                //     //         break;
                //     //     }
                //     // }
                //     // if(same == true){
                //     //     continue;
                //     // }else{
                //     //     result_count += count;
                //     //     cout<<"augmented:"<<count<<endl;
                //     // }
                //     result_count += count;
                // }
                hit_count += candidates[local_depth].size();
                extra_count += candidates[local_depth].size()*count;
                if(extra_count-count >= emb_limit_){
                    result_count = extra_count-count;
                    return true;
                }
                // cout<<"count:"<<result_count+hit_count*count<<endl;
                // result_count += count*candidates[local_depth].size();
#endif
                // result_count += count*(candidates[local_depth].size()-1); // -1 prevent counting redundantly
                // if(candidates[local_depth].size()-1 > 0){
                //     cout<<"augmented"<<endl;
                // }
                // if(result_count >= emb_limit_){
                //     return true;
                // }
                candidate_offset[local_depth] = candidates[local_depth].size();
                continue;
            }
            Vertex v = candidates[local_depth][candidate_offset[local_depth]];
            matched.insert(v);
            local_emb[zone->equivalent_depth_map[local_depth]] = v;
            candidate_offset[local_depth]++;
            local_depth++;
            candidate_offset[local_depth] = 0;
            candidates[local_depth].clear();
            for(auto c:swappable_candidates){
                if(matched.find(c) == matched.end()){
                    bool valid = true;
                    for(auto d : zone->candidates_to_validate[local_depth]){
                        if(data_graph_->checkEdgeExistence(c, local_emb[d]) == false){
                            valid = false;
                            break;
                        }
                    }
                    if(valid == true){
                        candidates[local_depth].push_back(c);
                    }
                }
            }
        }
        local_depth--;
        matched.erase(local_emb[zone->equivalent_depth_map[local_depth]]);
        if(local_depth==-1){
            break;
        }
    }
#if PRINT_RESULT == 0
    //remove redundancy
    if(hit_count > 1){
        result_count += (hit_count-1)*count;
        zone->recorder.insert(r);
    }
    // result_count -= count;
#else
    zone->recorder.insert(r);
#endif
    // MatchInfo info;
    // unordered_map<vector<Vertex>, MatchInfo, HashFunc_Vec>& recorder = zone->recorder_emb_count;
    // recorder.insert({r, info});
    return false;
}

// void SuccessorEquivalentCache::finish_insert_cache_result(uint32_t cur_depth, Vertex* embedding, bitset<MAX_QUERY_SIZE>& conflict_source_depth, bitset<MAX_QUERY_SIZE>& empty_failure_depth, bitset<MAX_QUERY_SIZE>& failing_set_depth, uint64_t& result_count, Vertex* visited_query_vertices, vector<vector<Vertex>>& results){
//     EquivalenceZone* zone = end_index_[cur_depth];
//     if(zone == NULL){
//         return;
//     }
//     vector<Vertex> r(embedding+zone->start_depth, embedding+zone->end_depth+1);

//     for(int x=0;x<zone->equivalent_relative_offset.size();++x){
//         int ex_depth = zone->equivalent_relative_offset[x];
//         for(int y=x+1;y<zone->equivalent_relative_offset.size();++y){
//             int ey_depth = zone->equivalent_relative_offset[y];
//             if(r[ex_depth] > r[ey_depth]){
//                 swap(r[ex_depth], r[ey_depth]);
//             }
//         }
//     }
//     unordered_map<vector<Vertex>, vector<vector<Vertex>>, HashFunc_Vec>& recorder_emb_result = zone->recorder_emb_result;
//     if(recorder_emb_result.size() > size_limit_){
//         recorder_emb_result.erase(recorder_emb_result.begin());
//     }
//     recorder_emb_result.insert({r, {}});
//     auto itf = recorder_emb_result.find(r);
//     auto itf_de = zone->recorder_emb_count.find(r);
//     for(int i=result_count_recorder_[cur_depth]; i<result_count; ++i){
//         itf->second.push_back(results[i]);
//     }
// }

void SuccessorEquivalentCache::clear(uint32_t cur_depth){
    EquivalenceZone* zone = start_index_[cur_depth];
    if(start_index_[cur_depth] == NULL){
        return;
    }
    zone->recorder.clear();
//     zone->recorder_emb_count.clear();
// #if PRINT_RESULT == 1
//     zone->recorder_emb_result.clear();
// #endif
}

SuccessorEquivalentCache::~SuccessorEquivalentCache(){
    for(auto zone : equivalent_zones_){
        delete zone;
    }
}

