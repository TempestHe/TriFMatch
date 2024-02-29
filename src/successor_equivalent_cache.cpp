// #include "successor_equivalent_cache.h"

// SuccessorEquivalentCache::SuccessorEquivalentCache(vector<uint32_t>& equivalent_set_depth){
//     equivalent_set_depth_ = equivalent_set_depth;
//     start_depth_ = *(equivalent_set_depth.rbegin());
//     end_depth_ = equivalent_set_depth[0];
// }

// void SuccessorEquivalentCache::insert_emb_count(Vertex* emb, FailingSet* fs, uint32_t emb_count){
//     vector<Vertex> r(emb+start_depth_, emb+end_depth_+1);
//     sort(r.begin(), r.end());
//         MatchInfo info;
//         info.emb_count = emb_count;
//         info.fs = *fs;
// #if SUBTREE_REDUCTION_FILTERING==1
//         for(int d=start_depth_;d<=end_depth_;++d){
//             if(fs->conflict_set_depth_.test(d) == true){
//                 info.conflicted_vertices.push_back(emb[d]);
//             }
//         }
// #endif
//         recoder_emb_count_.insert({r, info});
// }

// bool SuccessorEquivalentCache::validate_emb_count(Vertex* embedding, FailingSet*& fs, uint32_t& emb_count, Vertex* visited_query_vertices, vector<Vertex>& order_index, vector<bitset<MAX_QUERY_SIZE>>& ancestors_depth){
//     vector<Vertex> r(embedding+start_depth_, embedding+end_depth_+1);
//     sort(r.begin(), r.end());
//     auto itf = recoder_emb_count_.find(r);
//     if(itf != recoder_emb_count_.end()){
//         fs = &(itf->second.fs);
// #if SUBTREE_REDUCTION_FILTERING==1
//         for(int d=start_depth_;d<=end_depth_;++d){
//             fs->conflict_set_depth_.reset(d);
//         }
//         for(auto v : itf->second.conflicted_vertices){
//             fs->conflict_set_depth_.set(visited_query_vertices[v]);
//             fs->failing_set_ |= ancestors_depth[visited_query_vertices[v]];
//         }
// #endif
//         emb_count = itf->second.emb_count;
//         return true;
//     }
//     return false;
// }

// #if PRINT_RESULT==1
// void SuccessorEquivalentCache::insert_emb_result(Vertex* emb, vector<vector<Vertex>>& matches, uint32_t starting_emb_count, uint32_t emb_count){
//     vector<Vertex> r(emb+start_depth_, emb+end_depth_+1);
//     sort(r.begin(), r.end());
//     vector<vector<Vertex>>& result = recoder_emb_result_[r];
//     for(int i=starting_emb_count; i<emb_count; ++i){
//         result.push_back(matches[i]);
//     }
// }

// bool SuccessorEquivalentCache::validate_emb_result(Vertex* embedding, vector<vector<Vertex>>& matches){
//     vector<Vertex> r(embedding+start_depth_, embedding+end_depth_+1);
//     sort(r.begin(), r.end());
//     auto itf = recoder_emb_result_.find(r);

//     if(itf != recoder_emb_result_.end()){
//         for(auto res : itf->second){
//             memcpy(&(res[0]), embedding, sizeof(Vertex)*(end_depth_+1));

//             for(int i=1;i<res.size();++i){
//                 cout<<res[i]<<" ";
//             }
//             cout<<endl;
//             matches.push_back(res);
//         }
//         return true;
//     }
//     return false;
// }
// #endif

// void SuccessorEquivalentCache::clear(){
//     recoder_emb_count_.clear();
// #if PRINT_RESULT==1
//     recoder_emb_result_.clear();
// #endif
// }


// // =============================================================

// SuccessorEquivalentSet::SuccessorEquivalentSet(Graph* query_graph, vector<Vertex>& order, vector<Vertex>& order_index, vector<bitset<MAX_QUERY_SIZE>>& ancestors_depth){
//     order_ = order;
//     order_index_ = order_index;
//     ancestors_depth_ = ancestors_depth;
//     emb_count_recorder_ = vector<uint32_t>(order.size(), 0);

//     query_vertex_count_ = query_graph->getVerticesCount();
//     end_index_ = vector<SuccessorEquivalentCache*>(order.size(), NULL);
//     start_index_.resize(order.size());
    
//     vector<uint32_t> equivalent_set_depth;
//     unordered_set<uint32_t> explored;
//     for(int x=order.size()-2;x>=1;--x){
//         equivalent_set_depth.clear();
//         equivalent_set_depth.push_back(x);
        
//         for(int c = x-1;c>=1;--c){
//             equivalent_set_depth.push_back(c);
//             // validate whether their successors are equal
//             bool valid = true;
//             if(query_graph->getVertexLabel(order[c]) == query_graph->getVertexLabel(order[x])){
//                 // int frontier_depth = (equivalent_set_depth[0]==order.size()-1) ? equivalent_set_depth[0]-1 : equivalent_set_depth[0];
//                 for(int y=order.size()-1;y>x;--y){
//                     int connected = 0;
//                     for(auto m : equivalent_set_depth){
//                         if(connected == 0){
//                             if(query_graph->checkEdgeExistence(order[y], order[m]) == true){
//                                 connected = 1;
//                             }else{
//                                 connected = -1;
//                             }
//                         }else{
//                             if(query_graph->checkEdgeExistence(order[y], order[m]) == true && connected==-1){
//                                 valid = false;
//                                 break;
//                             }else if(query_graph->checkEdgeExistence(order[y], order[m]) == false && connected==1){
//                                 valid = false;
//                                 break;
//                             }
//                         }
//                     }
//                     if(valid == false){
//                         break;
//                     }
//                 }
//             }else{
//                 valid = false;
//             }
//             if(valid == false){
//                 equivalent_set_depth.pop_back();
//                 if(equivalent_set_depth.size() > 1){
//                     for(auto d : equivalent_set_depth){
//                         explored.insert(d);
//                     }
//                     if(equivalent_set_depth.size() < order_.size()-equivalent_set_depth[0]){
//                         // for(auto d : equivalent_set_depth){
//                         //     cout<<d<<" ";
//                         // }
//                         // cout<<endl;
//                         SuccessorEquivalentCache* cache = new SuccessorEquivalentCache(equivalent_set_depth);
//                         end_index_[equivalent_set_depth[0]] = cache;
//                         start_index_[*(equivalent_set_depth.rbegin())].push_back(cache);
//                     }
                    
//                 }
//                 break;
//             }
//         }
//         explored.insert(x);
//     }
// }

// void SuccessorEquivalentSet::clear_recording(int depth){
//     for(auto cache : start_index_[depth+1]){
//         cache->clear();
//     }
// }

// bool SuccessorEquivalentSet::validate_recording(Vertex* embedding, FailingSet*& fs, int depth, Vertex u, uint32_t& emb_count, Vertex* visited_query_vertices){
//     SuccessorEquivalentCache* cache = end_index_[depth];
//     if(cache != NULL){
//         return cache->validate_emb_count(embedding, fs, emb_count, visited_query_vertices, order_index_, ancestors_depth_);
//     }
//     return false;
// }

// void SuccessorEquivalentSet::start_recording(uint32_t emb_count, int depth, Vertex u){
//     emb_count_recorder_[depth] = emb_count;
// }

// void SuccessorEquivalentSet::stop_recording(Vertex* embedding, FailingSet* fs, uint32_t emb_count, int depth){
//     SuccessorEquivalentCache* cache = end_index_[depth];
//     if(cache != NULL){
//         cache->insert_emb_count(embedding, fs, emb_count-emb_count_recorder_[depth]);
//     }
// }

// #if PRINT_RESULT == 1
// void SuccessorEquivalentSet::recording_matches(Vertex* embedding, uint32_t emb_count, int depth, vector<vector<Vertex>>& matches){
//     SuccessorEquivalentCache* cache = end_index_[depth];
//     if(cache != NULL){
//         cache->insert_emb_result(embedding, matches, emb_count_recorder_[depth], emb_count);
//     }
// }

// void SuccessorEquivalentSet::validate_recording_matches(Vertex* embedding, int depth, Vertex u, vector<vector<Vertex>>& matches){
//     SuccessorEquivalentCache* cache = end_index_[depth];
//     if(cache != NULL){
//         cache->validate_emb_result(embedding, matches);
//     }
// }
// #endif

// SuccessorEquivalentSet::~SuccessorEquivalentSet(){
//     for(auto cache : end_index_){
//         if(cache != NULL){
//             delete cache;
//         }
//     }
// }
