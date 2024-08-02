#include "candidates_recorder.h"

CandidatesHistoryStack::CandidatesHistoryStack(Graph* query_graph, Graph* data_graph, vector<Vertex>& successors, uint32_t cur_depth){
    successors_ = successors;
    cur_depth_ = cur_depth;
    block_size_ = data_graph->getGraphMaxDegree();
    last_block_remaining_size_ = block_size_;
    current_block_ptr = new Vertex [block_size_];
    total_block_size_ = block_size_;

    if(successors.size() > 0){
        int max_successor_depth = *(max_element(successors.begin(), successors.end()));
    }
    block_content_.push_back(current_block_ptr);

}

void CandidatesHistoryStack::rebuild(Graph* query_graph, Graph* data_graph, vector<Vertex>& successors, vector<vector<Vertex>>& predecessors, uint32_t cur_depth){
    successors_ = successors;
    cur_depth_ = cur_depth;
    block_size_ = data_graph->getGraphMaxDegree();
    last_block_remaining_size_ = block_size_;
    current_block_ptr = new Vertex [block_size_];
    total_block_size_ = block_size_;

#if LOOKAHEAD_OVERLAP == 0
    local_predecessors_.clear();
    intersection_buffer_ = new Vertex [block_size_];
    intersection_buffer_tmp_ = new Vertex [block_size_];
#endif

    if(successors.size() > 0){
        int max_successor_depth = *(max_element(successors.begin(), successors.end()));
#if LOOKAHEAD_OVERLAP == 0
        local_predecessors_.resize(max_successor_depth+1);
        for(auto suc_depth : successors){
            for(auto pred_depth : predecessors[suc_depth]){
                if(pred_depth <= cur_depth){
                    local_predecessors_[suc_depth].push_back(pred_depth);
                }
            }
        }
#endif
    }
    block_content_.push_back(current_block_ptr);
}

void CandidatesHistoryStack::append_history(Vertex v, TrieEncoder* encoder, vector<vector<pair<Vertex*, uint32_t>>>& candidates_stack, vector<vector<bitset<MAX_QUERY_SIZE>>>& failing_set_stack
#if LOOKAHEAD_OVERLAP == 0
    , vector<Vertex>& embedding_depth
#endif
){
#if CD_FILTERING == 1
    history_root_vertices_.push_back(v);
    history_failing_set_.push_back(bitset<MAX_QUERY_SIZE>());
    history_candidates_.push_back({});
    vector<pair<Vertex*, uint32_t>>& scope = *(history_candidates_.rbegin());
#endif
    for(auto suc_depth : successors_){
        uint32_t count;
        Vertex* cans = encoder->get_edge_candidate(cur_depth_, suc_depth, v, count);
        if(candidates_stack[suc_depth].empty()){
            candidates_stack[suc_depth].push_back({cans, count});
#if CD_FILTERING == 1
            scope.push_back({cans, count});
#endif
#if FAILING_SET==1
            failing_set_stack[suc_depth].push_back(*(failing_set_stack[suc_depth].rbegin()) | *(failing_set_stack[cur_depth_].rbegin()));
#endif
        }else{
#if LOOKAHEAD_OVERLAP == 0
            if(count>= last_block_remaining_size_){
                current_block_ptr = new Vertex [block_size_];
                total_block_size_ += block_size_;
                last_block_remaining_size_ = block_size_;
                block_content_.push_back(current_block_ptr);
            }
            assert(local_predecessors_[suc_depth].size() > 0);
            if(local_predecessors_[suc_depth].size() == 1){
                uint32_t pred_depth = local_predecessors_[suc_depth][0];
                uint32_t count;
                Vertex* cans = encoder->get_edge_candidate(pred_depth, suc_depth, embedding_depth[pred_depth], count);
                memcpy(current_block_ptr, cans, sizeof(Vertex)*count);
                candidates_stack[suc_depth].push_back({current_block_ptr, count});
#if CD_FILTERING == 1
                scope.push_back({current_block_ptr, count});
#endif
                current_block_ptr += count;
                last_block_remaining_size_ -= count;
            }else{
                for(uint32_t x=1;x<local_predecessors_[suc_depth].size();++x){
                    uint32_t pred_depth = local_predecessors_[suc_depth][x];
                    uint32_t count;
                    Vertex* cans = encoder->get_edge_candidate(pred_depth, suc_depth, embedding_depth[pred_depth], count);
                    if(x == 1){
                        uint32_t prepre_depth = local_predecessors_[suc_depth][0];
                        uint32_t pre_count;
                        Vertex* cans_pre = encoder->get_edge_candidate(prepre_depth, suc_depth, embedding_depth[prepre_depth], pre_count);
                        ComputeSetIntersection::ComputeCandidates(cans, count, cans_pre, pre_count, intersection_buffer_, intersection_buffer_size_);
                    }else{
                        ComputeSetIntersection::ComputeCandidates(cans, count, intersection_buffer_, intersection_buffer_size_, intersection_buffer_tmp_, intersection_buffer_size_tmp_);
                        swap(intersection_buffer_, intersection_buffer_tmp_);
                        swap(intersection_buffer_size_, intersection_buffer_size_tmp_);
                    }
                }
                memcpy(current_block_ptr, intersection_buffer_, sizeof(Vertex)*intersection_buffer_size_);
                candidates_stack[suc_depth].push_back({current_block_ptr, intersection_buffer_size_});
#if CD_FILTERING == 1
                scope.push_back({current_block_ptr, intersection_buffer_size_});
#endif
                current_block_ptr += intersection_buffer_size_;
                last_block_remaining_size_ -= intersection_buffer_size_;
            }
#else
            auto& last_ele = *(candidates_stack[suc_depth].rbegin());
            if(last_ele.second >= last_block_remaining_size_){
                current_block_ptr = new Vertex [block_size_];
                total_block_size_ += block_size_;
                last_block_remaining_size_ = block_size_;
                block_content_.push_back(current_block_ptr);
            }
            uint32_t result_count;
            ComputeSetIntersection::ComputeCandidates(cans, count, last_ele.first, last_ele.second, current_block_ptr, result_count);
            candidates_stack[suc_depth].push_back({current_block_ptr, result_count});
#if CD_FILTERING == 1
            scope.push_back({current_block_ptr, result_count});
#endif
            current_block_ptr += result_count;
            last_block_remaining_size_ -= result_count;
#endif

#if FAILING_SET==1
            uint32_t can_s = candidates_stack[suc_depth].size();
            if(can_s > 1 && candidates_stack[suc_depth][can_s-1].second==candidates_stack[suc_depth][can_s-2].second){ // 
                failing_set_stack[suc_depth].push_back(*(failing_set_stack[suc_depth].rbegin()));
            }else{
                failing_set_stack[suc_depth].push_back(*(failing_set_stack[suc_depth].rbegin()) | *(failing_set_stack[cur_depth_].rbegin()));
            }
            assert(failing_set_stack[suc_depth].size()-1 == candidates_stack[suc_depth].size());
#endif

        }
    }
}




void CandidatesHistoryStack::reorder_searching_candidates(vector<Vertex>& ordered_cans, Vertex* cans, uint32_t cans_count, TrieEncoder* encoder, vector<vector<pair<Vertex*, uint32_t>>>& candidates_stack){
    vector<uint32_t> card_vec;
    for(int x=0;x<cans_count;++x){
        Vertex v = cans[x];
        ordered_cans.push_back(v);
        history_candidates_.push_back({});
        vector<pair<Vertex*, uint32_t>>& scope = *(history_candidates_.rbegin());
        uint32_t card = 0;
        for(auto suc_depth : successors_){
            uint32_t count;
            Vertex* cans = encoder->get_edge_candidate(cur_depth_, suc_depth, v, count);
            if(candidates_stack[suc_depth].empty()){
                scope.push_back({cans, count});
                card += count;
            }else{
                auto& last_ele = *(candidates_stack[suc_depth].rbegin());
                if(last_ele.second >= last_block_remaining_size_){
                    current_block_ptr = new Vertex [block_size_];
                    total_block_size_ += block_size_;
                    last_block_remaining_size_ = block_size_;
                    block_content_.push_back(current_block_ptr);
                }
                uint32_t result_count;
                ComputeSetIntersection::ComputeCandidates(cans, count, last_ele.first, last_ele.second, current_block_ptr, result_count);
                scope.push_back({current_block_ptr, result_count});
                current_block_ptr += result_count;
                last_block_remaining_size_ -= result_count;
                card += result_count;
            }
        }
        card_vec.push_back(card);
    }
    // reorder
    for(int x=0;x<cans_count;++x){
        for(int y=x;y<cans_count;++y){
            if(card_vec[x]<card_vec[y]){
                swap(card_vec[x], card_vec[y]);
                swap(ordered_cans[x], ordered_cans[y]);
                swap(history_candidates_[x], history_candidates_[y]);
            }
        }
    }
}


void CandidatesHistoryStack::reorder_searching_init(Vertex* cans, uint32_t cans_count, TrieEncoder* encoder){
    for(int x=0;x<cans_count;++x){
        Vertex v = cans[x];
        history_candidates_.push_back({});
        vector<pair<Vertex*, uint32_t>>& scope = *(history_candidates_.rbegin());
        for(auto suc_depth : successors_){
            uint32_t count;
            Vertex* cans = encoder->get_edge_candidate(cur_depth_, suc_depth, v, count);
            scope.push_back({cans, count});
        }
    }
}


void CandidatesHistoryStack::append_history_reorded_strict(Vertex v, uint32_t offset, vector<vector<pair<Vertex*, uint32_t>>>& candidates_stack, vector<vector<bitset<MAX_QUERY_SIZE>>>& failing_set_stack){
    history_root_vertices_.push_back(v);
    history_failing_set_.push_back(bitset<MAX_QUERY_SIZE>());
    int suc_offset = 0;
    for(auto suc_depth : successors_){
        candidates_stack[suc_depth].push_back(history_candidates_[offset][suc_offset]);

#if FAILING_SET==1
        uint32_t can_s = candidates_stack[suc_depth].size();
        if(can_s > 1 && candidates_stack[suc_depth][can_s-1].second==candidates_stack[suc_depth][can_s-2].second){ // 
            failing_set_stack[suc_depth].push_back(*(failing_set_stack[suc_depth].rbegin()));
        }else{
            failing_set_stack[suc_depth].push_back(*(failing_set_stack[suc_depth].rbegin()) | *(failing_set_stack[cur_depth_].rbegin()));
        }
        assert(failing_set_stack[suc_depth].size()-1 == candidates_stack[suc_depth].size());
#endif
        suc_offset ++;
    }

}

void CandidatesHistoryStack::clear(){
#if CD_FILTERING == 1
    history_root_vertices_.clear();
    history_failing_set_.clear();
    history_candidates_.clear();
    validate_exclusion_offsets_.clear();
    validate_completion_offset_.clear();
#endif
    if(block_content_.size() > 1){
        for(auto blk : block_content_){
            delete blk;
        }
        total_block_size_ = total_block_size_*2;
        last_block_remaining_size_ = total_block_size_;
        current_block_ptr = new Vertex [last_block_remaining_size_];
        block_content_.clear();
        block_content_.push_back(current_block_ptr);
    }else{
        current_block_ptr = block_content_[0];
        last_block_remaining_size_ = total_block_size_;
    }
}

CandidatesHistoryStack::~CandidatesHistoryStack(){
    for(auto block : block_content_){
        delete [] block;
    }
#if LOOKAHEAD_OVERLAP==0
    if(intersection_buffer_ != NULL)
        delete [] intersection_buffer_;
    if(intersection_buffer_tmp_ != NULL)
        delete [] intersection_buffer_tmp_;
#endif
}