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

void CandidatesHistoryStack::rebuild(Graph* query_graph, Graph* data_graph, vector<Vertex>& successors, uint32_t cur_depth){
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

void CandidatesHistoryStack::append_history(Vertex v, TrieEncoder* encoder, vector<vector<pair<Vertex*, uint32_t>>>& candidates_stack){
    history_root_vertices_.push_back(v);
    history_failing_set_.push_back(bitset<MAX_QUERY_SIZE>());
    // history_conflict_sources_.push_back(bitset<MAX_QUERY_SIZE>());
    history_candidates_.push_back({});
    // validation_marker_.push_back(0);
    vector<pair<Vertex*, uint32_t>>& scope = *(history_candidates_.rbegin());
    for(auto suc_depth : successors_){
        uint32_t count;
        Vertex* cans = encoder->get_edge_candidate(cur_depth_, suc_depth, v, count);
        if(candidates_stack[suc_depth].empty()){
            candidates_stack[suc_depth].push_back({cans, count});
            scope.push_back({cans, count});
        }else{
            auto& last_ele = *(candidates_stack[suc_depth].rbegin());
            if(last_ele.second >= last_block_remaining_size_){
                current_block_ptr = new Vertex [block_size_];
                total_block_size_ += block_size_;
                last_block_remaining_size_ = block_size_;
                // current_block_ptr = new Vertex [last_ele.second];
                // total_block_size_ += last_ele.second;
                // last_block_remaining_size_ = last_ele.second;
                block_content_.push_back(current_block_ptr);
            }
            uint32_t result_count;
            ComputeSetIntersection::ComputeCandidates(cans, count, last_ele.first, last_ele.second, current_block_ptr, result_count);
            candidates_stack[suc_depth].push_back({current_block_ptr, result_count});
            scope.push_back({current_block_ptr, result_count});
            current_block_ptr += result_count;
            last_block_remaining_size_ -= result_count;
        }
    }
}

void CandidatesHistoryStack::clear(){
    history_root_vertices_.clear();
    history_failing_set_.clear();
    // history_conflict_sources_.clear();
    history_candidates_.clear();
    // validation_marker_.clear();
    validate_exclusion_offsets_.clear();
    validate_completion_offset_.clear();
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
}