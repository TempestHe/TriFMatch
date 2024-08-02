#include "subgraph_enumeration.h"
#include "../utility/computesetintersection.h"
#include <chrono>
#include "signal.h"
#include "time.h"

timer_t id;
void set_stop_flag_enum(union sigval val){
    *((bool*)val.sival_ptr) = true;
    timer_delete(id);
}

void register_timer_enum(bool* stop_flag){
    struct timespec spec;
    struct sigevent ent;
    struct itimerspec value;
    struct itimerspec get_val;

    /* Init */
    memset(&ent, 0x00, sizeof(struct sigevent));
    memset(&get_val, 0x00, sizeof(struct itimerspec));

    int test_val = 0;
    /* create a timer */
    ent.sigev_notify = SIGEV_THREAD;
    ent.sigev_notify_function = set_stop_flag_enum;
    ent.sigev_value.sival_ptr = stop_flag;
    timer_create(CLOCK_MONOTONIC, &ent, &id);

    /* start a timer */
    value.it_value.tv_sec = TIME_LIMIT;
    value.it_value.tv_nsec = 0;
    value.it_interval.tv_sec = 0;
    value.it_interval.tv_nsec = 0;
    timer_settime(id, 0, &value, NULL);
}

#if PRINT_MEM_INFO == 1
inline int GetCurrentPid(){
    return getpid();
}

void thread_get_mem_info(float& peak_memory, bool& stop){
    peak_memory = 0.0; // MB
    float read_memory;
    int current_pid = GetCurrentPid();
    while (stop==false)
    {
        read_memory = GetMemoryUsage(current_pid);
        peak_memory = (peak_memory>read_memory) ? peak_memory : read_memory;
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
}

inline float GetMemoryUsage(int pid){
    char file_name[64] = { 0 };
    FILE* fd;
    char line_buff[512] = { 0 };
    sprintf(file_name, "/proc/%d/status", pid);

    fd = fopen(file_name, "r");
    if (nullptr == fd)
        return 0;

    char name[64];
    int vmrss = 0;
    for (int i = 0; i < VMRSS_LINE - 1; i++)
        fgets(line_buff, sizeof(line_buff), fd);

    fgets(line_buff, sizeof(line_buff), fd);
    sscanf(line_buff, "%s %d", name, &vmrss);
    fclose(fd);

    // cnvert VmRSS from KB to MB
    return vmrss / 1024.0;
}
#endif



SubgraphEnum::SubgraphEnum(Graph* data_graph){
    data_graph_ = data_graph;

    storage_ = NULL;
    pp_ = NULL;

}

bool debug(int depth, Vertex* emb, vector<Vertex>& order){
    //560,508855,97672,97662,97668,97652,508904,511522,97650,97648,190091,97655,510309,190092,515644,508857,80602,97664,511287,512715,97649,97653,522241,514799,16479,538317,80600,515790,522225,538839,80597,539943
    // vector<Vertex> content = {0, 13985,13928,10141,3325,13940,4298,13963,14001,13970,13931,10151};
    vector<Vertex> content = {0, 568,600,570,593,599,564,595,569};
    if(depth >= content.size()){
        return false;
    }
    for(int i=1;i<depth+1;++i){
        if(emb[i] != content[i]){
            return false;
        }
    }
    return true;
}

void SubgraphEnum::initialization(){
    // initialize the order_index
    order_index_ = vector<Vertex>(order_.size(), 0);
    for(int offset=0;offset<order_.size();++offset){
        order_index_[order_[offset]] = offset;
    }

    // initialize the ancestors
    uint32_t query_vertex_count = query_graph_->getVerticesCount();
    ancestors_depth_.resize(query_vertex_count+1);
    for(uint32_t i=0; i<query_vertex_count+1; ++i){
        ancestors_depth_[i].reset();
    }
    for(uint32_t i=1; i<=query_vertex_count; ++i){
        ancestors_depth_[i].set(i);

        uint32_t nbrs_count;
        const Vertex* nbrs = query_graph_->getVertexNeighbors(order_[i], nbrs_count);
        for(uint32_t j=0;j<nbrs_count;++j){
            Vertex n = nbrs[j];
            uint32_t depth = order_index_[n];
            if(depth > i){
                ancestors_depth_[depth].set(depth);
                ancestors_depth_[depth] = ancestors_depth_[i] | ancestors_depth_[depth];
            }
        }
    }

    // initialize the successors and predecessors
    successor_neighbors_in_depth_.clear();
    predecessor_neighbors_in_depth_.clear();
    successor_neighbors_in_depth_.resize(query_vertex_count+1);
    predecessor_neighbors_in_depth_.resize(query_vertex_count+1);
    for(uint32_t i=1; i<=query_vertex_count; ++i){
        uint32_t nbrs_count;
        const Vertex* nbrs = query_graph_->getVertexNeighbors(order_[i], nbrs_count);
        for(uint32_t j=0;j<nbrs_count;++j){
            Vertex n = nbrs[j];
            uint32_t depth = order_index_[n];
            if(depth > i){
                successor_neighbors_in_depth_[i].push_back(depth);
            }else{
                predecessor_neighbors_in_depth_[i].push_back(depth);
            }
        }
    }
    for(uint32_t i=1;i<=query_vertex_count;++i){
        sort(successor_neighbors_in_depth_[i].begin(), successor_neighbors_in_depth_[i].end());
    }

    // initialize the marker marks all the unmatched queries having candidates
    full_descendent_.set();
    bitset<MAX_QUERY_SIZE> touched_descendent;
    for(int x=1;x<query_vertex_count+1;x++){
        Vertex u = order_[x];
        for(auto suc_u_depth : successor_neighbors_in_depth_[x]){
            touched_descendent.set(suc_u_depth);
        }
        touched_descendent.set(x);
        bool full_touched = true;
        for(int d=x;d<=query_vertex_count;++d){
            if(touched_descendent.test(d) == false){
                full_touched = false;
                break;
            }
        }
        if(full_touched == true){
            break;
        }else{
            full_descendent_.reset(x);
        }
    }

    // intialize the parent map
    query_vertex_count_ = query_graph_->getVerticesCount();
    parent_failing_set_map_.clear();
    parent_failing_set_map_.resize(query_vertex_count_+1);
    for(int i=0;i<=query_vertex_count_; ++i){
        parent_failing_set_map_[i].resize(query_vertex_count_+1);
    }
    for(int d=1;d<=query_vertex_count_; ++d){
        for(auto suc_depth : successor_neighbors_in_depth_[d]){
            parent_failing_set_map_[d][suc_depth] = ancestors_depth_[d] | ancestors_depth_[suc_depth];
        }
    }

    // initialize the data for pghole conflict filter
#if PGHOLE_CONFLICT == 1
    check_neighbor_conflict_.reset();
    unmatched_group.clear();
    matched_group.clear();
    unmatched_group.resize(query_vertex_count+1);
    matched_group.resize(query_vertex_count+1);
    bitset<MAX_QUERY_SIZE> covered_depth;
    for(int d=1;d<=query_vertex_count;++d){
        Vertex u = order_[d];
        covered_depth.set(d);
        unordered_map<Label, set<Vertex>> group;
        Label l_u = query_graph_->getVertexLabel(u);
        group.insert({l_u, {}});
        for(auto suc_depth : successor_neighbors_in_depth_[d]){
            covered_depth.set(suc_depth);
            Label l = query_graph_->getVertexLabel(order_[suc_depth]);
            if(group.find(l) == group.end()){
                group.insert({l, {}});
            }
            group[l].insert(suc_depth);
        }
        for(int d_c = d+1; d_c<=query_vertex_count; ++d_c){
            if(covered_depth.test(d_c) == true && query_graph_->getVertexLabel(order_[d_c]) == l_u){
                group[l_u].insert(d_c);
            }
        }
        bool to_check = false;
        for(auto p : group){
            if(p.second.size() > 1){
                vector<Vertex> tmp;
                tmp.assign(p.second.begin(), p.second.end());
                unmatched_group[u].push_back(tmp);
                tmp.clear();
                for(int d_i=1; d_i<=d; ++d_i){
                    if(query_graph_->getVertexLabel(order_[d_i]) == p.first){
                        tmp.push_back(d_i);
                    }
                }
                matched_group[u].push_back(tmp);
                to_check = true;
            }
        }
        if(to_check == true){
            check_neighbor_conflict_.set(d);
        }
    }
#endif


    failing_set_stack_.clear();
    failing_set_stack_.resize(query_vertex_count+1);
    for(int d=1;d<=query_vertex_count;++d){
        bitset<MAX_QUERY_SIZE> tmp;
        tmp.set(d);
        failing_set_stack_[d].push_back(tmp);
    }

    same_label_vertices_.clear();
    same_label_vertices_.resize(query_vertex_count+1);
    set<uint32_t> frontiers;
    for(int d=1;d<=query_vertex_count;++d){
        frontiers.erase(d);
        for(auto suc : successor_neighbors_in_depth_[d]){
            frontiers.insert(suc);
        }
        for(auto f : frontiers){
            Label l = query_graph_->getVertexLabel(order_[f]);
            if(same_label_vertices_[d].find(l) == same_label_vertices_[d].end()){
                same_label_vertices_[d].insert({l, {}});
            }
        }
        for(int d_p=1;d_p<=d;++d_p){
            Label l = query_graph_->getVertexLabel(order_[d_p]);
            if(same_label_vertices_[d].find(l) != same_label_vertices_[d].end()){
                same_label_vertices_[d][l].push_back(d_p);
            }
        }
        for(auto f : frontiers){
            Label l = query_graph_->getVertexLabel(order_[f]);
            same_label_vertices_[d][l].push_back(f);
        }
    }

}


bool SubgraphEnum::find_conflict(uint32_t cur_depth, uint32_t depth){
    if(depth <= cur_depth){
        Vertex v = embedding_depth_[depth];
        if(used_set_.find(v) == used_set_.end()){
            used_set_.insert(v);
            if(reverse_match_.find(v) == reverse_match_.end() || find_conflict(cur_depth, reverse_match_[v])){
                reverse_match_[v] = depth;
                return true;
            }
        }
        conflict_vec_.push_back(depth);
        return false;
        // reverse_match_[embedding_depth_[depth]] = depth;
        // return false;
    }
    int count = candidates_stack_[depth].rbegin()->second;
    Vertex* cans = candidates_stack_[depth].rbegin()->first;
    for(int i=0;i<count;++i){
        Vertex v = cans[i];
        // if(visited_query_depth_[v] > 0){
        //     continue;
        // }
        if(used_set_.find(v) == used_set_.end()){
            used_set_.insert(v);
            if(reverse_match_.find(v) == reverse_match_.end() || find_conflict(cur_depth, reverse_match_[v])){
                reverse_match_[v] = depth;
                return true;
            }
        }
    }
    conflict_vec_.push_back(depth);
    return false;
}

bool SubgraphEnum::check_conflict(uint32_t cur_depth){
    for(auto p : same_label_vertices_[cur_depth]){
        reverse_match_.clear();
        for(auto d : p.second){
            used_set_.clear();
            conflict_vec_.clear();
            if(find_conflict(cur_depth, d) == false){
                reverse_match_[MAX_QUERY_SIZE] = d;
                return false;
            }
        }
    }
    return true;
}



bool validate_correctness(Graph* query_graph, Graph* data_graph, vector<Vertex>& order_index, vector<Vertex>& emb_depth){
    for(Vertex u=0;u<query_graph->getVerticesCount();++u){
        uint32_t nbrs_count;
        const Vertex* nbrs = query_graph->getVertexNeighbors(u, nbrs_count);
        for(int i=0;i<nbrs_count;++i){
            Vertex u_n = nbrs[i];
            Vertex v = emb_depth[order_index[u]];
            Vertex v_n = emb_depth[order_index[u_n]];
            if(data_graph->checkEdgeExistence(v,v_n) == false){
                return false;
            }
        }
    }
    return true;
}

void SubgraphEnum::debug_print(int depth){
    for(int d=depth; d<order_.size();d++){
        if(candidates_stack_[d].size() > 0){
            Vertex* can = candidates_stack_[d].rbegin()->first;
            uint32_t count = candidates_stack_[d].rbegin()->second;
            cout<<d<<"("<<query_graph_->getVertexLabel(order_[d])<<"):";
            for(int i=0;i<count;++i){
                if(visited_query_depth_[can[i]] == 0){
                    cout<<can[i]<<",";
                }else{
                    cout<<can[i]<<"(x),";
                }
            }
            cout<<endl;
        }
    }
}


void MDE(Graph *query_graph, vector<Vertex>& order_, catalog *storage) {
    bool is_isolated[260];
    int query_vertex_number = query_graph->getVerticesCount();
    // order.resize(query_vertex_number*2);
    Vertex* order = new Vertex[query_vertex_number*2];
    bool used[260];
    bool extendable[260];
    for (int i = 0; i < query_vertex_number; ++i) {
        extendable[i] = false;
        is_isolated[i] = false;
        used[i] = false;
    }

    uint32_t* candidates_count = new uint32_t[query_vertex_number];
    for(Vertex u=0;u<query_vertex_number;++u){
        candidates_count[u] = storage->num_candidates_[u];
    }

    int start_vertex = 0;
    for (int i = 1; i < query_vertex_number; ++i) {
        if (query_graph->getVertexDegree(i) > query_graph->getVertexDegree(start_vertex))
            start_vertex = i;
        else if (query_graph->getVertexDegree(i) == query_graph->getVertexDegree(start_vertex)) {
            if (candidates_count[i] < candidates_count[start_vertex])
                start_vertex = i;
        }
    }
    order[0] = start_vertex;
    used[start_vertex] = true;

    Vertex edge_count;
    const Vertex* neighbors = query_graph -> getVertexNeighbors(start_vertex, edge_count);
    for (int i = 0; i < edge_count; ++i)
        extendable[neighbors[i]] = true;
    for (int loop = 1; loop < query_vertex_number; ++loop) {
        int best = -1;
        int scoreA;
        for (int i = 0; i < query_vertex_number; ++i) if (extendable[i]) {
            int ScoreB=0;
            Vertex edge_count;
            const Vertex* neighbors = query_graph -> getVertexNeighbors(i, edge_count);
            for (int p = 0; p < edge_count; ++p) {
                int un = neighbors[p];
                if (!used[un])
                    ++ScoreB;
            }
            if (ScoreB == 0) {
                best = i;
                scoreA = ScoreB;
                break;
            }

            if (best == -1) {
                best = i;
                scoreA = ScoreB;
                continue;
            }
            if (ScoreB > scoreA || (ScoreB == scoreA && candidates_count[i] < candidates_count[best])) {
                best = i;
                scoreA = ScoreB;
            }
        }
        if (scoreA == 0) 
            is_isolated[best] = true;
        assert(extendable[best]);
        extendable[best] = false;
        order[loop] = best;
        used[best] = true;
        Vertex edge_count;
        const Vertex* neighbors = query_graph -> getVertexNeighbors(best, edge_count);
        for (int i = 0; i < edge_count; ++i)
            if (!used[neighbors[i]])
                extendable[neighbors[i]] = true;
    }

    vector<int> tmp(query_vertex_number);
    int l = 0, r = query_vertex_number;
    for (int i = 0; i < query_vertex_number;++i) {
        int v = order[i];
        if (!is_isolated[v])
            tmp[l++] = v;
        else tmp[--r] = v;
    }
    assert(l == r);
    for (int i = 0; i < query_vertex_number;++i) 
        order[i] = tmp[i];
    // printf("iso: %d\n", r);
    for (int i = 0; i < query_vertex_number;++i) 
        order_[i] = order[i];
    delete candidates_count;
}


#if PROACTIVE_CANDIDATE_COMPUTING == 1
void SubgraphEnum::match(Graph* query_graph, string ordering_method, uint64_t count_limit, uint32_t time_limit){
    query_graph_ = query_graph;
    // Execute Preprocessor
    query_time_ = 0;
    intersection_count_ = 0;

#if PRINT_MEM_INFO == 1
    bool stop_thread = false;
    int current_pid = GetCurrentPid();
    thread mem_info_thread(thread_get_mem_info, ref(peak_memory_), ref(stop_thread));
    float starting_memory_cost = GetMemoryUsage(current_pid);
#endif

    pp_ = new preprocessor();
    storage_ = new catalog(query_graph_, data_graph_);
    pp_->execute(query_graph, data_graph_, storage_, true);
    preprocessing_time_ = NANOSECTOSEC(pp_->preprocess_time_);
    query_time_ += preprocessing_time_;

    // Generate Query Plan
    std::vector<std::vector<uint32_t>> spectrum;
    query_plan_generator::generate_query_plan_with_nd(query_graph, storage_, spectrum);
    ordering_time_ = NANOSECTOSEC(query_plan_generator::ordering_time_);

    // Order Adjustment
    auto start = std::chrono::high_resolution_clock::now();
    order_ = spectrum[0];
    MDE(query_graph_, order_, storage_);
    // order_ = {30,46,52,31,5,34,49,22,20,13,37,0,27,38,18,23,25,42,39,24,53,9,32,41,15,8,21,33,6,35,29,55,45,43,12,10,54,40,3,50,36,16,14,26,48,19,11,28,44,51,47,7,4,2,1,17};
    order_.insert(order_.begin(), 0); // padding
    // cout<<"original_order:";
    // for(auto n : order_){
    //     cout<<n<<", ";
    // }
    // cout<<endl;
    // order_adjustment();
    initialization();
    // for(auto n : order_){
    //     cout<<n<<", ";
    // }
    // cout<<endl;
    auto end = std::chrono::high_resolution_clock::now();
    order_adjust_time_ = NANOSECTOSEC(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
    query_time_ += order_adjust_time_;

    // encoding the relations
    start = std::chrono::high_resolution_clock::now();
    TrieEncoder* encoder = new TrieEncoder(storage_, order_, order_index_, query_graph_, data_graph_, successor_neighbors_in_depth_);
    end = std::chrono::high_resolution_clock::now();
    query_time_ += NANOSECTOSEC(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

    matches_.clear();
#if PRINT_LEAF_STATE == 1
    leaf_states_counter_ = vector<uint64_t>(20, 0);
#endif

#if SUCCESSOR_EQUIVALENT_SET == 1
    SuccessorEquivalentCache sccache(order_, order_index_, query_graph_, data_graph_, successor_neighbors_in_depth_, ancestors_depth_, count_limit);
    compressed_ = sccache.compressed_;
    compressed_neq_ = sccache.compressed_neq_;
#endif

    // start enumeration
    // start timer
    stop_ = false;
    register_timer_enum(&stop_);
    start = std::chrono::high_resolution_clock::now();

    query_vertex_count_ = query_graph_->getVerticesCount();
    data_vertex_count_ = data_graph_->getVerticesCount();

    // initializing candidates stack
    history_candidates_ = new CandidatesHistoryStack [query_vertex_count_+1];
    for(int d=1;d<=query_vertex_count_;++d){
        history_candidates_[d].rebuild(query_graph_, data_graph_, successor_neighbors_in_depth_[d], predecessor_neighbors_in_depth_, d);
    }
    candidates_stack_.clear();
    candidates_stack_.resize(query_vertex_count_+1);
    searching_candidates_.clear();
    searching_candidates_.resize(query_vertex_count_+1);

#if CD_EXCLUSION_FILTER==1 || CD_COMPLETE_FILTER == 1
    ContainmentFilter CFilter(query_graph_, data_graph_, successor_neighbors_in_depth_, predecessor_neighbors_in_depth_, order_, order_index_, ancestors_depth_, encoder);
#endif

    // initializing the candidates
    encoder->get_candidates(1, searching_candidates_[1]);
    candidates_stack_[1].push_back({&(searching_candidates_[1][0]), searching_candidates_[1].size()});

    // history_candidates_[1].reorder_searching_init(&(searching_candidates_[1][0]), searching_candidates_[1].size(), encoder);

    embedding_depth_.clear();
    embedding_depth_.resize(query_vertex_count_+1);

    visited_query_depth_ = new Vertex [data_vertex_count_];
    candidates_offset_ = new Vertex [query_vertex_count_+1];
    memset(visited_query_depth_, 0, sizeof(Vertex)*data_vertex_count_);
    memset(candidates_offset_, 0, sizeof(Vertex)*(query_vertex_count_+1));

    int cur_depth = 1;
    candidates_offset_[cur_depth] = 0;
    find_matches_ = new bool [query_vertex_count_+1];
    find_matches_[query_vertex_count_] = true;

    vector<float> candidate_score;
    Vertex u, v;
    // Vertex conflicted_query_vertex;
    uint32_t conflicted_depth;
    state_count_ = 0;
    emb_count_ = 0;

    vector<bitset<MAX_QUERY_SIZE>> conflict_source_depth;
    vector<bitset<MAX_QUERY_SIZE>> failing_set_depth;
    conflict_source_depth.resize(query_vertex_count_+1);
    failing_set_depth.resize(query_vertex_count_+1);
    vector<uint32_t> max_search_depth = vector<uint32_t>(query_vertex_count_+1, 0);

    // bool debug_print=false;
    Vertex* local_candidates = &(searching_candidates_[cur_depth][0]);
    uint32_t local_size = searching_candidates_[cur_depth].size();
    while(true){
        while(
#if CANDIDATE_REORDER == 1
            candidates_offset_[cur_depth]<searching_candidates_[cur_depth].size()
#else
            candidates_offset_[cur_depth]<candidates_stack_[cur_depth].rbegin()->second
#endif
            ){
            if(stop_==true){
                for(int i=1;i<=cur_depth;++i){
                    cout<<embedding_depth_[i]<<" ";
                }
                cout<<endl;
                goto EXIT;
            }

            // vector<Vertex>& current_candidates = searching_candidates_[cur_depth];
#if CANDIDATE_REORDER == 1
            local_candidates = &(searching_candidates_[cur_depth][0]);
            local_size = searching_candidates_[cur_depth].size();
#else
            local_candidates = candidates_stack_[cur_depth].rbegin()->first;
            local_size = candidates_stack_[cur_depth].rbegin()->second;
#endif
            state_count_ ++;
            find_matches_[cur_depth] = false;
            u = order_[cur_depth];
            uint32_t offset = candidates_offset_[cur_depth];
            v = local_candidates[offset];
            embedding_depth_[cur_depth] = v;

            // if(cur_depth<=18){
            //     cout<<"cur:"<<cur_depth<<":"<<v<<":"<<candidates_offset_[cur_depth]<<endl;
            //     for(int i=1;i<=cur_depth;++i){
            //         cout<<embedding_depth_[i]<<" ";
            //     }
            //     cout<<endl;
            //     // candidates_offset_[cur_depth] ++;
            //     // continue;
            // }
            // cout<<cur_depth<<":";
            // for(int i=1;i<=cur_depth;++i){
            //     cout<<embedding_depth_[i]<<" ";
            // }
            // cout<<endl;

            // cout<<"searching:"<<cur_depth<<":"<<u<<":"<<v<<":"<<candidates_offset_[cur_depth]<<endl;
#if FAILING_SET == 1
            failing_set_depth[cur_depth].reset();
#endif
#if SUCCESSOR_EQUIVALENT_SET == 1
                sccache.clear(cur_depth+1);
#endif
#if CD_EXCLUSION_FILTER==1 || CD_COMPLETE_FILTER == 1
            conflict_source_depth[cur_depth].reset();
#endif

            if(cur_depth == query_vertex_count_){
                // max_search_depth[cur_depth-1] = (max_search_depth[cur_depth-1]<max_search_depth[cur_depth]) ? max_search_depth[cur_depth] : max_search_depth[cur_depth-1];
                for(int i=0;i<local_size; ++i){
                    v = local_candidates[i];
                    state_count_ ++;
                    if(visited_query_depth_[v] == 0){
                        find_matches_[cur_depth-1] = true;
                        find_matches_[cur_depth] = true;
                        emb_count_ ++;
                        embedding_depth_[cur_depth] = v;

#if PRINT_RESULT==1
                        matches_.push_back(embedding_depth_);
                        for(Vertex z=1;z<=query_vertex_count_; ++z){
                            cout<<embedding_depth_[z]<<" ";
                        }
                        // cout<<validate_correctness(query_graph_, data_graph_, order_index_, embedding_depth_);
                        cout<<endl;
                        // exit(0);
#endif
#if FAILING_SET==1
                        // belongs to embedding-class
                        failing_set_depth[cur_depth].reset();
#endif
                        if(emb_count_ >= count_limit){
                            goto EXIT;
                        }
                    }else{

                        conflicted_depth = visited_query_depth_[v];
#if FAILING_SET==1
                        // failing_set_depth[cur_depth] = ancestors_depth_[cur_depth] | ancestors_depth_[conflicted_depth];
                        failing_set_depth[cur_depth] = *(failing_set_stack_[cur_depth].rbegin()) | *(failing_set_stack_[conflicted_depth].rbegin());
#endif
#if CD_EXCLUSION_FILTER==1 || CD_COMPLETE_FILTER == 1
                        conflict_source_depth[conflicted_depth].set(cur_depth);
#endif
                    }
                    
                    // updating the failing set of the parent
                    if(find_matches_[cur_depth] == true){
#if FAILING_SET == 1
                        failing_set_depth[cur_depth-1].reset();
#endif
                    }else if(find_matches_[cur_depth-1] == false){
#if FAILING_SET == 1
                        failing_set_depth[cur_depth-1] |= failing_set_depth[cur_depth];
#endif
                    }
                    // cout<<"conflict:"<<cur_depth<<":"<<v<<":"<<visited_query_depth_[v]<<":"<<state_count_<<endl;
                }
                candidates_offset_[cur_depth] = candidates_stack_[cur_depth].rbegin()->second;
            }else{
                // normal enumeration
                candidates_offset_[cur_depth] ++;
                history_candidates_[cur_depth+1].clear();

                // start intersection
                // history_candidates_[cur_depth].append_history_reorded_strict(v, offset, candidates_stack_, failing_set_stack_);
                history_candidates_[cur_depth].append_history(v, encoder, candidates_stack_, failing_set_stack_
#if LOOKAHEAD_OVERLAP == 0
                    , embedding_depth_
#endif
                );
                // intersection_count_ += successor_neighbors_in_depth_[cur_depth].size();

                // applying containment filtering
#if CD_EXCLUSION_FILTER==1 || CD_COMPLETE_FILTER == 1
                if(offset > 0){
#if CD_EXCLUSION_FILTER==1
                    if(CFilter.validate_exclusion_containment(cur_depth, candidates_stack_, v, history_candidates_[cur_depth]) == true){
                        // cout<<"ff2-exclusion:"<<cur_depth<<":"<<u<<":"<<state_count_<<endl;
                        // if(debug_print){
                        // cout<<"ff2-standard:"<<cur_depth<<":"<<u<<endl;
                        // for(int x=1;x<=cur_depth;++x){
                        //     cout<<embedding_depth_[x]<<", ";
                        // }
                        // cout<<endl;
                        // }
#if PRINT_LEAF_STATE==1
                        leaf_states_counter_[CD_STATE] ++;
#endif
                        for(auto suc_depth : successor_neighbors_in_depth_[cur_depth]){
                            candidates_stack_[suc_depth].pop_back();
#if FAILING_SET == 1
                            failing_set_stack_[suc_depth].pop_back();
                            assert(candidates_stack_[suc_depth].size() == failing_set_stack_[suc_depth].size()-1);
#endif                        
                        }
// #if FAILING_SET == 1
//                         cur_sets.failing_set_ = result_fs->failing_set_;
// #endif
//                         cur_sets.conflict_set_depth_ = result_fs->conflict_set_depth_;
                        max_search_depth[cur_depth-1] = (max_search_depth[cur_depth-1]<max_search_depth[cur_depth]) ? max_search_depth[cur_depth] : max_search_depth[cur_depth-1];
                        continue;
                    }
#endif
#if CD_COMPLETE_FILTER == 1
                    if(CFilter.validate_complete_containment(cur_depth, candidates_stack_, v, history_candidates_[cur_depth], embedding_depth_) == true){
                        // cout<<"ff2-complete:"<<cur_depth<<":"<<u<<":"<<state_count_<<endl;
                        // if(debug_print){
                        // cout<<"ff2-extended:"<<cur_depth<<":"<<u<<endl;
                        // for(int x=1;x<=cur_depth;++x){
                        //     cout<<embedding_depth_[x]<<", ";
                        // }
                        // cout<<endl;
                        // }
                        for(auto suc_depth : successor_neighbors_in_depth_[cur_depth]){
                            candidates_stack_[suc_depth].pop_back();
#if FAILING_SET==1
                            failing_set_stack_[suc_depth].pop_back();
                            assert(candidates_stack_[suc_depth].size() == failing_set_stack_[suc_depth].size()-1);
#endif                        
                        }
#if PRINT_LEAF_STATE==1
                        leaf_states_counter_[CD_STATE]++;
#endif
                        // max_search_depth[cur_depth-1] = (max_search_depth[cur_depth-1]<max_search_depth[cur_depth]) ? max_search_depth[cur_depth] : max_search_depth[cur_depth-1];
                        continue;
                    }
#endif
                }
#endif

// #if FAILING_SET==1
                failing_set_depth[cur_depth].reset();
// #endif

#if PGHOLE_EMPTYSET == 1
                bool empty_candidates = false;
                for(auto suc_depth : successor_neighbors_in_depth_[cur_depth]){
                    if(candidates_stack_[suc_depth].rbegin()->second == 0){
                        if(find_matches_[cur_depth-1] == false){
                            failing_set_depth[cur_depth-1] |= *(failing_set_stack_[suc_depth].rbegin());
                            empty_candidates = true;
                            break;
                        }
                    }
                }
                if(empty_candidates == true){
#if PRINT_LEAF_STATE==1
                    leaf_states_counter_[FD_STATE] ++;
#endif
                    for(auto suc_depth : successor_neighbors_in_depth_[cur_depth]){
                        candidates_stack_[suc_depth].pop_back();
#if FAILING_SET==1
                        failing_set_stack_[suc_depth].pop_back();
                        assert(candidates_stack_[suc_depth].size() == failing_set_stack_[suc_depth].size()-1);
#endif                    
                    }
                    continue;
                }
//                 bitset<MAX_QUERY_SIZE> failing_set_for_empty;
//                 bool empty_candidates = false;
//                 for(auto suc_depth : successor_neighbors_in_depth_[cur_depth]){
//                     if(candidates_stack_[suc_depth].rbegin()->second == 0){
//                         if(empty_candidates == false){
//                             failing_set_for_empty = ancestors_depth_[suc_depth];
//                             empty_candidates = true;
//                             break;
//                         }else{
//                             failing_set_for_empty &= ancestors_depth_[suc_depth];
//                         }
//                         // cout<<cur_depth<<":"<<suc_depth<<endl;
//                     }
//                 }
//                 if(empty_candidates == true){
//                     for(auto suc_depth : successor_neighbors_in_depth_[cur_depth]){
//                         candidates_stack_[suc_depth].pop_back();
//                     }
// #if FAILING_SET==1
//                     failing_set_depth[cur_depth] = failing_set_for_empty;
//                     if(find_matches_[cur_depth-1] == false){
//                         failing_set_depth[cur_depth-1] |= failing_set_depth[cur_depth];
//                     }
// #endif
// #if PRINT_LEAF_STATE==1
//                     leaf_states_counter_[FD_STATE] ++;
// #endif
//                     cout<<"emptyset:"<<cur_depth<<":"<<u<<":"<<v<<":"<<state_count_<<endl;
//                     // max_search_depth[cur_depth-1] = (max_search_depth[cur_depth-1]<max_search_depth[cur_depth]) ? max_search_depth[cur_depth] : max_search_depth[cur_depth-1];
//                     continue;
//                 }
#endif


                if(visited_query_depth_[v] > 0){
                    conflicted_depth = visited_query_depth_[v];
                    for(auto suc_depth : successor_neighbors_in_depth_[cur_depth]){
                        candidates_stack_[suc_depth].pop_back();
#if FAILING_SET==1
                        failing_set_stack_[suc_depth].pop_back();
                        assert(candidates_stack_[suc_depth].size() == failing_set_stack_[suc_depth].size()-1);
#endif                    
                    }
#if FAILING_SET==1
                    // belongs to conflict-class
                    if(find_matches_[cur_depth-1] == false){
                        // failing_set_depth[cur_depth-1] |= ancestors_depth_[cur_depth] | ancestors_depth_[conflicted_depth];
                        failing_set_depth[cur_depth-1] |= *(failing_set_stack_[cur_depth].rbegin()) | *(failing_set_stack_[conflicted_depth].rbegin());
                    }
#endif
#if CD_EXCLUSION_FILTER==1 || CD_COMPLETE_FILTER == 1
                    conflict_source_depth[conflicted_depth].set(cur_depth);
#endif               
                    // cout<<"conflict:"<<cur_depth<<":"<<v<<":"<<visited_query_depth_[v]<<":"<<state_count_<<endl;
                    // max_search_depth[cur_depth-1] = (max_search_depth[cur_depth-1]<max_search_depth[cur_depth]) ? max_search_depth[cur_depth] : max_search_depth[cur_depth-1];
                    continue;
                }


#if PGHOLE_CONFLICT == 1
                if(check_conflict(cur_depth) == false){
                    // cout<<"DYM_conflict:";
                    for(auto d : conflict_vec_){
                        // cout<<d<<", ";
                        if(find_matches_[cur_depth-1] == false){
                            failing_set_depth[cur_depth-1] |= *(failing_set_stack_[d].rbegin());
                        }
                        if(d <= cur_depth){
                            conflict_source_depth[d].set(1);
                        }
                    }
                    // cout<<endl;

                    for(auto suc_depth : successor_neighbors_in_depth_[cur_depth]){
                        candidates_stack_[suc_depth].pop_back();
#if FAILING_SET==1
                        failing_set_stack_[suc_depth].pop_back();
                        assert(candidates_stack_[suc_depth].size() == failing_set_stack_[suc_depth].size()-1);
#endif                    
                    }
#if PRINT_LEAF_STATE==1
                    leaf_states_counter_[FD_STATE] ++;
#endif
                    continue;
                }
#endif

#if PGHOLE_CONFLICT == 1
//                 visited_query_depth_[v] = cur_depth;
//                 if(check_neighbor_conflict_.test(cur_depth) == true){
//                     bool whether_conflict = true;
//                     vector<Vertex> local_conflict_depths;
//                     vector<pair<int, int>> conflict_pairs; // for conflict source
//                     for(int y=0; y<unmatched_group[u].size(); ++y){
//                         vector<Vertex>& candidate_depths = unmatched_group[u][y];
//                         vector<Vertex>& matched_depths = matched_group[u][y];
//                         conflict_checking_order.clear();
//                         for(auto depth : candidate_depths){
//                             conflict_checking_order.push_back(depth);
//                             conflict_checking_order.push_back(candidates_stack_[depth].rbegin()->second);
//                         }
//                         // reorder
//                         for(int m=0;m<conflict_checking_order.size();m+=2){
//                             for(int n=m+2;n<conflict_checking_order.size();n+=2){
//                                 if(conflict_checking_order[m+1] > conflict_checking_order[n+1]){
//                                     swap(conflict_checking_order[m], conflict_checking_order[n]);
//                                     swap(conflict_checking_order[m+1], conflict_checking_order[n+1]);
//                                 }
//                             }
//                         }
//                         // start validating
//                         neighbor_candidates_.clear();
//                         conflict_pairs.clear();
//                         int m = 0;
//                         bool is_empty = true;
//                         for(;m<conflict_checking_order.size();m+=2){
//                             Vertex* candidate_set = candidates_stack_[conflict_checking_order[m]].rbegin()->first;
//                             uint32_t candidate_set_size = candidates_stack_[conflict_checking_order[m]].rbegin()->second;
                            
//                             is_empty = true;
//                             for(int x=0;x<candidate_set_size;++x){
//                                 Vertex can = candidate_set[x];
//                                 if(visited_query_depth_[can] == 0){
//                                     neighbor_candidates_.insert(can);
//                                     is_empty = false;
//                                 }else{
//                                     conflict_pairs.push_back({visited_query_depth_[can], conflict_checking_order[m]});
//                                     // cout<<"pgconflit:"<<cur_depth<<":"<<v<<":"<<visited_query_depth_[can]<<"-"<<conflict_checking_order[m]<<endl;
//                                 }
//                             }
// #if PGHOLE_EMPTYSET == 0
//                             if(candidate_set_size == 0){
//                                 is_empty = false;
//                             }
// #endif
//                             if(neighbor_candidates_.size() >= candidate_depths.size()){
//                                 whether_conflict = false;
//                                 break;
//                             }else if(neighbor_candidates_.size() < m/2+1){
//                                 whether_conflict = true;
//                                 break;
//                             }
//                             if(is_empty == true){
//                                 whether_conflict = true;
//                                 break;
//                             }
//                         }
//                         if(whether_conflict == true){
//                             if(is_empty == false){
//                                 for(int i=0;i<=m && i!=conflict_checking_order.size();i+=2){
//                                     // local_conflict_depths.push_back(conflict_checking_order[i]);
//                                     failing_set_depth[cur_depth] |= ancestors_depth_[conflict_checking_order[i]];
//                                     // cout<<"pgconflit:"<<cur_depth<<":"<<v<<":"<<conflict_checking_order[i]<<endl;
//                                 }
//                             }
//                             break;
//                         }
//                     }
//                     if(whether_conflict == true){
//                         cout<<"---------------------pghole_conflict:"<<cur_depth<<":"<<v<<endl;
//                         // for(auto conflict_depth : local_conflict_depths){
//                         //     cur_sets.failing_set_ |= ancestors_depth_[conflict_depth];
//                         // }
//                         for(auto& p : conflict_pairs){
//                             conflict_source_depth[p.first].set(p.second);
//                             failing_set_depth[cur_depth] |= ancestors_depth_[p.first];
//                             failing_set_depth[cur_depth] |= ancestors_depth_[p.second];
//                         }
//                         if(find_matches_[cur_depth-1] == false){
//                             failing_set_depth[cur_depth-1] |= failing_set_depth[cur_depth];
//                         }

//                         for(auto suc_depth : successor_neighbors_in_depth_[cur_depth]){
//                             candidates_stack_[suc_depth].pop_back();
//                             failing_set_stack_[suc_depth].pop_back();
//                         }


// #if PRINT_LEAF_STATE==1
//                         leaf_states_counter_[FD_STATE] ++;
// #endif
//                         visited_query_depth_[v] = 0;
//                         // max_search_depth[cur_depth-1] = (max_search_depth[cur_depth-1]<max_search_depth[cur_depth]) ? max_search_depth[cur_depth] : max_search_depth[cur_depth-1];
//                         continue;
//                     }
//                 }
#endif

#if SUCCESSOR_EQUIVALENT_SET == 1
                if(sccache.validate_cache_count(cur_depth, &(embedding_depth_[0])) == true){
                    find_matches_[cur_depth-1] = true;
#if PRINT_RESULT == 1
                    assert(matches_.size() == emb_count_);
#endif
                    for(auto suc_depth : successor_neighbors_in_depth_[cur_depth]){
                        candidates_stack_[suc_depth].pop_back();
                        failing_set_stack_[suc_depth].pop_back();
                    }
                    visited_query_depth_[v] = 0;
                    continue;
                }
                sccache.start_insert_cache(cur_depth, emb_count_);
#endif
                // start extending
                visited_query_depth_[v] = cur_depth;
                // if(debug(cur_depth, &(embedding_depth_[0]), order_)){
                //     cout<<"search:"<<cur_depth<<":"<<u<<":"<<v<<endl;
                // }
                uint32_t next_depth = cur_depth+1;
                // Vertex next_u = order_[next_depth];
                candidates_offset_[next_depth] = 0;
                cur_depth ++;
#if CANDIDATE_REORDER == 1
                // searching_candidates_[next_depth].clear();
                // uint32_t num_candidates = candidates_stack_[next_depth].rbegin()->second;
                // Vertex* cans = candidates_stack_[next_depth].rbegin()->first;
                // history_candidates_[cur_depth].reorder_searching_candidates(searching_candidates_[next_depth], cans, num_candidates, encoder, candidates_stack_);

                if(next_depth != query_vertex_count_
#if PGHOLE_EMPTYSET == 0
                    && candidates_stack_[next_depth].rbegin()->second>0
#endif
                ){
                    uint32_t num_candidates = candidates_stack_[next_depth].rbegin()->second;
                    Vertex* cans = candidates_stack_[next_depth].rbegin()->first;
                    searching_candidates_[next_depth] = vector<Vertex>(cans, cans+num_candidates);
                    candidate_score.resize(num_candidates);
                    memset(&candidate_score[0], 0, sizeof(float)*num_candidates);
                    vector<Vertex>& candidates = searching_candidates_[next_depth];

                    for(uint32_t i=0;i<num_candidates;++i){
                        for(auto suc_depth : successor_neighbors_in_depth_[next_depth]){
                            uint32_t count;
                            encoder->get_edge_candidate(next_depth, suc_depth, candidates[i], count);
                            candidate_score[i] += count;
                        }
                    }

                    // reorder the candidates
                    for(int i=0;i<num_candidates-1;++i){
                        for(int j=i+1;j<num_candidates;++j){
                            if(candidate_score[i]<candidate_score[j]){
                                swap(candidate_score[i],candidate_score[j]);
                                swap(candidates[i], candidates[j]);
                            }
                        }
                    }
                }else{
                    uint32_t num_candidates = candidates_stack_[next_depth].rbegin()->second;
                    Vertex* cans = candidates_stack_[next_depth].rbegin()->first;
                    searching_candidates_[next_depth] = vector<Vertex>(cans, cans+num_candidates);
                }
#else
                // local_candidates = candidates_stack_[next_depth].rbegin()->first;
                // local_size = candidates_stack_[next_depth].rbegin()->second;
                // uint32_t num_candidates = candidates_stack_[next_depth].rbegin()->second;
                // Vertex* cans = candidates_stack_[next_depth].rbegin()->first;
                // searching_candidates_[next_depth] = vector<Vertex>(cans, cans+num_candidates);
#endif
            }

        }
#if PGHOLE_EMPTYSET==0
        if(candidates_stack_[cur_depth].rbegin()->second == 0){
#if FAILING_SET == 1
            failing_set_depth[cur_depth] = ancestors_depth_[cur_depth];
#endif
            if(find_matches_[cur_depth-1] == false){
#if FAILING_SET==1
                failing_set_depth[cur_depth-1] |= failing_set_depth[cur_depth];
#endif
            }
#if PRINT_LEAF_STATE==1
            leaf_states_counter_[EMPTYSET] ++;
            leaf_states_depth_[EMPTYSET] += cur_depth;
#endif
            // cout<<"emptyset:"<<cur_depth<<":"<<u<<endl;
        }
#endif 
        cur_depth --;
        if(cur_depth == 0){
            break;
        }

        v = embedding_depth_[cur_depth];
        visited_query_depth_[v] = 0;
        // max_search_depth[cur_depth-1] = (max_search_depth[cur_depth-1]<max_search_depth[cur_depth]) ? max_search_depth[cur_depth] : max_search_depth[cur_depth-1];

#if SUCCESSOR_EQUIVALENT_SET == 1
        if(find_matches_[cur_depth] == true){
            bool full = sccache.finish_insert_cache_count(cur_depth, &(embedding_depth_[0]), emb_count_, visited_query_depth_
#if PRINT_RESULT==1
            , matches_
#endif
            );
            if(full == true){
                goto EXIT;
            }
#if PRINT_RESULT==1
            assert(matches_.size() == emb_count_);
#endif
        }
#endif

#if FAILING_SET==1
        find_matches_[cur_depth-1] |= find_matches_[cur_depth];
        // update the failing_set of the parent
        {
            if(find_matches_[cur_depth] == true){
                failing_set_depth[cur_depth].reset();
            }else if(find_matches_[cur_depth-1]==false){
                failing_set_depth[cur_depth-1] |= failing_set_depth[cur_depth];
            }
        }
#endif

#if CD_EXCLUSION_FILTER==1 || CD_COMPLETE_FILTER == 1
        if(find_matches_[cur_depth] == false){
            *(history_candidates_[cur_depth].history_failing_set_.rbegin()) = failing_set_depth[cur_depth];
            // *(history_candidates_[cur_depth].history_conflict_sources_.rbegin()) = conflict_source_depth[cur_depth];
            // *(history_candidates_[cur_depth].validation_marker_.rbegin()) = 1;
            uint32_t cur_offset = candidates_offset_[cur_depth]-1;
            if(conflict_source_depth[cur_depth].any() == false){
                history_candidates_[cur_depth].validate_exclusion_offsets_.push_back(cur_offset);
            }
            history_candidates_[cur_depth].validate_completion_offset_.push_back(cur_offset);
        }
#endif

        for(auto suc_depth : successor_neighbors_in_depth_[cur_depth]){
            candidates_stack_[suc_depth].pop_back();
#if FAILING_SET==1
            failing_set_stack_[suc_depth].pop_back();
            assert(candidates_stack_[suc_depth].size() == failing_set_stack_[suc_depth].size()-1);
#endif
        }

#if FAILING_SET_PRUNING==1
        if(find_matches_[cur_depth] == false && failing_set_depth[cur_depth].test(cur_depth) == false){
#if PRINT_LEAF_STATE==1
            leaf_states_counter_[FAILINGSET_STATE] ++; // += candidates_stack_[cur_depth].rbegin()->second - candidates_offset_[cur_depth];
#endif
            // filtering the redundant siblings
            candidates_offset_[cur_depth] = candidates_stack_[cur_depth].rbegin()->second;
            failing_set_depth[cur_depth-1] = failing_set_depth[cur_depth];
            // debug.set(cur_depth);

            // cout<<"filtering set:"<<cur_depth<<":"<<state_count_<<endl;
        }
#endif
    }
EXIT:
    end = std::chrono::high_resolution_clock::now();
    enumeration_time_ = NANOSECTOSEC(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
    query_time_ += enumeration_time_;

    if(stop_ == false){
        timer_delete(id);
    }

#if PRINT_MEM_INFO == 1
    stop_thread = true;
    mem_info_thread.join();
    // sleep(1);
    float cur_mem = GetMemoryUsage(current_pid);
    peak_memory_ = (peak_memory_ > cur_mem) ? peak_memory_ : cur_mem;
    peak_memory_ -= starting_memory_cost;
#endif

    delete pp_;
    delete storage_;
    delete candidates_offset_;
    delete visited_query_depth_;
    delete find_matches_;

#if CD_EXCLUSION_FILTER==1 || CD_COMPLETE_FILTER == 1
    delete [] history_candidates_;
#endif

    pp_ = NULL;
    storage_ = NULL;
}





#else


void SubgraphEnum::match(Graph* query_graph, string ordering_method, uint64_t count_limit, uint32_t time_limit){
    query_graph_ = query_graph;
    // Execute Preprocessor
    query_time_ = 0;
    intersection_count_ = 0;

#if PRINT_MEM_INFO == 1
    bool stop_thread = false;
    int current_pid = GetCurrentPid();
    thread mem_info_thread(thread_get_mem_info, ref(peak_memory_), ref(stop_thread));
    float starting_memory_cost = GetMemoryUsage(current_pid);
#endif

    pp_ = new preprocessor();
    storage_ = new catalog(query_graph_, data_graph_);
    pp_->execute(query_graph, data_graph_, storage_, true);
    preprocessing_time_ = NANOSECTOSEC(pp_->preprocess_time_);
    query_time_ += preprocessing_time_;

    // Generate Query Plan
    std::vector<std::vector<uint32_t>> spectrum;
    query_plan_generator::generate_query_plan_with_nd(query_graph, storage_, spectrum);
    ordering_time_ = NANOSECTOSEC(query_plan_generator::ordering_time_);

    // Order Adjustment
    auto start = std::chrono::high_resolution_clock::now();
    order_ = spectrum[0];
    MDE(query_graph_, order_, storage_);
    // order_ = {1,0,3,15,5,10,7,11,9,12,14,4,2,6,8,13};
    order_.insert(order_.begin(), 0); // padding
    // cout<<"original_order:";
    // for(auto n : order_){
    //     cout<<n<<", ";
    // }
    // cout<<endl;
    // order_adjustment();
    initialization();
    // for(auto n : order_){
    //     cout<<n<<", ";
    // }
    // cout<<endl;
    auto end = std::chrono::high_resolution_clock::now();
    order_adjust_time_ = NANOSECTOSEC(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
    query_time_ += order_adjust_time_;

    // encoding the relations
    start = std::chrono::high_resolution_clock::now();
    TrieEncoder* encoder = new TrieEncoder(storage_, order_, order_index_, query_graph_, data_graph_, successor_neighbors_in_depth_);
    end = std::chrono::high_resolution_clock::now();
    query_time_ += NANOSECTOSEC(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());


#if PRINT_LEAF_STATE == 1
    leaf_states_counter_ = vector<uint64_t>(20, 0);
#endif

    // start enumeration
    // start timer
    stop_ = false;
    register_timer_enum(&stop_);
    start = std::chrono::high_resolution_clock::now();

    query_vertex_count_ = query_graph_->getVerticesCount();
    data_vertex_count_ = data_graph_->getVerticesCount();


    searching_candidates_.clear();
    searching_candidates_.resize(query_vertex_count_+1);
    // initializing the candidates
    encoder->get_candidates(1, searching_candidates_[1]);

    embedding_depth_.clear();
    embedding_depth_.resize(query_vertex_count_+1);

    visited_query_depth_ = new Vertex [data_vertex_count_];
    candidates_offset_ = new Vertex [query_vertex_count_+1];
    memset(visited_query_depth_, 0, sizeof(Vertex)*data_vertex_count_);
    memset(candidates_offset_, 0, sizeof(Vertex)*(query_vertex_count_+1));

    int cur_depth = 1;
    candidates_offset_[cur_depth] = 0;
    find_matches_ = new bool [query_vertex_count_+1];
    find_matches_[query_vertex_count_] = true;

    vector<float> candidate_score;
    Vertex u, v;
    // Vertex conflicted_query_vertex;
    uint32_t conflicted_depth;
    state_count_ = 0;
    emb_count_ = 0;

    vector<bitset<MAX_QUERY_SIZE>> failing_set_depth;
    failing_set_depth.resize(query_vertex_count_+1);

    // bool debug_print=false;
    Vertex* intersection_buffer = new Vertex[data_graph_->getGraphMaxDegree()];
    uint32_t intersection_buffer_size = 0;
    Vertex* intersection_buffer_tmp = new Vertex[data_graph_->getGraphMaxDegree()];
    uint32_t intersection_buffer_size_tmp = 0;
    while(true){
        while(candidates_offset_[cur_depth]<searching_candidates_[cur_depth].size()){
            if(stop_==true){
                for(int i=1;i<=cur_depth;++i){
                    cout<<embedding_depth_[i]<<" ";
                }
                cout<<endl;
                goto EXIT;
            }

            vector<Vertex>& current_candidates = searching_candidates_[cur_depth];
            state_count_ ++;
            find_matches_[cur_depth] = false;
            u = order_[cur_depth];
            uint32_t offset = candidates_offset_[cur_depth];
            v = current_candidates[offset];
            embedding_depth_[cur_depth] = v;

            // cout<<"searching:"<<cur_depth<<":"<<u<<":"<<v<<":"<<candidates_offset_[cur_depth]<<endl;
#if FAILING_SET == 1
            failing_set_depth[cur_depth].reset();
#endif

            if(cur_depth == query_vertex_count_){
                for(int i=0;i<current_candidates.size(); ++i){
                    v = current_candidates[i];
                    state_count_ ++;
                    if(visited_query_depth_[v] == 0){
                        find_matches_[cur_depth-1] = true;
                        find_matches_[cur_depth] = true;
                        emb_count_ ++;
                        embedding_depth_[cur_depth] = v;

#if PRINT_RESULT==1
                        matches_.push_back(embedding_depth_);
                        for(Vertex z=1;z<=query_vertex_count_; ++z){
                            cout<<embedding_depth_[z]<<" ";
                        }
                        // cout<<validate_correctness(query_graph_, data_graph_, order_index_, embedding_depth_);
                        cout<<endl;
                        // exit(0);
#endif
#if FAILING_SET==1
                        // belongs to embedding-class
                        failing_set_depth[cur_depth].reset();
#endif
                        if(emb_count_ >= count_limit){
                            goto EXIT;
                        }
                    }else{

                        conflicted_depth = visited_query_depth_[v];
#if FAILING_SET==1
                        failing_set_depth[cur_depth] = ancestors_depth_[cur_depth] | ancestors_depth_[conflicted_depth];
#endif
                    }
                    
                    // updating the failing set of the parent
                    if(find_matches_[cur_depth] == true){
#if FAILING_SET == 1
                        failing_set_depth[cur_depth-1].reset();
#endif
                    }else if(find_matches_[cur_depth-1] == false){
#if FAILING_SET == 1
                        failing_set_depth[cur_depth-1] |= failing_set_depth[cur_depth];
#endif
                    }
                    // cout<<"conflict:"<<cur_depth<<":"<<v<<":"<<visited_query_depth_[v]<<":"<<state_count_<<endl;
                }
                candidates_offset_[cur_depth] = searching_candidates_[cur_depth].size();
            }else{
                // normal enumeration
                candidates_offset_[cur_depth] ++;

                
                failing_set_depth[cur_depth].reset();

                if(visited_query_depth_[v] > 0){
                    conflicted_depth = visited_query_depth_[v];
#if FAILING_SET==1
                    // belongs to conflict-class
                    if(find_matches_[cur_depth-1] == false){
                        failing_set_depth[cur_depth-1] |= ancestors_depth_[cur_depth] | ancestors_depth_[conflicted_depth];
                    }
#endif              

                    // cout<<"conflict:"<<cur_depth<<":"<<v<<":"<<visited_query_depth_[v]<<":"<<state_count_<<endl;
                    continue;
                }

                // start extending
                visited_query_depth_[v] = cur_depth;
                // if(debug(cur_depth, &(embedding_depth_[0]), order_)){
                //     // if(cur_depth == 4)
                //     //     debug_print = true;
                //     cout<<"search:"<<cur_depth<<":"<<u<<":"<<v<<endl;
                // }

                // start intersection
                uint32_t next_depth = cur_depth+1;
                intersection_count_ += predecessor_neighbors_in_depth_[next_depth].size();
                // if(predecessor_neighbors_in_depth_[next_depth].size() == 1){
                //     uint32_t pred_depth = predecessor_neighbors_in_depth_[next_depth][0];
                //     uint32_t count;
                //     Vertex* cans = encoder->get_edge_candidate(pred_depth, next_depth, embedding_depth_[pred_depth], count);
                //     searching_candidates_[next_depth] = vector<Vertex>(cans, cans+count);
                // }else{
                //     for(uint32_t x=1;x<predecessor_neighbors_in_depth_[next_depth].size();++x){
                //         uint32_t pred_depth = predecessor_neighbors_in_depth_[next_depth][x];
                //         uint32_t count;
                //         Vertex* cans = encoder->get_edge_candidate(pred_depth, next_depth, embedding_depth_[pred_depth], count);
                //         if(x == 1){
                //             uint32_t prepre_depth = predecessor_neighbors_in_depth_[next_depth][0];
                //             uint32_t pre_count;
                //             Vertex* cans_pre = encoder->get_edge_candidate(prepre_depth, next_depth, embedding_depth_[prepre_depth], pre_count);
                //             ComputeSetIntersection::ComputeCandidates(cans, count, cans_pre, pre_count, intersection_buffer, intersection_buffer_size);
                //         }else{
                //             ComputeSetIntersection::ComputeCandidates(cans, count, intersection_buffer, intersection_buffer_size, intersection_buffer_tmp, intersection_buffer_size_tmp);
                //             swap(intersection_buffer, intersection_buffer_tmp);
                //             swap(intersection_buffer_size, intersection_buffer_size_tmp);
                //         }
                //     }
                //     searching_candidates_[next_depth] = vector<Vertex>(intersection_buffer, intersection_buffer+intersection_buffer_size);
                // }
                // intersection_count_ += predecessor_neighbors_in_depth_[next_depth].size();
                for(uint32_t x=0;x<predecessor_neighbors_in_depth_[next_depth].size();++x){
                    uint32_t pred_depth = predecessor_neighbors_in_depth_[next_depth][x];
                    uint32_t count;
                    Vertex* cans = encoder->get_edge_candidate(pred_depth, next_depth, embedding_depth_[pred_depth], count);
                    if(x == 0){
                        memcpy(intersection_buffer, cans, count*sizeof(Vertex));
                        intersection_buffer_size = count;
                    }else{
                        swap(intersection_buffer, intersection_buffer_tmp);
                        swap(intersection_buffer_size, intersection_buffer_size_tmp);
                        ComputeSetIntersection::ComputeCandidates(cans, count, intersection_buffer_tmp, intersection_buffer_size_tmp, intersection_buffer, intersection_buffer_size);
                    }
                }
                searching_candidates_[next_depth] = vector<Vertex>(intersection_buffer, intersection_buffer+intersection_buffer_size);
                cur_depth++;
                candidates_offset_[cur_depth] = 0;
            }
        }
        if(searching_candidates_[cur_depth].size() == 0){
#if FAILING_SET == 1
            failing_set_depth[cur_depth] = ancestors_depth_[cur_depth];
#endif
            if(find_matches_[cur_depth-1] == false){
#if FAILING_SET==1
                failing_set_depth[cur_depth-1] |= failing_set_depth[cur_depth];
#endif
            }
            // cout<<"emptyset:"<<cur_depth<<":"<<u<<endl;
        }
        cur_depth --;
        if(cur_depth == 0){
            break;
        }

        v = embedding_depth_[cur_depth];
        visited_query_depth_[v] = 0;

#if FAILING_SET==1
        find_matches_[cur_depth-1] |= find_matches_[cur_depth];
        // update the failing_set of the parent
        {
            if(find_matches_[cur_depth] == true){
                failing_set_depth[cur_depth].reset();
            }else if(find_matches_[cur_depth-1]==false){
                failing_set_depth[cur_depth-1] |= failing_set_depth[cur_depth];
            }
        }
#endif

#if FAILING_SET_PRUNING==1
        if(find_matches_[cur_depth] == false && failing_set_depth[cur_depth].test(cur_depth) == false){
#if PRINT_LEAF_STATE==1
            leaf_states_counter_[FAILINGSET_STATE] ++; // += candidates_stack_[cur_depth].rbegin()->second - candidates_offset_[cur_depth];
#endif
            // filtering the redundant siblings
            candidates_offset_[cur_depth] = searching_candidates_[cur_depth].size();
            failing_set_depth[cur_depth-1] = failing_set_depth[cur_depth];
            // cout<<"filtering set:"<<cur_depth<<":"<<state_count_<<endl;
        }
#endif
    }
EXIT:
    end = std::chrono::high_resolution_clock::now();
    enumeration_time_ = NANOSECTOSEC(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
    query_time_ += enumeration_time_;

    if(stop_ == false){
        timer_delete(id);
    }

#if PRINT_MEM_INFO == 1
    stop_thread = true;
    mem_info_thread.join();
    // sleep(1);
    float cur_mem = GetMemoryUsage(current_pid);
    peak_memory_ = (peak_memory_ > cur_mem) ? peak_memory_ : cur_mem;
    peak_memory_ -= starting_memory_cost;
#endif

    delete pp_;
    delete storage_;
    delete candidates_offset_;
    delete visited_query_depth_;
    delete find_matches_;
    delete intersection_buffer;
    delete intersection_buffer_tmp;

    pp_ = NULL;
    storage_ = NULL;
}

#endif




void SubgraphEnum::match_homo(Graph* query_graph, string ordering_method, uint64_t count_limit, uint32_t time_limit){
    query_graph_ = query_graph;
    // Execute Preprocessor
    query_time_ = 0;
    intersection_count_ = 0;

#if PRINT_MEM_INFO == 1
    bool stop_thread = false;
    int current_pid = GetCurrentPid();
    thread mem_info_thread(thread_get_mem_info, ref(peak_memory_), ref(stop_thread));
    float starting_memory_cost = GetMemoryUsage(current_pid);
#endif

    pp_ = new preprocessor();
    storage_ = new catalog(query_graph_, data_graph_);
    pp_->execute(query_graph, data_graph_, storage_, true);
    preprocessing_time_ = NANOSECTOSEC(pp_->preprocess_time_);
    query_time_ += preprocessing_time_;

    // Generate Query Plan
    std::vector<std::vector<uint32_t>> spectrum;
    query_plan_generator::generate_query_plan_with_nd(query_graph, storage_, spectrum);
    ordering_time_ = NANOSECTOSEC(query_plan_generator::ordering_time_);

    // Order Adjustment
    auto start = std::chrono::high_resolution_clock::now();
    order_ = spectrum[0];
    // order_ = {1,0,3,15,5,10,7,11,9,12,14,4,2,6,8,13};
    order_.insert(order_.begin(), 0); // padding
    // cout<<"original_order:";
    // for(auto n : order_){
    //     cout<<n<<", ";
    // }
    // cout<<endl;
    // order_adjustment();
    initialization();
    // for(auto n : order_){
    //     cout<<n<<", ";
    // }
    // cout<<endl;
    auto end = std::chrono::high_resolution_clock::now();
    order_adjust_time_ = NANOSECTOSEC(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
    query_time_ += order_adjust_time_;

    // encoding the relations
    start = std::chrono::high_resolution_clock::now();
    TrieEncoder* encoder = new TrieEncoder(storage_, order_, order_index_, query_graph_, data_graph_, successor_neighbors_in_depth_);
    end = std::chrono::high_resolution_clock::now();
    query_time_ += NANOSECTOSEC(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

    matches_.clear();
#if PRINT_LEAF_STATE == 1
    leaf_states_counter_ = vector<uint64_t>(20, 0);
#endif

#if SUCCESSOR_EQUIVALENT_SET == 1
    SuccessorEquivalentCache sccache(order_, order_index_, query_graph_, data_graph_, successor_neighbors_in_depth_, ancestors_depth_, count_limit);
    compressed_ = sccache.compressed_;
    compressed_neq_ = sccache.compressed_neq_;
#endif

    // start enumeration
    // start timer
    stop_ = false;
    register_timer_enum(&stop_);
    start = std::chrono::high_resolution_clock::now();

    query_vertex_count_ = query_graph_->getVerticesCount();
    data_vertex_count_ = data_graph_->getVerticesCount();

    // initializing candidates stack
    history_candidates_ = new CandidatesHistoryStack [query_vertex_count_+1];
    for(int d=1;d<=query_vertex_count_;++d){
        history_candidates_[d].rebuild(query_graph_, data_graph_, successor_neighbors_in_depth_[d], predecessor_neighbors_in_depth_, d);
    }
    candidates_stack_.clear();
    candidates_stack_.resize(query_vertex_count_+1);
    searching_candidates_.clear();
    searching_candidates_.resize(query_vertex_count_+1);

#if CD_EXCLUSION_FILTER==1 || CD_COMPLETE_FILTER == 1
    ContainmentFilter CFilter(query_graph_, data_graph_, successor_neighbors_in_depth_, predecessor_neighbors_in_depth_, order_, order_index_, ancestors_depth_, encoder);
#endif

    // initializing the candidates
    encoder->get_candidates(1, searching_candidates_[1]);
    candidates_stack_[1].push_back({&(searching_candidates_[1][0]), searching_candidates_[1].size()});

    embedding_depth_.clear();
    embedding_depth_.resize(query_vertex_count_+1);

    visited_query_depth_ = new Vertex [data_vertex_count_];
    candidates_offset_ = new Vertex [query_vertex_count_+1];
    memset(visited_query_depth_, 0, sizeof(Vertex)*data_vertex_count_);
    memset(candidates_offset_, 0, sizeof(Vertex)*(query_vertex_count_+1));

    int cur_depth = 1;
    candidates_offset_[cur_depth] = 0;
    find_matches_ = new bool [query_vertex_count_+1];
    find_matches_[query_vertex_count_] = true;

    vector<float> candidate_score;
    Vertex u, v;
    // Vertex conflicted_query_vertex;
    uint32_t conflicted_depth;
    state_count_ = 0;
    emb_count_ = 0;

    vector<bitset<MAX_QUERY_SIZE>> conflict_source_depth;
    vector<bitset<MAX_QUERY_SIZE>> failing_set_depth;
    conflict_source_depth.resize(query_vertex_count_+1);
    failing_set_depth.resize(query_vertex_count_+1);
    vector<uint32_t> max_search_depth = vector<uint32_t>(query_vertex_count_+1, 0);

    // bool debug_print=false;
    Vertex* local_candidates = &(searching_candidates_[cur_depth][0]);
    uint32_t local_size = searching_candidates_[cur_depth].size();
    while(true){
        while(
#if CANDIDATE_REORDER == 1
            candidates_offset_[cur_depth]<searching_candidates_[cur_depth].size()
#else
            candidates_offset_[cur_depth]<candidates_stack_[cur_depth].rbegin()->second
#endif
            ){
            if(stop_==true){
                for(int i=1;i<=cur_depth;++i){
                    cout<<embedding_depth_[i]<<" ";
                }
                cout<<endl;
                goto EXIT;
            }

            

            // vector<Vertex>& current_candidates = searching_candidates_[cur_depth];
#if CANDIDATE_REORDER == 1
            local_candidates = &(searching_candidates_[cur_depth][0]);
            local_size = searching_candidates_[cur_depth].size();
#else
            local_candidates = candidates_stack_[cur_depth].rbegin()->first;
            local_size = candidates_stack_[cur_depth].rbegin()->second;
#endif

            // if(cur_depth <=10){
            //     cout<<cur_depth<<":"<<candidates_offset_[cur_depth]<<":"<<local_size<<endl;
            // }
            state_count_ ++;
            find_matches_[cur_depth] = false;
            u = order_[cur_depth];
            uint32_t offset = candidates_offset_[cur_depth];
            // v = current_candidates[offset];
            v = local_candidates[offset];
            embedding_depth_[cur_depth] = v;
            
#if FAILING_SET == 1
            failing_set_depth[cur_depth].reset();
#endif

            if(cur_depth == query_vertex_count_){
                for(int i=0;i<local_size; ++i){
                    find_matches_[cur_depth] = true;
                    v = local_candidates[i];
                    state_count_ ++;

                    emb_count_ ++;
#if PRINT_RESULT==1
                        matches_.push_back(embedding_depth_);
                        for(Vertex z=1;z<=query_vertex_count_; ++z){
                            cout<<embedding_depth_[z]<<" ";
                        }
                        // cout<<validate_correctness(query_graph_, data_graph_, order_index_, embedding_depth_);
                        cout<<endl;
                        // exit(0);
#endif
                    failing_set_depth[cur_depth].reset();
                    if(emb_count_ >= count_limit){
                        goto EXIT;
                    }
                    // updating the failing set of the parent

                    // cout<<"conflict:"<<cur_depth<<":"<<v<<":"<<visited_query_depth_[v]<<":"<<state_count_<<endl;
                }
                if(find_matches_[cur_depth] == true){
#if FAILING_SET == 1
                    failing_set_depth[cur_depth-1].reset();
#endif
                }else if(find_matches_[cur_depth-1] == false){
#if FAILING_SET == 1
                    failing_set_depth[cur_depth-1] |= failing_set_depth[cur_depth];
#endif
                }
                find_matches_[cur_depth-1] |= find_matches_[cur_depth];
                candidates_offset_[cur_depth] = candidates_stack_[cur_depth].rbegin()->second;
            }else{
                // normal enumeration
                candidates_offset_[cur_depth] ++;
                history_candidates_[cur_depth+1].clear();

                // start intersection
                history_candidates_[cur_depth].append_history(v, encoder, candidates_stack_, failing_set_stack_
#if LOOKAHEAD_OVERLAP == 0
                    , embedding_depth_
#endif
                );
                // intersection_count_ += successor_neighbors_in_depth_[cur_depth].size();

                // applying containment filtering
#if CD_EXCLUSION_FILTER==1 || CD_COMPLETE_FILTER == 1
                if(offset > 0){
#if CD_EXCLUSION_FILTER==1
                    if(CFilter.validate_containment_for_homo(cur_depth, candidates_stack_, v, history_candidates_[cur_depth]) == true){
                        // cout<<"ff2-exclusion:"<<cur_depth<<":"<<u<<":"<<state_count_<<endl;
                        // if(debug_print){
                        // cout<<"ff2-standard:"<<cur_depth<<":"<<u<<endl;
                        // for(int x=1;x<=cur_depth;++x){
                        //     cout<<embedding_depth_[x]<<", ";
                        // }
                        // cout<<endl;
                        // }
#if PRINT_LEAF_STATE==1
                        leaf_states_counter_[CD_STATE] ++;
#endif
                        for(auto suc_depth : successor_neighbors_in_depth_[cur_depth]){
                            candidates_stack_[suc_depth].pop_back();
                        }
// #if FAILING_SET == 1
//                         cur_sets.failing_set_ = result_fs->failing_set_;
// #endif
//                         cur_sets.conflict_set_depth_ = result_fs->conflict_set_depth_;
                        // max_search_depth[cur_depth-1] = (max_search_depth[cur_depth-1]<max_search_depth[cur_depth]) ? max_search_depth[cur_depth] : max_search_depth[cur_depth-1];
                        continue;
                    }
#endif
#if CD_COMPLETE_FILTER == 1
                    if(CFilter.validate_complete_containment(cur_depth, candidates_stack_, v, history_candidates_[cur_depth], embedding_depth_) == true){
                        // cout<<"ff2-complete:"<<cur_depth<<":"<<u<<":"<<state_count_<<endl;
                        // if(debug_print){
                        // cout<<"ff2-extended:"<<cur_depth<<":"<<u<<endl;
                        // for(int x=1;x<=cur_depth;++x){
                        //     cout<<embedding_depth_[x]<<", ";
                        // }
                        // cout<<endl;
                        // }
                        for(auto suc_depth : successor_neighbors_in_depth_[cur_depth]){
                            candidates_stack_[suc_depth].pop_back();
                        }
#if PRINT_LEAF_STATE==1
                        leaf_states_counter_[CD_STATE]++;
#endif
                        // max_search_depth[cur_depth-1] = (max_search_depth[cur_depth-1]<max_search_depth[cur_depth]) ? max_search_depth[cur_depth] : max_search_depth[cur_depth-1];
                        continue;
                    }
#endif
                }
#endif

// #if FAILING_SET==1
                failing_set_depth[cur_depth].reset();
// #endif

#if PGHOLE_EMPTYSET == 1
                bitset<MAX_QUERY_SIZE> failing_set_for_empty;
                bool empty_candidates = false;
                for(auto suc_depth : successor_neighbors_in_depth_[cur_depth]){
                    if(candidates_stack_[suc_depth].rbegin()->second == 0){
                        if(empty_candidates == false){
                            failing_set_for_empty = ancestors_depth_[suc_depth];
                            empty_candidates = true;
                            break;
                        }else{
                            failing_set_for_empty &= ancestors_depth_[suc_depth];
                        }
                        // cout<<cur_depth<<":"<<suc_depth<<endl;
                    }
                }
                if(empty_candidates == true){
                    for(auto suc_depth : successor_neighbors_in_depth_[cur_depth]){
                        candidates_stack_[suc_depth].pop_back();
                    }
#if FAILING_SET==1
                    failing_set_depth[cur_depth] = failing_set_for_empty;
                    if(find_matches_[cur_depth-1] == false){
                        failing_set_depth[cur_depth-1] |= failing_set_depth[cur_depth];
                    }
#endif
#if PRINT_LEAF_STATE==1
                    leaf_states_counter_[FD_STATE] ++;
#endif
                    // cout<<"emptyset:"<<cur_depth<<":"<<u<<":"<<v<<":"<<state_count_<<endl;
                    // max_search_depth[cur_depth-1] = (max_search_depth[cur_depth-1]<max_search_depth[cur_depth]) ? max_search_depth[cur_depth] : max_search_depth[cur_depth-1];
                    continue;
                }
#endif

                // start extending
                visited_query_depth_[v] = cur_depth;
                // if(debug(cur_depth, &(embedding_depth_[0]), order_)){
                //     // if(cur_depth == 4)
                //     //     debug_print = true;
                //     cout<<"search:"<<cur_depth<<":"<<u<<":"<<v<<endl;
                // }
                uint32_t next_depth = cur_depth+1;
                // Vertex next_u = order_[next_depth];
                candidates_offset_[next_depth] = 0;
                cur_depth ++;
#if CANDIDATE_REORDER == 1
                if(next_depth != query_vertex_count_
#if PGHOLE_EMPTYSET == 0
                    && candidates_stack_[next_depth].rbegin()->second>0
#endif
                ){
                    uint32_t num_candidates = candidates_stack_[next_depth].rbegin()->second;
                    Vertex* cans = candidates_stack_[next_depth].rbegin()->first;
                    searching_candidates_[next_depth] = vector<Vertex>(cans, cans+num_candidates);
                    candidate_score.resize(num_candidates);
                    memset(&candidate_score[0], 0, sizeof(float)*num_candidates);
                    vector<Vertex>& candidates = searching_candidates_[next_depth];

                    for(uint32_t i=0;i<num_candidates;++i){
                        for(auto suc_depth : successor_neighbors_in_depth_[next_depth]){
                            uint32_t count;
                            encoder->get_edge_candidate(next_depth, suc_depth, candidates[i], count);
                            candidate_score[i] += count;
                        }
                    }

                    // reorder the candidates
                    for(int i=0;i<num_candidates-1;++i){
                        for(int j=i+1;j<num_candidates;++j){
                            if(candidate_score[i]<candidate_score[j]){
                                swap(candidate_score[i],candidate_score[j]);
                                swap(candidates[i], candidates[j]);
                            }
                        }
                    }
                }else{
                    uint32_t num_candidates = candidates_stack_[next_depth].rbegin()->second;
                    Vertex* cans = candidates_stack_[next_depth].rbegin()->first;
                    searching_candidates_[next_depth] = vector<Vertex>(cans, cans+num_candidates);
                }
#else
                // local_candidates = candidates_stack_[next_depth].rbegin()->first;
                // local_size = candidates_stack_[next_depth].rbegin()->second;
                // uint32_t num_candidates = candidates_stack_[next_depth].rbegin()->second;
                // Vertex* cans = candidates_stack_[next_depth].rbegin()->first;
                // searching_candidates_[next_depth] = vector<Vertex>(cans, cans+num_candidates);
#endif
            }

        }
#if PGHOLE_EMPTYSET==0
        if(candidates_stack_[cur_depth].rbegin()->second == 0){
#if FAILING_SET == 1
            failing_set_depth[cur_depth] = ancestors_depth_[cur_depth];
#endif
            if(find_matches_[cur_depth-1] == false){
#if FAILING_SET==1
                failing_set_depth[cur_depth-1] |= failing_set_depth[cur_depth];
#endif
            }
#if PRINT_LEAF_STATE==1
            leaf_states_counter_[EMPTYSET] ++;
            leaf_states_depth_[EMPTYSET] += cur_depth;
#endif
            // cout<<"emptyset:"<<cur_depth<<":"<<u<<endl;
        }
#endif 
        cur_depth --;
        if(cur_depth == 0){
            break;
        }

        v = embedding_depth_[cur_depth];
        visited_query_depth_[v] = 0;
        // max_search_depth[cur_depth-1] = (max_search_depth[cur_depth-1]<max_search_depth[cur_depth]) ? max_search_depth[cur_depth] : max_search_depth[cur_depth-1];

#if FAILING_SET==1
        find_matches_[cur_depth-1] |= find_matches_[cur_depth];
        // update the failing_set of the parent
        {
            if(find_matches_[cur_depth] == true){
                failing_set_depth[cur_depth].reset();
            }else if(find_matches_[cur_depth-1]==false){
                failing_set_depth[cur_depth-1] |= failing_set_depth[cur_depth];
            }
        }
#endif

#if CD_EXCLUSION_FILTER==1
        if(find_matches_[cur_depth] == false){
            *(history_candidates_[cur_depth].history_failing_set_.rbegin()) = failing_set_depth[cur_depth];
            // *(history_candidates_[cur_depth].history_conflict_sources_.rbegin()) = conflict_source_depth[cur_depth];
            // *(history_candidates_[cur_depth].validation_marker_.rbegin()) = 1;
            uint32_t cur_offset = candidates_offset_[cur_depth]-1;
            history_candidates_[cur_depth].validate_completion_offset_.push_back(cur_offset);
        }
#endif

        for(auto suc_depth : successor_neighbors_in_depth_[cur_depth]){
            candidates_stack_[suc_depth].pop_back();
        }

#if FAILING_SET_PRUNING==1
        if(find_matches_[cur_depth] == false && failing_set_depth[cur_depth].test(cur_depth) == false){
#if PRINT_LEAF_STATE==1
            leaf_states_counter_[FAILINGSET_STATE] ++; // += candidates_stack_[cur_depth].rbegin()->second - candidates_offset_[cur_depth];
#endif
            // filtering the redundant siblings
            candidates_offset_[cur_depth] = candidates_stack_[cur_depth].rbegin()->second;
            failing_set_depth[cur_depth-1] = failing_set_depth[cur_depth];
            // debug.set(cur_depth);

            // cout<<"filtering set:"<<cur_depth<<":"<<state_count_<<endl;
        }
#endif
    }
EXIT:
    end = std::chrono::high_resolution_clock::now();
    enumeration_time_ = NANOSECTOSEC(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
    query_time_ += enumeration_time_;

    if(stop_ == false){
        timer_delete(id);
    }

#if PRINT_MEM_INFO == 1
    stop_thread = true;
    mem_info_thread.join();
    // sleep(1);
    float cur_mem = GetMemoryUsage(current_pid);
    peak_memory_ = (peak_memory_ > cur_mem) ? peak_memory_ : cur_mem;
    peak_memory_ -= starting_memory_cost;
#endif

    delete pp_;
    delete storage_;
    delete candidates_offset_;
    delete visited_query_depth_;
    delete find_matches_;

#if CD_EXCLUSION_FILTER==1 || CD_COMPLETE_FILTER == 1
    delete [] history_candidates_;
#endif

    pp_ = NULL;
    storage_ = NULL;
}




