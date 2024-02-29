// #pragma once
// #include "../graph/graph.h"
// #include "../utility/utils.h"

// struct MatchInfo{
//     uint32_t emb_count;
//     FailingSet fs;
//     vector<Vertex> conflicted_vertices;
// };

// struct HashFunc_Vec{
//     size_t operator() (const vector<Vertex>& key) const{
//         std::hash<int> hasher;
//         size_t seed = 0;
//         for(auto v : key){
//             seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
//         }
//         return seed;
//     }
// };

// class SuccessorEquivalentCache{
// public:
//     bitset<MAX_QUERY_SIZE> mask_;
//     vector<uint32_t> equivalent_set_depth_;
//     unordered_map<vector<Vertex>, MatchInfo, HashFunc_Vec> recoder_emb_count_;
//     uint32_t start_depth_;
//     uint32_t end_depth_;

//     SuccessorEquivalentCache(vector<uint32_t>& equivalent_set_depth);
// #if PRINT_RESULT == 1
//     unordered_map<vector<Vertex>, vector<vector<Vertex>>, HashFunc_Vec> recoder_emb_result_;
// #endif

//     void insert_emb_count(Vertex* emb, FailingSet* fs, uint32_t emb_count);
//     bool validate_emb_count(Vertex* embedding, FailingSet*& fs, uint32_t& emb_count, Vertex* visited_query_vertices, vector<Vertex>& order_index, vector<bitset<MAX_QUERY_SIZE>>& ancestors_depth);
//     void insert_emb_result(Vertex* emb, vector<vector<Vertex>>& matches, uint32_t starting_emb_count, uint32_t emb_count);
//     bool validate_emb_result(Vertex* embedding, vector<vector<Vertex>>& matches);
//     void clear();

// };

// class SuccessorEquivalentSet{
// public:
//     vector<Vertex> order_;
//     vector<uint32_t> order_index_;
//     vector<bitset<MAX_QUERY_SIZE>> ancestors_depth_;

//     vector<SuccessorEquivalentCache*> end_index_;
//     vector<vector<SuccessorEquivalentCache*>> start_index_;

//     vector<uint32_t> emb_count_recorder_;
//     uint32_t query_vertex_count_;

//     SuccessorEquivalentSet(Graph* query_graph, vector<Vertex>& order, vector<Vertex>& order_index, vector<bitset<MAX_QUERY_SIZE>>& ancestors_depth);

//     void clear_recording(int depth);
//     bool validate_recording(Vertex* embedding, FailingSet*& fs, int depth, Vertex u, uint32_t& emb_count, Vertex* visited_query_vertices);
//     void start_recording(uint32_t emb_count, int depth, Vertex u);
//     void stop_recording(Vertex* embedding, FailingSet* fs, uint32_t emb_count, int depth);

// #if PRINT_RESULT == 1
//     void recording_matches(Vertex* embedding, uint32_t emb_count, int depth, vector<vector<Vertex>>& matches);
//     void validate_recording_matches(Vertex* embedding, int depth, Vertex u, vector<vector<Vertex>>& matches);
// #endif

//     ~SuccessorEquivalentSet();
// };

// struct SuccessorEquivalentPair{
//     Vertex parent;
//     Vertex u;
// };

