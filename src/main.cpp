#include <iostream>
#include <sys/types.h>
#include <dirent.h>
#include <getopt.h>
#include <unistd.h>
#include <unordered_set>
#include <limits>
#include "../graph/graph.h"
#include "preprocessor.h"
#include "query_plan_generator.h"
#include "subgraph_enumeration.h"
// #include "../utility/utils.h"
// #include "../model/model.h"
// #include "../graph_decomposition/decomposition.h"

using namespace std;

void getFiles(string path, vector<string>& filenames){
    DIR *pDir;
    struct dirent* ptr;
    if(!(pDir = opendir(path.c_str()))){
        cout<<"Folder does not exist"<<endl;
        return;
    }
    while((ptr = readdir(pDir))!=0){
        if(strcmp(ptr->d_name, ".")!=0 && strcmp(ptr->d_name, "..") != 0){
            filenames.push_back(path+"/"+ptr->d_name);
        }
    }
    closedir(pDir);
}


struct Param{
    string query_path;
    string data_file;
    uint64_t num;
    uint32_t offset;
};

static struct Param parsed_input_para;

static const struct option long_options[] = {
    {"query", required_argument, NULL, 'q'},
    {"data", required_argument, NULL, 'd'},
    {"num", required_argument, NULL, 'n'},
    {"offset", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, '?'},
};

void parse_args(int argc, char** argv){
    int opt;
    int options_index=0;
    string suffix;
    parsed_input_para.num = std::numeric_limits<uint64_t>::max();
    parsed_input_para.offset = 0;
    while((opt=getopt_long_only(argc, argv, "q:d:n:?", long_options, &options_index)) != -1){
        switch (opt)
        {
        case 0:
            break;
        case 'q':
            parsed_input_para.query_path = string(optarg);
            break;
        case 'd':
            parsed_input_para.data_file = string(optarg);
            break;
        case 'o':
            parsed_input_para.offset = atoi(optarg);
            break;
        case 'n':
            if(string(optarg).compare(string("MAX")) != 0){
                parsed_input_para.num = atoi(optarg);
            }
            break;
        case '?':
            cout<<"------------------ args list ------------------------"<<endl;
            cout<<"--query\tpath of the query graph"<<endl;
            cout<<"--data\tpath of the data graph"<<endl;
            cout<<"--num\tnumber of results to be found"<<endl;
            break;
        default:
            break;
        }
    }
}


void print_input_args(){
    cout<<"input_data_graph:"<<parsed_input_para.data_file<<endl;
    cout<<"input_query: "<<parsed_input_para.query_path<<endl;
    cout<<"maximun number of results intend to find: "<<parsed_input_para.num<<endl;
}

int main(int argc, char** argv){
    parse_args(argc, argv);

    cout<<"Loading Graphs"<<endl;
    auto start = std::chrono::high_resolution_clock::now();
    Graph* data_graph = new Graph(true);
    Graph* query_graph = new Graph(true);
    data_graph->loadGraphFromFile(parsed_input_para.data_file);
    query_graph->loadGraphFromFile(parsed_input_para.query_path);
    query_graph->buildCoreTable();
    auto end = std::chrono::high_resolution_clock::now();

    cout<<"Graph Loading Time (s):"<<NANOSECTOSEC(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count())<<endl;
    
    // std::cout << "Query Graph Meta Information" << std::endl;
    // query_graph->printGraphMetaData();
    // std::cout << "-----" << std::endl;
    // std::cout << "Data Graph Meta Information" << std::endl;
    // data_graph->printGraphMetaData();

    cout<<"Start Querying"<<endl;
    // long result_count;
    SubgraphEnum subgraph_enum(data_graph);
    double enumeration_time, preprocessing_time, ordering_time;
    long long state_count=0;
    subgraph_enum.match(query_graph, string("nd"), parsed_input_para.num, 300);
    enumeration_time = subgraph_enum.enumeration_time_;
    preprocessing_time = subgraph_enum.preprocessing_time_;
    ordering_time = subgraph_enum.ordering_time_;
    state_count = subgraph_enum.state_count_;
    long result_count = subgraph_enum.emb_count_;

    cout<<"========== RESULT INFO ============"<<endl;
    cout<<"Results #:"<<result_count<<endl;
    cout<<"Total Query Time (s):\t"<<subgraph_enum.query_time_<<endl;
    cout<<"Enumeration Time (s):\t"<<subgraph_enum.enumeration_time_<<endl;
    cout<<"Preprocessing Time (s):\t"<<subgraph_enum.preprocessing_time_<<endl;
    cout<<"Ordering Time (s):\t"<<subgraph_enum.ordering_time_<<endl;
    cout<<"State Count:\t"<<subgraph_enum.state_count_<<endl;
#if PRINT_MEM_INFO == 1
    cout<<"Peak Memory (MB):\t"<<subgraph_enum.peak_memory_<<endl;
#endif

#if PRINT_LEAF_STATE == 1
    cout<<"============ LEAF STATE DISTRIBUTION ============"<<endl;
    cout<<"State # Filtered by Containment-Driven filtering:"<<subgraph_enum.leaf_states_counter_[CD_STATE]<<endl;
    cout<<"State # Filtered by Failure-Driven filtering:"<<subgraph_enum.leaf_states_counter_[FD_STATE]<<endl;
    cout<<"State # Filtered by Failing set pruning:"<<subgraph_enum.leaf_states_counter_[FAILINGSET_STATE]<<endl;
#endif

    return 0;
}




// int main(int argc, char** argv){
//     parse_args(argc, argv);
//     Graph* data_graph = new Graph(true);
    
//     data_graph->loadGraphFromFile(parsed_input_para.data_file);

//     vector<string> query_files;
//     getFiles(parsed_input_para.query_path, query_files);
    
//     SubgraphEnum subgraph_enum(data_graph);
//     // query_files = {"../../../dataset/wordnet/query_graph/query_sparse_16_108.graph"};
//     // query_files = {"../../../dataset/extended_queries/youtube/sparse/query_64_123"};
//     query_files = {"/home/yixin/jiezhonghe/project/TriFMatch/sss/result/query_64_223"};
    
//     // query_files = {"../../../dataset/yeast/query_graph/query_dense_8_3.graph"};
//     // "../../../dataset/dblp/extended/sparse/query_64_10"
//     // };
//     // query_files.clear();
//     // vector<string> files = {
//     // "query_64_35","query_48_199","query_64_108","query_64_179","query_48_127","query_64_37",
//     // "query_64_26","query_64_75","query_40_143","query_48_179","query_56_112","query_40_31",
//     // "query_64_42","query_64_193","query_64_161","query_48_39","query_64_150","query_56_88",
//     // "query_64_175","query_56_46","query_40_185","query_40_99","query_48_98","query_64_16",
//     // "query_64_101","query_48_87","query_56_1","query_56_12","query_64_11","query_56_33",
//     // "query_40_58","query_56_178","query_64_147","query_48_90","query_56_54","query_56_123",
//     // "query_48_125","query_48_124","query_56_30","query_64_173","query_56_20","query_56_172",
//     // "query_64_141","query_64_122","query_48_115","query_64_52","query_56_192","query_64_92",
//     // "query_64_10","query_64_5","query_56_117","query_64_6","query_48_41","query_56_18",
//     // "query_48_157","query_56_25","query_48_73","query_64_143","query_64_109","query_64_58",
//     // "query_56_28","query_48_24","query_64_40","query_64_48","query_48_64","query_56_44",
//     // "query_48_49","query_64_31","query_48_133","query_48_4","query_56_125","query_64_146",
//     // "query_48_42","query_48_58","query_64_67","query_48_19","query_64_144","query_48_126",
//     // "query_64_83","query_40_16","query_64_80"
//     // };
//     // vector<string> files = {
//     // "query_48_176","query_48_46","query_56_11","query_64_87","query_64_25","query_56_189",
//     // "query_48_7","query_48_47","query_56_45","query_56_66","query_56_129","query_64_43",
//     // "query_64_103","query_48_17","query_64_57","query_64_100","query_56_65","query_64_116",
//     // "query_64_198","query_64_110","query_56_156","query_64_192","query_64_164","query_64_107",
//     // "query_56_165","query_64_138","query_40_105","query_64_77","query_48_77","query_56_91",
//     // "query_56_119","query_64_165","query_56_113","query_64_64","query_64_86","query_56_35",
//     // "query_56_160","query_48_10","query_56_154","query_56_53","query_48_170","query_56_32",
//     // "query_56_198","query_64_9","query_64_28","query_56_140","query_56_196","query_56_151",
//     // "query_64_105","query_64_20","query_64_69","query_64_0","query_56_106","query_64_88",
//     // "query_48_123","query_64_71","query_64_181","query_56_197","query_64_53","query_64_171",
//     // "query_40_98","query_64_104","query_64_38","query_64_21","query_56_94","query_64_51",
//     // "query_56_5","query_56_37","query_64_169","query_56_146","query_56_13","query_64_183",
//     // "query_64_135","query_48_135","query_56_103","query_64_170","query_64_159","query_48_74",
//     // "query_56_163","query_64_156","query_64_79","query_64_7","query_64_186","query_64_13",
//     // "query_64_91","query_56_29","query_64_2","query_64_114","query_64_24","query_56_48",
//     // "query_56_134","query_64_158","query_56_14","query_56_50","query_64_68","query_64_180",
//     // "query_48_172","query_56_104","query_64_59","query_56_128","query_64_151","query_64_140",
//     // "query_64_191","query_56_133","query_64_47","query_48_13","query_64_27","query_56_21",
//     // "query_64_8","query_56_181","query_48_105","query_56_83","query_64_4","query_56_39",
//     // "query_64_166","query_48_185","query_64_39","query_56_155","query_64_177","query_56_141",
//     // "query_56_142","query_48_150","query_56_43","query_48_83","query_64_137","query_56_60"
//     // };
//     // for(auto q : files){
//     //     query_files.push_back(string("../../../dataset/dblp/extended/dense/")+q);
//     // }

//     int file_id = 0;
//     // uint32_t compressed_total = 0;
//     // uint32_t compressed_neq_total = 0;
//     // uint32_t query_vertices_total = 0;
//     for(auto file : query_files){
//         Graph* query_graph = new Graph(true);
//         query_graph->loadGraphFromFile(file);
//         query_graph->buildCoreTable();
//         if(file_id<parsed_input_para.offset){
//             file_id ++;
//             continue;
//         }
//         // if(file_id!=1038){
//         //     file_id ++;
//         //     continue;
//         // }
//         // cout<<file<<endl;
//         double enumeration_time, preprocessing_time, ordering_time;
//         long long state_count=0;
// #if HOMOMORPHISM == 0
//         subgraph_enum.match(query_graph, string("nd"), parsed_input_para.num, 300);
// #else
//         subgraph_enum.match_homo(query_graph, string("nd"), parsed_input_para.num, 300);
// #endif
//         enumeration_time = subgraph_enum.enumeration_time_;
//         preprocessing_time = subgraph_enum.preprocessing_time_;
//         ordering_time = subgraph_enum.ordering_time_;
//         state_count = subgraph_enum.state_count_;
//         long result_count = subgraph_enum.emb_count_;
//         cout<<file<<":"<<file_id<<": results:"<<result_count<<" query_time:"<<subgraph_enum.query_time_<<" enumeration_time:"<<enumeration_time<<" preprocessing_time:"<<preprocessing_time<<" ordering_time:"<<ordering_time<<" order_adjust_time:"<<subgraph_enum.order_adjust_time_<<" state_count:"<<state_count
// #if SUCCESSOR_EQUIVALENT_SET == 1
//         <<" compressed:"<<subgraph_enum.compressed_
//         <<" compressed_neq:"<<subgraph_enum.compressed_neq_
// #endif
// #if PRINT_MEM_INFO == 1
//         <<" peak_memory:"<<subgraph_enum.peak_memory_
// #endif  
// #if PRINT_LEAF_STATE == 1
//         <<" SUBTREE_REDUCTION:"<<subgraph_enum.leaf_states_counter_[CD_STATE]
//         <<" PGHOLE_FILTERING:"<<subgraph_enum.leaf_states_counter_[FD_STATE]
//         <<" FAILING_SETS:"<<subgraph_enum.leaf_states_counter_[FAILINGSET_STATE]
// #endif
//         <<" intersection_count:"<<subgraph_enum.intersection_count_
//         <<endl;
//         // if(file_id == 6){
//         //     exit(0);
//         // }
//         file_id++;
//         delete query_graph;
//     }
//     // cout<<"compressed_total:"<<compressed_total<<" compressed_neq_total:"<<compressed_neq_total<<" query_vertices_total:"<<query_vertices_total<<endl;
//     // cout<<compressed_total/(float)query_vertices_total<<":"<<compressed_neq_total/(float)query_vertices_total<<endl;
    
//     // std::cout << "Query Graph Meta Information" << std::endl;
//     // query_graph->printGraphMetaData();
//     // std::cout << "-----" << std::endl;
//     // std::cout << "Data Graph Meta Information" << std::endl;
//     // data_graph->printGraphMetaData();

//     // long result_count;
    
//     return 0;
// }
