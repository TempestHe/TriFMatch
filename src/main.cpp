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
