#pragma once
#include <unordered_map>
#include <iostream>
#include <vector>
#include "../utility/sparsepp/spp.h"
#include "../configuration/config.h"


using spp::sparse_hash_map;
class Graph {
private:
    bool enable_label_offset_;

    uint32_t vertices_count_;
    uint32_t edges_count_;
    uint32_t labels_count_;
    uint32_t max_degree_;
    uint32_t max_label_frequency_;

    uint32_t* offsets_;
    Vertex * neighbors_;
    Label* labels_;
    uint32_t* reverse_index_offsets_;
    uint32_t* reverse_index_;

    int* core_table_;
    uint32_t core_length_;

    std::unordered_map<Label, uint32_t> labels_frequency_;
    sparse_hash_map<uint64_t, std::vector<edge>* >* edge_index_;

    uint32_t* labels_offsets_;
    sparse_hash_map<uint32_t, uint32_t>* nlf_;

private:
    void BuildReverseIndex();

    void BuildNLF();
    void BuildLabelOffset();

public:
    Graph(const bool enable_label_offset) {
        enable_label_offset_ = enable_label_offset;

        vertices_count_ = 0;
        edges_count_ = 0;
        labels_count_ = 0;
        max_degree_ = 0;
        max_label_frequency_ = 0;
        core_length_ = 0;

        offsets_ = NULL;
        neighbors_ = NULL;
        labels_ = NULL;
        reverse_index_offsets_ = NULL;
        reverse_index_ = NULL;
        core_table_ = NULL;
        labels_frequency_.clear();
        edge_index_ = NULL;
        labels_offsets_ = NULL;
        nlf_ = NULL;
    }

    ~Graph() {
        delete[] offsets_;
        delete[] neighbors_;
        delete[] labels_;
        delete[] reverse_index_offsets_;
        delete[] reverse_index_;
        delete[] core_table_;
        delete edge_index_;
        delete[] labels_offsets_;
        delete[] nlf_;
    }

public:
    void loadGraphFromFile(const std::string& file_path);
    void loadGraphFromFileCompressed(const std::string& degree_path, const std::string& edge_path,
                                     const std::string& label_path);
    void storeComparessedGraph(const std::string& degree_path, const std::string& edge_path,
                               const std::string& label_path);
    void printGraphMetaData();
public:
    const uint32_t getLabelsCount() const {
        return labels_count_;
    }

    const uint32_t getVerticesCount() const {
        return vertices_count_;
    }

    const uint32_t getEdgesCount() const {
        return edges_count_;
    }

    const uint32_t getGraphMaxDegree() const {
        return max_degree_;
    }

    const uint32_t getGraphMaxLabelFrequency() const {
        return max_label_frequency_;
    }

    const uint32_t getVertexDegree(const Vertex id) const {
        return offsets_[id + 1] - offsets_[id];
    }

    const uint32_t getLabelsFrequency(const Label label) const {
        return labels_frequency_.find(label) == labels_frequency_.end() ? 0 : labels_frequency_.at(label);
    }

    const uint32_t getCoreValue(const Vertex id) const {
        return core_table_[id];
    }

    const uint32_t get2CoreSize() const {
        return core_length_;
    }
    const Label getVertexLabel(const Vertex id) const {
        return labels_[id];
    }

    const Vertex * getVertexNeighbors(const Vertex id, uint32_t& count) const {
        count = offsets_[id + 1] - offsets_[id];
        return neighbors_ + offsets_[id];
    }

    const sparse_hash_map<uint64_t, std::vector<edge>*>* getEdgeIndex() const {
        return edge_index_;
    }

    const uint32_t * getVerticesByLabel(const Label id, uint32_t& count) const {
        count = reverse_index_offsets_[id + 1] - reverse_index_offsets_[id];
        return reverse_index_ + reverse_index_offsets_[id];
    }

    const Vertex * getNeighborsByLabel(const Vertex id, const Label label, uint32_t& count) const {
        uint32_t offset = id * labels_count_ + label;
        count = labels_offsets_[offset + 1] - labels_offsets_[offset];
        return neighbors_ + labels_offsets_[offset];
    }

    const sparse_hash_map<uint32_t, uint32_t>* getVertexNLF(const Vertex id) const {
        return nlf_ + id;
    }

    bool checkEdgeExistence(const Vertex u, const Vertex v, const Label u_label) const {
        uint32_t count = 0;
        const Vertex* neighbors = getNeighborsByLabel(v, u_label, count);
        int begin = 0;
        int end = count - 1;
        while (begin <= end) {
            int mid = begin + ((end - begin) >> 1);
            if (neighbors[mid] == u) {
                return true;
            }
            else if (neighbors[mid] > u)
                end = mid - 1;
            else
                begin = mid + 1;
        }

        return false;
    }

    bool checkEdgeExistence(Vertex u, Vertex v) const {
        if (getVertexDegree(u) < getVertexDegree(v)) {
            std::swap(u, v);
        }
        uint32_t count = 0;
        const Vertex* neighbors =  getVertexNeighbors(v, count);

        int begin = 0;
        int end = count - 1;
        while (begin <= end) {
            int mid = begin + ((end - begin) >> 1);
            if (neighbors[mid] == u) {
                return true;
            }
            else if (neighbors[mid] > u)
                end = mid - 1;
            else
                begin = mid + 1;
        }

        return false;
    }

    void buildCoreTable();

    void buildEdgeIndex();
};

