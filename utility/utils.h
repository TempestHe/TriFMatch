#pragma once
#include <string>
#include <vector>
#include <unordered_set>
#include <fstream>
#include <iostream>
#include <unistd.h>

#include "../configuration/config.h"

double get_time(timeval& st, timeval& et);

uint64_t compress_vertex_pair(Vertex conflicted_query_vertex, Vertex conflict_origin);
uint64_t compress_vertex_pair_perm(Vertex source, Vertex target);


class DynamicBitmap{
private:
    uint32_t bit_size_;
    uint8_t char_size_;
    char* content_;
public:
    DynamicBitmap();

    DynamicBitmap(uint32_t size);

    DynamicBitmap(const DynamicBitmap& bitmap);

    DynamicBitmap& operator=(const DynamicBitmap& bitmap);

    DynamicBitmap& operator|=(DynamicBitmap& fs);

    DynamicBitmap& operator&=(DynamicBitmap& fs);

    void resize(uint32_t size);

    void set();

    void reset();

    void set(uint32_t i);

    void reset(uint32_t i);

    bool test(uint32_t i);

    void swap_content(DynamicBitmap& bitmap);

    ~DynamicBitmap();
};

class FailingSet{
public:
    bitset<MAX_QUERY_SIZE> failing_set_;

    FailingSet(const FailingSet& fs);
    FailingSet& operator=(FailingSet& fs);

    FailingSet();
};

void intersect_candidates(vector<Vertex>& vec1, vector<Vertex>& vec2, vector<Vertex>& result);
bool validate_vector_containment(vector<Vertex>& vec1, vector<Vertex>& vec2); // true if vec1 contains vec2
bool validate_vector_element_existence(vector<Vertex>& vec1, Vertex v);