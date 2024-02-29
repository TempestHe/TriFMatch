#include "utils.h"
#include <string.h>

double get_time(timeval& st, timeval& et){
    return (double)(et.tv_sec-st.tv_sec+(double)(et.tv_usec-st.tv_usec)/1000000);
}

uint64_t compress_vertex_pair(Vertex conflicted_query_vertex, Vertex conflict_origin){
    return (uint64_t)conflicted_query_vertex<<32|conflict_origin;
}

uint64_t compress_vertex_pair_perm(Vertex source, Vertex target){
    if(source > target){
        return (uint64_t)target<<32|source;
    }else{
        return (uint64_t)source<<32|target;
    }
}


DynamicBitmap::DynamicBitmap(){
    content_ = NULL;
}

DynamicBitmap::DynamicBitmap(uint32_t size){
    char_size_ = sizeof(char);
    bit_size_ = size/char_size_+1;
    content_ = new char [bit_size_];
    memset(content_, 0, char_size_*bit_size_);
}

void DynamicBitmap::resize(uint32_t size){
    char_size_ = sizeof(char);
    bit_size_ = size/char_size_+1;
    content_ = new char [bit_size_];
    memset(content_, 0, char_size_*bit_size_);
}

DynamicBitmap::DynamicBitmap(const DynamicBitmap& bitmap){
    bit_size_ = bitmap.bit_size_;
    content_ = new char [bit_size_];
    char_size_ = sizeof(char);
    memcpy(content_, bitmap.content_, char_size_*bit_size_);
}

DynamicBitmap& DynamicBitmap::operator=(const DynamicBitmap& bitmap){
    bit_size_ = bitmap.bit_size_;
    content_ = new char [bit_size_];
    char_size_ = sizeof(char);
    memcpy(content_, bitmap.content_, char_size_*bit_size_);
    return *this;
}

void DynamicBitmap::set(){
    memset(content_, 1, char_size_*bit_size_);
}

void DynamicBitmap::reset(){
    memset(content_, 0, char_size_*bit_size_);
}

void DynamicBitmap::set(uint32_t i){
    content_[i/char_size_] |= 1<<(i%char_size_);
}

void DynamicBitmap::reset(uint32_t i){
    content_[i/char_size_] &= ~(1<<(i%char_size_));
}

bool DynamicBitmap::test(uint32_t i){
    if(content_[i/char_size_] & 1<<(i%char_size_)>0){
        return true;
    }
    return false;
}

DynamicBitmap& DynamicBitmap::operator|=(DynamicBitmap& fs){
    for(int i=0;i<char_size_;++i){
        content_[i] |= fs.content_[i];
    }
    return *this;
}

DynamicBitmap& DynamicBitmap::operator&=(DynamicBitmap& fs){
    for(int i=0;i<char_size_;++i){
        content_[i] &= fs.content_[i];
    }
    return *this;
}

void DynamicBitmap::swap_content(DynamicBitmap& bitmap){
    swap(bit_size_, bitmap.bit_size_);
    swap(char_size_, bitmap.char_size_);
    swap(content_, bitmap.content_);
}

DynamicBitmap::~DynamicBitmap(){
    if(content_ != NULL)
        delete [] content_;
}



FailingSet::FailingSet(){
    failing_set_.reset();
}

FailingSet::FailingSet(const FailingSet& fs){
    failing_set_ = fs.failing_set_;
}

FailingSet& FailingSet::operator=(FailingSet& fs){
    failing_set_ = fs.failing_set_;
    return *this;
}

void intersect_candidates(vector<Vertex>& vec1, vector<Vertex>& vec2, vector<Vertex>& result){
    for(int i=0,j=0; i<vec1.size() && j<vec2.size();){
        if(vec1[i]==vec2[j]){
            result.push_back(vec1[i]);
            ++i;++j;
        }else if(vec1[i]>vec2[j]){
            ++j;
        }else{
            ++i;
        }
    }
}

bool validate_vector_containment(vector<Vertex>& vec1, vector<Vertex>& vec2){
    if(vec1.size()<vec2.size()){
        return false;
    }
    int i=0,j=0;
    for(; i<vec1.size() && j<vec2.size();){
        if(vec1[i]==vec2[j]){
            ++i;++j;
        }else if(vec1[i]>vec2[j]){
            ++j;
            return false;
        }else{
            ++i;
        }
    }
    if(j<vec2.size()){
        return false;
    }
    return true;
}

bool validate_vector_element_existence(vector<Vertex>& vec1, Vertex v){
    for(int i=0;i<vec1.size();++i){
        if(vec1[i] == v){
            return true;
        }else if(vec1[i] > v){
            return false;
        }
    }
    return false;
}