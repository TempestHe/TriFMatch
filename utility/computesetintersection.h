#ifndef SUBGRAPHMATCHING_COMPUTE_SET_INTERSECTION_H
#define SUBGRAPHMATCHING_COMPUTE_SET_INTERSECTION_H

#include "../configuration/config.h"
#include <immintrin.h>
#include <x86intrin.h>

/*
 * Because the set intersection is designed for computing common neighbors, the target is uieger.
 */

class ComputeSetIntersection {
public:
    static size_t galloping_cnt_;
    static size_t merge_cnt_;

    static void ComputeCandidates(const Vertex* larray, uint32_t l_count, const Vertex* rarray,
                                  uint32_t r_count, Vertex* cn, uint32_t &cn_count);
    static void ComputeCandidates(const Vertex* larray, uint32_t l_count, const Vertex* rarray,
                                  uint32_t r_count, uint32_t &cn_count);

#if SI == 0
    static void ComputeCNGallopingAVX2(const Vertex* larray, uint32_t l_count,
                                       const Vertex* rarray, uint32_t r_count, Vertex* cn,
                                       uint32_t &cn_count);
    static void ComputeCNGallopingAVX2(const Vertex* larray, uint32_t l_count,
                                       const Vertex* rarray, uint32_t r_count, uint32_t &cn_count);

    static void ComputeCNMergeBasedAVX2(const Vertex* larray, uint32_t l_count, const Vertex* rarray,
                                        uint32_t r_count, Vertex* cn, uint32_t &cn_count);
    static void ComputeCNMergeBasedAVX2(const Vertex* larray, uint32_t l_count, const Vertex* rarray,
                                        uint32_t r_count, uint32_t &cn_count);
    static const uint32_t BinarySearchForGallopingSearchAVX2(const Vertex*  array, uint32_t offset_beg, uint32_t offset_end, uint32_t val);
    static const uint32_t GallopingSearchAVX2(const Vertex*  array, uint32_t offset_beg, uint32_t offset_end, uint32_t val);
#elif SI == 1

    static void ComputeCNGallopingAVX512(const Vertex* larray, const uint32_t l_count,
                                         const Vertex* rarray, const uint32_t r_count, Vertex* cn,
                                         uint32_t &cn_count);
    static void ComputeCNGallopingAVX512(const Vertex* larray, const uint32_t l_count,
                                         const Vertex* rarray, const uint32_t r_count, uint32_t &cn_count);

    static void ComputeCNMergeBasedAVX512(const Vertex* larray, const uint32_t l_count, const Vertex* rarray,
                                          const uint32_t r_count, Vertex* cn, uint32_t &cn_count);
    static void ComputeCNMergeBasedAVX512(const Vertex* larray, const uint32_t l_count, const Vertex* rarray,
                                          const uint32_t r_count, uint32_t &cn_count);

#elif SI == 2

    static void ComputeCNNaiveStdMerge(const Vertex* larray, uint32_t l_count, const Vertex* rarray,
                                       uint32_t r_count, Vertex* cn, uint32_t &cn_count);
    static void ComputeCNNaiveStdMerge(const Vertex* larray, uint32_t l_count, const Vertex* rarray,
                                       uint32_t r_count, uint32_t &cn_count);

    static void ComputeCNGalloping(const Vertex * larray, uint32_t l_count, const Vertex * rarray,
                                   uint32_t r_count, Vertex * cn, ui& cn_count);
    static void ComputeCNGalloping(const Vertex * larray, uint32_t l_count, const Vertex * rarray,
                                   uint32_t r_count, ui& cn_count);
    static const uint32_t GallopingSearch(const Vertex *src, uint32_t begin, uint32_t end, uint32_t target);
    static const uint32_t BinarySearch(const Vertex *src, uint32_t begin, uint32_t end, uint32_t target);

#endif
};


#endif //FSE_COMPUTESETINTERSECTION_H
