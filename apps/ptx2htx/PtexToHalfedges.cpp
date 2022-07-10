#include <stdint.h>
#include <vector>

#include "Ptexture.h"
#include "CatmullClark.h"
#include "PtexMesh.h"

#ifndef CC_ASSERT
#    include <assert.h>
#    define CC_ASSERT(x) assert(x)
#endif

#ifndef CC_LOG
#    include <stdio.h>
#    define CC_LOG(format, ...) do { fprintf(stdout, format "\n", ##__VA_ARGS__); fflush(stdout); } while(0)
#endif

#ifndef CC_MALLOC
#    include <stdlib.h>
#    define CC_MALLOC(x) (malloc(x))
#    define CC_FREE(x) (free(x))
#else
#    ifndef CC_FREE
#        error CC_MALLOC defined without CC_FREE
#    endif
#endif

#ifndef CC_MEMCPY
#    include <string.h>
#    define CC_MEMCPY(dest, src, count) memcpy(dest, src, count)
#endif

#ifndef CC_MEMSET
#    include <string.h>
#    define CC_MEMSET(ptr, value, num) memset(ptr, value, num)
#endif

#ifndef _OPENMP
#   ifndef CC_ATOMIC
#       define CC_ATOMIC
#   endif
#   ifndef CC_PARALLEL_FOR
#       define CC_PARALLEL_FOR
#   endif
#   ifndef CC_BARRIER
#       define CC_BARRIER
#   endif
#else
#   if defined(_WIN32)
#       ifndef CC_ATOMIC
#           define CC_ATOMIC          __pragma("omp atomic" )
#       endif
#       ifndef CC_PARALLEL_FOR
#           define CC_PARALLEL_FOR    __pragma("omp parallel for")
#       endif
#       ifndef CC_BARRIER
#           define CC_BARRIER         __pragma("omp barrier")
#       endif
#   else
#       ifndef CC_ATOMIC
#           define CC_ATOMIC          _Pragma("omp atomic" )
#       endif
#       ifndef CC_PARALLEL_FOR
#           define CC_PARALLEL_FOR    _Pragma("omp parallel for")
#       endif
#       ifndef CC_BARRIER
#           define CC_BARRIER         _Pragma("omp barrier")
#       endif
#   endif
#endif

/*******************************************************************************
 * NextPowerOfTwo -- Returns the upper power of two value
 *
 * if the input is already a power of two, its value is returned.
 *
 */
uint32_t cbf_NextPowerOfTwo(uint32_t x)
{
    x--;
    x|= x >> 1;
    x|= x >> 2;
    x|= x >> 4;
    x|= x >> 8;
    x|= x >> 16;
    x++;

    return x;
}


/*******************************************************************************
 * FindMSB -- Returns the position of the most significant bit
 *
 */
static inline int32_t cbf_FindMSB(uint32_t x)
{
    int32_t msb = 0;

    while (x > 1u) {
        ++msb;
        x = x >> 1;
    }

    return msb;
}

uint32_t cbf_BitCapacity(const uint32_t *cbf)
{
    return cbf[0];
}

uint32_t cbf_BitCount(const uint32_t *cbf)
{
    return cbf[1];
}

void cbf_SetBit(uint32_t *cbf, uint32_t bitID, uint32_t bitValue)
{
    const uint32_t capacity = cbf_BitCapacity(cbf);

    cbf[capacity + bitID] = bitValue;
}

uint32_t cbf_GetBit(const uint32_t *cbf, uint32_t bitID)
{
    const uint32_t capacity = cbf_BitCapacity(cbf);

    return cbf[capacity + bitID];
}

void cbf_Reduce(uint32_t *cbf)
{
    const uint32_t capacity = cbf_BitCapacity(cbf);
    int32_t depth = cbf_FindMSB(capacity);

    while (--depth >= 0) {
        uint32_t minHeapID = 1u << depth;
        uint32_t maxHeapID = 2u << depth;

        for (uint32_t heapID = minHeapID; heapID < maxHeapID; ++heapID) {
            cbf[heapID] = cbf[2u * heapID + 0] + cbf[2u * heapID + 1];
        }
    }
}

uint32_t *cbf_Create(uint32_t size)
{
    const uint32_t bitCount = cbf_NextPowerOfTwo(size);
    uint32_t *cbf = (uint32_t *)CC_MALLOC(sizeof(*cbf) * 2 * bitCount);

    cbf[0] = bitCount;

    for (uint32_t bitID = 0; bitID < bitCount; ++bitID) {
        cbf[bitCount + bitID] = 0u;
    }

    cbf_Reduce(cbf);

    return cbf;
}

void cbf_Release(uint32_t *cbf)
{
    CC_FREE(cbf);
}

uint32_t cbf_EncodeBit(const uint32_t *cbf, uint32_t bitID)
{
    uint32_t bitFieldSize = cbf_BitCapacity(cbf);
    uint32_t heapID = bitID + bitFieldSize;
    uint32_t handle = 0;

    while (heapID > 1u) {
        uint32_t siblingID = (uint)heapID & (~1u);
        uint32_t bitCount = cbf[siblingID];

        handle+= (heapID & 1u) * bitCount;
        heapID/= 2u;
    }

    return handle;
}

uint32_t cbf_DecodeBit(const uint32_t *cbf, uint32_t handle)
{
    uint32_t bitID = 1;
    uint32_t bitFieldSize = cbf_BitCapacity(cbf);

    while (bitID < bitFieldSize) {
        uint32_t heapValue = cbf[bitID * 2u];
        uint32_t b = handle < heapValue ? 0 : 1;

        bitID = 2u * bitID | b;
        handle -= heapValue * b;
    }

    return (bitID ^ bitFieldSize);
}

/*******************************************************************************
 * ComputeTwins -- Computes the twin of each half edge
 *
 * This routine is what effectively converts a traditional "indexed mesh"
 * into a halfedge mesh (in the case where all the primitives are the same).
 *
 */
typedef struct {
    int halfedgeID, hashID;
} TwinComputationData;

static int TwinComputationCompareCallback(const void *a, const void *b)
{
    const TwinComputationData *d1 = (const TwinComputationData *)a;
    const TwinComputationData *d2 = (const TwinComputationData *)b;

    if (d1->hashID > d2->hashID) {
        return 1;
    } else if (d1->hashID < d2->hashID) {
        return -1;
    } else {
        return 0;
    }
}

static int
BinarySearch(
    const TwinComputationData *data,
    int32_t hashID,
    int32_t beginID,
    int32_t endID
) {
    int32_t midID;

    if (beginID > endID)
       return -1; // not found

    midID = (beginID + endID) / 2;

    if (data[midID].hashID == hashID) {
        return data[midID].halfedgeID;
    } else if (hashID > data[midID].hashID) {
        return BinarySearch(data, hashID, midID + 1, endID);
    } else {
        return BinarySearch(data, hashID, beginID, midID - 1);
    }
}

static void ComputeTwins(cc_Mesh *mesh)
{
    const int32_t halfedgeCount = ccm_HalfedgeCount(mesh);
    const int32_t vertexCount = ccm_VertexCount(mesh);
    TwinComputationData *table = (TwinComputationData *)CC_MALLOC(halfedgeCount * sizeof(*table));

    CC_PARALLEL_FOR
    for (int32_t halfedgeID = 0; halfedgeID < halfedgeCount; ++halfedgeID) {
        const int32_t nextID = ccm_HalfedgeNextID(mesh, halfedgeID);
        const int32_t v0 = ccm_HalfedgeVertexID(mesh, halfedgeID);
        const int32_t v1 = ccm_HalfedgeVertexID(mesh, nextID);

        table[halfedgeID].halfedgeID = halfedgeID;
        table[halfedgeID].hashID = v0 + vertexCount * v1;
    }
    CC_BARRIER

    qsort(table, halfedgeCount, sizeof(table[0]), &TwinComputationCompareCallback);

    CC_PARALLEL_FOR
    for (int32_t halfedgeID = 0; halfedgeID < halfedgeCount; ++halfedgeID) {
        const int32_t nextID = ccm_HalfedgeNextID(mesh, halfedgeID);
        const int32_t v0 = ccm_HalfedgeVertexID(mesh, halfedgeID);
        const int32_t v1 = ccm_HalfedgeVertexID(mesh, nextID);
        const int32_t hashID = v1 + vertexCount * v0;
        const int32_t twinID = BinarySearch(table, hashID, 0, halfedgeCount - 1);

        mesh->halfedges[halfedgeID].twinID = twinID;
    }
    CC_BARRIER

    CC_FREE(table);
}


/*******************************************************************************
 * ComputeCreaseNeighbors -- Computes the neighbors of each crease
 *
 */
static void ComputeCreaseNeighbors(cc_Mesh *mesh)
{
    const int32_t edgeCount = ccm_EdgeCount(mesh);

CC_PARALLEL_FOR
    for (int32_t edgeID = 0; edgeID < edgeCount; ++edgeID) {
        const float sharpness = ccm_CreaseSharpness(mesh, edgeID);

        if (sharpness > 0.0f) {
            const int32_t halfedgeID = ccm_EdgeToHalfedgeID(mesh, edgeID);
            const int32_t nextID = ccm_HalfedgeNextID(mesh, halfedgeID);
            int32_t prevCreaseCount = 0;
            int32_t prevCreaseID = -1;
            int32_t nextCreaseCount = 0;
            int32_t nextCreaseID = -1;
            int32_t halfedgeIt;

            for (halfedgeIt = ccm_NextVertexHalfedgeID(mesh, halfedgeID);
                 halfedgeIt != halfedgeID && halfedgeIt >= 0;
                 halfedgeIt = ccm_NextVertexHalfedgeID(mesh, halfedgeIt)) {
                const float s = ccm_HalfedgeSharpness(mesh, halfedgeIt);

                if (s > 0.0f) {
                    prevCreaseID = ccm_HalfedgeEdgeID(mesh, halfedgeIt);
                    ++prevCreaseCount;
                }
            }

            if (prevCreaseCount == 1 && halfedgeIt == halfedgeID) {
                mesh->creases[edgeID].prevID = prevCreaseID;
            }

            if (ccm_HalfedgeSharpness(mesh, nextID) > 0.0f) {
                nextCreaseID = ccm_HalfedgeEdgeID(mesh, nextID);
                ++nextCreaseCount;
            }

            for (halfedgeIt = ccm_NextVertexHalfedgeID(mesh, nextID);
                 halfedgeIt != nextID && halfedgeIt >= 0;
                 halfedgeIt = ccm_NextVertexHalfedgeID(mesh, halfedgeIt)) {
                const float s = ccm_HalfedgeSharpness(mesh, halfedgeIt);
                const int32_t twinID = ccm_HalfedgeTwinID(mesh, halfedgeIt);

                // twin check is to avoid counting for halfedgeID
                if (s > 0.0f && twinID != halfedgeID) {
                    nextCreaseID = ccm_HalfedgeEdgeID(mesh, halfedgeIt);
                    ++nextCreaseCount;
                }
            }

            if (nextCreaseCount == 1 && halfedgeIt == nextID) {
                mesh->creases[edgeID].nextID = nextCreaseID;
            }
        }
    }
CC_BARRIER
}


/*******************************************************************************
 * MakeBoundariesSharp -- Tags boundary edges as sharp
 *
 * Following the Pixar standard, we tag boundary halfedges as sharp.
 * See "Subdivision Surfaces in Character Animation" by DeRose et al.
 * Note that we tag the sharpness value to 16 as subdivision can't go deeper
 * without overflowing 32-bit integers.
 *
 */
static void MakeBoundariesSharp(cc_Mesh *mesh)
{
    const int32_t edgeCount = ccm_EdgeCount(mesh);

CC_PARALLEL_FOR
    for (int32_t edgeID = 0; edgeID < edgeCount; ++edgeID) {
        const int32_t halfedgeID = ccm_EdgeToHalfedgeID(mesh, edgeID);
        const int32_t twinID = ccm_HalfedgeTwinID(mesh, halfedgeID);

        if (twinID < 0) {
            mesh->creases[edgeID].sharpness = 16.0f;
        }
    }
CC_BARRIER
}


/*******************************************************************************
 * LoadFaceMappings -- Computes the mappings for the faces of the mesh
 *
 */
static int32_t FaceScroll(int32_t id, int32_t direction, int32_t maxValue)
{
    const int32_t n = maxValue - 1;
    const int32_t d = direction;
    const int32_t u = (d + 1) >> 1; // in [0, 1]
    const int32_t un = u * n; // precomputation

    return (id == un) ? (n - un) : (id + d);
}

static int32_t
ScrollFaceHalfedgeID(
    int32_t halfedgeID,
    int32_t halfedgeFaceBeginID,
    int32_t halfedgeFaceEndID,
    int32_t direction
) {
    const int32_t faceHalfedgeCount = halfedgeFaceEndID - halfedgeFaceBeginID;
    const int32_t localHalfedgeID = halfedgeID - halfedgeFaceBeginID;
    const int32_t nextHalfedgeID = FaceScroll(localHalfedgeID,
                                              direction,
                                              faceHalfedgeCount);

    return halfedgeFaceBeginID + nextHalfedgeID;
}

static void LoadFaceMappings(cc_Mesh *mesh, const uint32_t *faceIterator)
{
    const int32_t halfedgeCount = ccm_HalfedgeCount(mesh);
    const int32_t faceCount = cbf_BitCount(faceIterator) - 1;

    mesh->faceToHalfedgeIDs = (int32_t *)CC_MALLOC(sizeof(int32_t) * faceCount);
    mesh->faceCount = faceCount;

CC_PARALLEL_FOR
    for (int32_t halfedgeID = 0; halfedgeID  < halfedgeCount; ++halfedgeID) {
        const int32_t tmp = cbf_EncodeBit(faceIterator, halfedgeID);
        const int32_t faceID = tmp - (cbf_GetBit(faceIterator, halfedgeID) ^ 1);

        mesh->halfedges[halfedgeID].faceID = faceID;
    }
CC_BARRIER

CC_PARALLEL_FOR
    for (int32_t faceID = 0; faceID < faceCount; ++faceID) {
        mesh->faceToHalfedgeIDs[faceID] = cbf_DecodeBit(faceIterator, faceID);
    }
CC_BARRIER


CC_PARALLEL_FOR
    for (int32_t halfedgeID = 0; halfedgeID  < halfedgeCount; ++halfedgeID) {
        const int32_t faceID = mesh->halfedges[halfedgeID].faceID;
        const int32_t beginID = cbf_DecodeBit(faceIterator, faceID);
        const int32_t endID = cbf_DecodeBit(faceIterator, faceID + 1);
        const int32_t nextID = ScrollFaceHalfedgeID(halfedgeID, beginID, endID, +1);
        const int32_t prevID = ScrollFaceHalfedgeID(halfedgeID, beginID, endID, -1);

        mesh->halfedges[halfedgeID].nextID = nextID;
        mesh->halfedges[halfedgeID].prevID = prevID;
    }
CC_BARRIER
}


/*******************************************************************************
 * LoadEdgeMappings -- Computes the mappings for the edges of the mesh
 *
 * Catmull-Clark subdivision requires access to the edges of an input mesh.
 * Since we are dealing with a half-edge representation, we virtually
 * have to iterate the half-edges in a sparse way (an edge is a pair of
 * neighboring half-edges in the general case, except for boundary edges
 * where it only consists of a single half-edge).
 * This function builds a data-structure that allows to do just that:
 * for each halfedge pair, we only consider the one that has the largest
 * halfedgeID. This allows to treat boundary and regular edges seamlessly.
 *
 */
static void LoadEdgeMappings(cc_Mesh *mesh)
{
    const int32_t halfedgeCount = ccm_HalfedgeCount(mesh);
    cc_Halfedge *halfedges = mesh->halfedges;
    uint32_t *edgeIterator = cbf_Create(halfedgeCount);
    int32_t edgeCount;

CC_PARALLEL_FOR
    for (int32_t halfedgeID = 0; halfedgeID < halfedgeCount; ++halfedgeID) {
        const int32_t twinID = ccm_HalfedgeTwinID(mesh, halfedgeID);
        const int32_t bitValue = halfedgeID > twinID ? 1 : 0;

        cbf_SetBit(edgeIterator, halfedgeID, bitValue);
    }
CC_BARRIER

    cbf_Reduce(edgeIterator);
    edgeCount = cbf_BitCount(edgeIterator);

    mesh->edgeToHalfedgeIDs = (int32_t *)CC_MALLOC(sizeof(int32_t) * edgeCount);
    mesh->edgeCount = edgeCount;

CC_PARALLEL_FOR
    for (int32_t halfedgeID = 0; halfedgeID  < halfedgeCount; ++halfedgeID) {
        const int32_t twinID = ccm_HalfedgeTwinID(mesh, halfedgeID);
        const int32_t bitID = halfedgeID > twinID ? halfedgeID : twinID;

        halfedges[halfedgeID].edgeID = cbf_EncodeBit(edgeIterator, bitID);
    }
CC_BARRIER

CC_PARALLEL_FOR
    for (int32_t edgeID = 0; edgeID < edgeCount; ++edgeID) {
        mesh->edgeToHalfedgeIDs[edgeID] = cbf_DecodeBit(edgeIterator, edgeID);
    }
CC_BARRIER

    cbf_Release(edgeIterator);
}


/*******************************************************************************
 * LoadVertexHalfedges -- Computes an iterator over one half-edge per vertex
 *
 * Catmull-Clark subdivision requires access to the half-edges that surround
 * the vertices of an input mesh.
 * This function determines a half-edge ID that starts from a
 * given vertex within that vertex. We distinguish two cases:
 * 1- If the vertex is a lying on a boundary, we stored the halfedge that
 * allows for iteration in the forward sense.
 * 2- Otherwise we store the largest half-edge ID.
 *
 */
static void LoadVertexHalfedges(cc_Mesh *mesh)
{
    const int32_t halfedgeCount = ccm_HalfedgeCount(mesh);
    const int32_t vertexCount = ccm_VertexCount(mesh);

    mesh->vertexToHalfedgeIDs = (int32_t *)CC_MALLOC(sizeof(int32_t) * vertexCount);

CC_PARALLEL_FOR
    for (int32_t halfedgeID = 0; halfedgeID < halfedgeCount; ++halfedgeID) {
        const int32_t vertexID = ccm_HalfedgeVertexID(mesh, halfedgeID);
        int32_t maxHalfedgeID = halfedgeID;
        int32_t boundaryHalfedgeID = halfedgeID;
        int32_t iterator;

        for (iterator = ccm_NextVertexHalfedgeID(mesh, halfedgeID);
             iterator >= 0 && iterator != halfedgeID;
             iterator = ccm_NextVertexHalfedgeID(mesh, iterator)) {
            maxHalfedgeID = maxHalfedgeID > iterator ? maxHalfedgeID : iterator;
            boundaryHalfedgeID = iterator;
        }

        // affect max half-edge ID to vertex
        if /*boundary involved*/ (iterator < 0) {
            if (halfedgeID == boundaryHalfedgeID) {
                mesh->vertexToHalfedgeIDs[vertexID] = boundaryHalfedgeID;
            }
        } else {
            if (halfedgeID == maxHalfedgeID) {
                mesh->vertexToHalfedgeIDs[vertexID] = maxHalfedgeID;
            }
        }
    }
CC_BARRIER
}

void
LoadPtexData(const PtexMesh *ptex, cc_Mesh *mesh, uint32_t *faceIterator)
{
    const int32_t faceCount = ptex->faceVertCounts.size();
    cc_Halfedge *halfedges = mesh->halfedges;
    cc_VertexPoint *vertexPoints = mesh->vertexPoints;
    int32_t bitID = 0;

    // load halfedge data
    for (int32_t halfedgeID = 0; halfedgeID < mesh->halfedgeCount; ++halfedgeID) {
        halfedges[halfedgeID].vertexID = ptex->faceVertIndices[halfedgeID];
    }

    // load vertex data
    for (int32_t vertexID = 0; vertexID < mesh->vertexCount; ++vertexID) {
        vertexPoints[vertexID].x = ptex->vertPositions[3 * vertexID + 0];
        vertexPoints[vertexID].y = ptex->vertPositions[3 * vertexID + 1];
        vertexPoints[vertexID].z = ptex->vertPositions[3 * vertexID + 2];
    }

    // build face iterator
    for (int32_t faceID = 0; faceID < faceCount; ++faceID) {
        cbf_SetBit(faceIterator, bitID, 1u);
        bitID+= ptex->faceVertCounts[faceID];
    }
    cbf_SetBit(faceIterator, bitID, 1u);
    cbf_Reduce(faceIterator);
    LoadFaceMappings(mesh, faceIterator);
}

cc_Mesh *PtexToHalfedges(const PtexMesh *ptex)
{
    const int32_t vertexCount = ptex->vertPositions.size() / 3u;
    const int32_t halfedgeCount = ptex->faceVertIndices.size();
    uint32_t *faceIterator;
    cc_Mesh *mesh;

    mesh = (cc_Mesh *)CC_MALLOC(sizeof(*mesh));
    mesh->halfedgeCount = halfedgeCount;
    mesh->halfedges = (cc_Halfedge *)CC_MALLOC(sizeof(cc_Halfedge) * halfedgeCount);
    mesh->vertexCount = vertexCount;
    mesh->vertexPoints = (cc_VertexPoint *)CC_MALLOC(sizeof(cc_VertexPoint) * vertexCount);
    mesh->uvCount = 0u;
    mesh->uvs = NULL;
    faceIterator = cbf_Create(halfedgeCount + 1);

    LoadPtexData(ptex, mesh, faceIterator);
    ComputeTwins(mesh);
    LoadEdgeMappings(mesh);
    LoadVertexHalfedges(mesh);

    if (true) {
        const int32_t creaseCount = ccm_EdgeCount(mesh);

        mesh->creases = (cc_Crease *)CC_MALLOC(sizeof(cc_Crease) * creaseCount);

CC_PARALLEL_FOR
        for (int32_t creaseID = 0; creaseID < creaseCount; ++creaseID) {
            mesh->creases[creaseID].nextID = creaseID;
            mesh->creases[creaseID].prevID = creaseID;
            mesh->creases[creaseID].sharpness = 0.0f;
        }
CC_BARRIER

        MakeBoundariesSharp(mesh);
    }

    ComputeCreaseNeighbors(mesh);

    cbf_Release(faceIterator);

    return mesh;
}
