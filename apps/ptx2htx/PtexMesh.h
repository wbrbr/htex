#ifndef PTEX_MESH
#define PTEX_MESH

// ptex mesh data
typedef struct {
    std::vector<float> vertPositions;
    std::vector<int32_t> faceVertIndices;
    std::vector<int32_t> faceVertCounts;
    std::vector<int32_t> adjFaces;
    std::vector<int32_t> adjEdges;
} PtexMesh;

// load ptex
PtexMesh *LoadPtexMesh(const char *pathToPtexFile, Ptex::String &error);
PtexMesh *LoadPtexMesh(PtexTexture *ptex);

#endif
