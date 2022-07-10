#include <vector>
#include <stdint.h>

#include "Ptexture.h"
#include "PtexMesh.h"

template <typename T>
bool
LoadMetaData(
    const char *name,
    PtexMetaData *metaData,
    std::vector<T> &values
) {
    bool metaDataFound = false;

    for (int32_t keyID = 0; keyID < metaData->numKeys(); keyID++) {
        const char* key;
        Ptex::MetaDataType type;
        metaData->getKey(keyID, key, type);

        if (!strcmp(key, name)) {
            int32_t valueCount = 0;
            const T *value;

            metaData->getValue(key, value, valueCount);
            values.resize(valueCount);

            for (int32_t valueID = 0; valueID < valueCount; ++valueID) {
                values[valueID] = value[valueID];
            }

            metaDataFound = true;
        }
    }

    return metaDataFound;
}

PtexMesh *LoadPtexMesh(PtexTexture *ptex)
{
    PtexMesh *mesh = new PtexMesh();
    PtexMetaData *metaData = ptex->getMetaData();

    if (!LoadMetaData("PtexVertPositions", metaData, mesh->vertPositions)) {
        delete mesh;

        return NULL;
    }

    if (!LoadMetaData("PtexFaceVertCounts", metaData, mesh->faceVertCounts)) {
        delete mesh;

        return NULL;
    }

    if (!LoadMetaData("PtexFaceVertIndices", metaData, mesh->faceVertIndices)) {
        delete mesh;

        return NULL;
    }

    mesh->adjFaces.resize(4*ptex->numFaces());
    mesh->adjEdges.resize(4*ptex->numFaces());

    for (int faceid = 0; faceid < ptex->numFaces(); faceid++) {
        Ptex::FaceInfo faceInfo = ptex->getFaceInfo(faceid);

        for (int eid = 0; eid < 4; eid++) {
            mesh->adjFaces[4*faceid+eid] = faceInfo.adjfaces[eid];
            mesh->adjEdges[4*faceid+eid] = faceInfo.adjedge(eid);
        }
    }

    return mesh;
}

PtexMesh *LoadPtexMesh(const char *pathToPtexFile, Ptex::String &error)
{
    PtexPtr<PtexTexture> ptex(PtexTexture::open(pathToPtexFile, error));

    if (ptex) {
        return LoadPtexMesh(ptex);
    }

    return NULL;
}
