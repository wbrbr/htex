#define STB_IMAGE_IMPLEMENTATION
#include <vector>
#include "Htexture.h"
#include "Ptexture.h"
#include "glm/glm.hpp"
#include "PtexToHalfedges.h"

struct Color {
    uint8_t r,g,b,a;
};

int min(int a, int b)
{
    return a < b ? a : b;
}

int max(int a, int b)
{
    return a > b ? a : b;
}

glm::vec2 BarycentricInterpolation(float u, float v, glm::vec2 v0, glm::vec2 v1, glm::vec2 v2)
{
    return (1-u-v)*v0 + u*v1 + v*v2;
}

int LocalHalfedgeIndex(const cc_Mesh* mesh, int halfedgeID)
{
    return halfedgeID - ccm_FaceToHalfedgeID(mesh, ccm_HalfedgeFaceID(mesh, halfedgeID));
}

void generateTexture(const cc_Mesh* mesh, int quadID, PtexTexture* ptex, std::vector<uint8_t>& texels, uint8_t& logu, uint8_t& logv) {
    glm::vec2 corners[] = {
            glm::vec2(0,0),
            glm::vec2(1,0),
            glm::vec2(1,1),
            glm::vec2(0,1)
    };

    int halfedge1 = ccm_EdgeToHalfedgeID(mesh, quadID);
    int halfedge2 = ccm_HalfedgeTwinID(mesh, halfedge1);

    int halfedgeMax = halfedge1 > halfedge2 ? halfedge1 : halfedge2;
    int halfedgeMin = halfedge1 < halfedge2 ? halfedge1 : halfedge2;

    int faceMax = ccm_HalfedgeFaceID(mesh, halfedgeMax);
    int faceMin = ccm_HalfedgeFaceID(mesh, halfedgeMin);

    int indexHalfedgeMax = LocalHalfedgeIndex(mesh, halfedgeMax);
    int indexHalfedgeMin = LocalHalfedgeIndex(mesh, halfedgeMin);

    glm::vec2 max_uv0(0.5, 0.5);
    glm::vec2 max_uv1 = corners[indexHalfedgeMax];
    glm::vec2 max_uv2 = corners[(indexHalfedgeMax+1)%4];

    glm::vec2 min_uv0(0.5, 0.5);
    glm::vec2 min_uv1 = corners[indexHalfedgeMin];
    glm::vec2 min_uv2 = corners[(indexHalfedgeMin+1)%4];

    PtexFilter* filter = PtexFilter::getFilter(ptex, PtexFilter::Options(PtexFilter::f_point));

    Ptex::Res res1 = ptex->getFaceInfo(faceMax).res;
    Ptex::Res res2 = ptex->getFaceInfo(faceMin).res;
    uint8_t maxres = max(max(max(res1.ulog2, res1.vlog2), res2.ulog2), res2.vlog2);

    logu = maxres - 1;
    logv = maxres - 1;

    int width = 1 << logu;
    int height = 1 << logv;

    texels.resize(4*width*height);

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            float u = ((float)j + .5f) / (float)width;
            float v = ((float)i + .5f) / (float)height;

            int px = i*width+j;

            glm::vec2 texcoords(0);
            int faceid = -1;
            if (u+v <= 1.f) {
                texcoords = BarycentricInterpolation(u, v, max_uv0, max_uv1, max_uv2);
                faceid = faceMax;
            } else {
                if (halfedgeMin >= 0) {
                    texcoords = BarycentricInterpolation(1-u, 1-v, min_uv0, min_uv1, min_uv2);
                    faceid = faceMin;
                } else {
                    texcoords = BarycentricInterpolation(1-u, 1-v, max_uv0, max_uv2, max_uv1);
                    faceid = faceMax;
                }
            }

            float data[4];
            // TODO: compute filter footprint
            filter->eval(data, 0, 4, faceid, texcoords.x, texcoords.y, 0, 0, 0, 0);
            texels[4*px+0] = (uint8_t)(data[0]*255.f);
            texels[4*px+1] = (uint8_t)(data[1]*255.f);
            texels[4*px+2] = (uint8_t)(data[2]*255.f);
            texels[4*px+3] = (uint8_t)(data[3]*255.f);
        }
    }
}

int main(int argc, char** argv) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input ptex> <output path>\n", argv[0]);
        return 1;
    }

    const char* inputPath = argv[1];
    const char* outputPath = argv[2];

    Ptex::String err;
    PtexTexture* ptex = PtexTexture::open(inputPath, err);
    if (!ptex) {
        fprintf(stderr, "Failed to open Ptex file: %s\n", err.c_str());
        return 1;
    }


    PtexMesh* ptex_mesh = LoadPtexMesh(ptex);
    cc_Mesh* halfedge_mesh = PtexToHalfedges(ptex_mesh);

    HtexWriter* writer = HtexWriter::open(outputPath, halfedge_mesh, Htex::mt_quad, Htex::dt_uint8, 4, 3, err);
    if (!writer) {
        fprintf(stderr, "Failed to create HtexWriter: %s\n", err.c_str());
        exit(1);
    }

    std::vector<uint8_t> texels;
    for (int quadID = 0; quadID < halfedge_mesh->edgeCount; quadID++) {
        uint8_t logu, logv;
        generateTexture(halfedge_mesh, quadID, ptex, texels, logu, logv);
        Htex::QuadInfo quadInfo{Htex::Res(logu, logv), quadID};
        if (!writer->writeQuad(quadID, quadInfo, texels.data())) {
            writer->close(err);
            fprintf(stderr, "Failed to write quad %d: %s\n", quadID, err.c_str());
            exit(1);
        }
        if (quadID % 100 == 0) {
            printf("%i / %i\n", quadID, halfedge_mesh->edgeCount);
        }
    }
    printf("Done resampling\n");

    if (!writer->close(err)) {
        fprintf(stderr, "Failed to write Htex file: %s\n", err.c_str());
        exit(1);
    }
}
