#include "embree3/rtcore.h"
#include "Htexture.h"
#include "glm/glm.hpp"
#include <vector>

glm::vec3 FaceBarycenter(cc_Mesh* mesh, int faceID)
{
    int startHalfedge = ccm_FaceToHalfedgeID(mesh, faceID);
    int currentHalfedge = startHalfedge;

    glm::vec3 sum(0);
    int n = 0;
    do {
        cc_VertexPoint vp = ccm_HalfedgeVertexPoint(mesh, currentHalfedge);
        glm::vec3 v(vp.x, vp.y, vp.z);
        sum += v;
        n++;

        currentHalfedge = ccm_HalfedgeNextID(mesh, currentHalfedge);
    } while (currentHalfedge != startHalfedge);

    sum /= n;
    return sum;
}

glm::vec3 ComputeVertexNormal(cc_Mesh* mesh, int vertexID, const std::vector<glm::vec3>& halfedgeNormals)
{
    int startHalfedge = ccm_VertexPointToHalfedgeID(mesh, vertexID);
    glm::vec3 n = glm::vec3(0);
    int currentHEdge = startHalfedge;
    do {
        n += halfedgeNormals[currentHEdge];
        currentHEdge = ccm_HalfedgeTwinID(mesh, currentHEdge);
        if (currentHEdge < 0) break;
        currentHEdge = ccm_HalfedgeNextID(mesh, currentHEdge);
    } while (currentHEdge != startHalfedge);

    // boundary, we do a backwards traversal
    if (currentHEdge != startHalfedge) {
        currentHEdge = ccm_HalfedgeTwinID(mesh, ccm_HalfedgePrevID(mesh, startHalfedge));
        while (currentHEdge >= 0) {
            n += halfedgeNormals[currentHEdge];
            currentHEdge = ccm_HalfedgeTwinID(mesh, ccm_HalfedgePrevID(mesh, currentHEdge));
        }
    }

    return glm::normalize(n);
}

glm::vec3 ComputeFacePointNormal(cc_Mesh* mesh, int faceID, const std::vector<glm::vec3>& halfedgeNormals)
{
    int startHalfedge = ccm_FaceToHalfedgeID(mesh, faceID);
    glm::vec3 n = glm::vec3(0);
    int currentHEdge = startHalfedge;
    do {
        n += ComputeVertexNormal(mesh, ccm_HalfedgeVertexID(mesh, currentHEdge), halfedgeNormals);
        currentHEdge = ccm_HalfedgeNextID(mesh, currentHEdge);
    } while(currentHEdge != startHalfedge);
    return glm::normalize(n);
}

glm::vec3 RandomUnitVectorHemisphere(glm::vec3 n)
{
    glm::vec3 res;
    do {
        res.x = (float)drand48()*2-1;
        res.y = (float)drand48()*2-1;
        res.z = (float)drand48()*2-1;
    } while (glm::dot(res, res) > 1 || glm::dot(res, n) < 0);
    return glm::normalize(res);
}

void GenerateTexture(cc_Mesh* mesh, int halfedgeID, uint8_t* quad_texels, int width, int height, bool is_max_halfedge, RTCScene scene, RTCGeometry geom)
{
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width-i; j++) {
            float u = ((float)j + .5f) / (float)width;
            float v = ((float)i + .5f) / (float)height;

            int px = i*width+j;
            if (!is_max_halfedge) {
                px = (height-1-i)*width + width-1-j;
            }

            float P[3];
            float N[3];

            rtcInterpolate0(geom, halfedgeID, u, v, RTC_BUFFER_TYPE_VERTEX, 0, P, 3);
            rtcInterpolate0(geom, halfedgeID, u, v, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, N, 3);

            glm::vec3 n = glm::normalize(glm::vec3(N[0], N[1], N[2]));

            const int sample_count = 16;
            int visible = 0;

            RTCIntersectContext context;
            rtcInitIntersectContext(&context);
            context.flags = RTC_INTERSECT_CONTEXT_FLAG_COHERENT;


            glm::vec3 dir_sum(0);
//#pragma omp parallel for default(none) shared(sample_count, dir_sum, P, n, scene, visible)
            for (int k = 0; k < sample_count / 16; k++) {
                int valid[16];
                for (int l = 0; l < 16; l++) {
                    valid[l] = 16*k+l < sample_count ? -1 : 0;
                }

                RTCRay16 rays;
                for (int l = 0; l < 16; l++) {
                    glm::vec3 dir = RandomUnitVectorHemisphere(n);

                    rays.org_x[l] = P[0];
                    rays.org_y[l] = P[1];
                    rays.org_z[l] = P[2];
                    rays.dir_x[l] = dir.x;
                    rays.dir_y[l] = dir.y;
                    rays.dir_z[l] = dir.z;
                    rays.tnear[l] = 1e-5;
                    rays.tfar[l] = INFINITY;
                }

                rtcOccluded16(valid, scene, &context, &rays);

                for (int l = 0; l < 16; l++) {
                    if (rays.tfar[l] != -INFINITY) {
                        visible++;
                    }
                }
            }



            /*if (rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID) {
                quad_texels[4*px+0] = 0x00;
                quad_texels[4*px+1] = 0x00;
                quad_texels[4*px+2] = 0xff;
                quad_texels[4*px+3] = 0xff;
            } else {
                quad_texels[4*px+0] = 0xff;
                quad_texels[4*px+1] = 0x00;
                quad_texels[4*px+2] = 0x00;
                quad_texels[4*px+3] = 0xff;
            }*/

            dir_sum /= sample_count;
            float ao = (float)visible / (float)sample_count;
#if 1
            quad_texels[4*px+0] = (uint8_t)(ao*255.f);
            quad_texels[4*px+1] = (uint8_t)(ao*255.f);
            quad_texels[4*px+2] = (uint8_t)(ao*255.f);
            quad_texels[4*px+3] = 0xff;
#elif 0
            quad_texels[4*px+0] = (uint8_t)(fabs(dir_sum.x)*255.f);
            quad_texels[4*px+1] = (uint8_t)(fabs(dir_sum.y)*255.f);
            quad_texels[4*px+2] = (uint8_t)(fabs(dir_sum.z)*255.f);
            quad_texels[4*px+3] = 0xff;
#elif 0
            quad_texels[4*px+0] = 0xff;
            quad_texels[4*px+1] = 0x00;
            quad_texels[4*px+2] = 0x00;
            quad_texels[4*px+3] = 0xff;
#elif 0
            quad_texels[4*px+0] = (uint8_t)glm::clamp(p.x * 255, 0.f, 255.f);
            quad_texels[4*px+1] = (uint8_t)glm::clamp(p.y * 255, 0.f, 255.f);
            quad_texels[4*px+2] = (uint8_t)glm::clamp(p.z * 255, 0.f, 255.f);
            quad_texels[4*px+3] = 0xff;
#elif 0
            quad_texels[4*px+0] = (uint8_t)(u*255.f);
            quad_texels[4*px+1] = (uint8_t)(v*255.f);
            quad_texels[4*px+2] = 0;
            quad_texels[4*px+3] = 0xff;
#elif 0
            quad_texels[4*px+0] = (uint8_t)glm::clamp(N[0]*255.f, 0.f, 255.f);
            quad_texels[4*px+1] = (uint8_t)glm::clamp(N[1]*255.f, 0.f, 255.f);
            quad_texels[4*px+2] = (uint8_t)glm::clamp(N[2]*255.f, 0.f, 255.f);
            quad_texels[4*px+3] = 0xff;
#endif
        }
    }
}

void embree_error(void* userPtr, RTCError code, const char* str)
{
    fprintf(stderr,  "Embree error: %s\n", str);
    exit(1);
}

int main(int argc, char** argv) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <ccm file> <log2 resolution> <output path>\nYou can use the obj_to_ccm program to generate a .ccm from a .obj file\n", argv[0]);
        return 1;
    }

    const char* ccmPath = argv[1];
    int log2_res = atoi(argv[2]);
    if (log2_res <= 0 || log2_res > 255) {
        fprintf(stderr, "Invalid resolution\n");
        return 1;
    }
    const char* outputPath = argv[3];

    cc_Mesh* halfedge_mesh = ccm_Load(ccmPath);

    Htex::String err;
    HtexWriter* writer = HtexWriter::open(outputPath, halfedge_mesh, Htex::mt_quad, Htex::dt_uint8, 4, 3, err);
    if (!writer) {
        fprintf(stderr, "Failed to create HtexWriter: %s\n", err.c_str());
        exit(1);
    }

    const int width = 1 << log2_res;
    const int height = 1 << log2_res;

    printf("Using a %ix%i texture for each quad\n", width, height);

    RTCDevice device = rtcNewDevice(NULL);
    rtcSetDeviceErrorFunction(device, embree_error, NULL);
    RTCScene scene = rtcNewScene(device);
    RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);
    float* vertexBuffer = (float*) rtcSetNewGeometryBuffer(geom,
                                                           RTC_BUFFER_TYPE_VERTEX,
                                                           0,
                                                           RTC_FORMAT_FLOAT3,
                                                           3*sizeof(float),
                                                           ccm_VertexCount(halfedge_mesh) + ccm_FaceCount(halfedge_mesh));

    for (int vid = 0; vid < ccm_VertexCount(halfedge_mesh); vid++) {
        cc_VertexPoint vp = ccm_VertexPoint(halfedge_mesh, vid);
        vertexBuffer[3*vid+0] = vp.x;
        vertexBuffer[3*vid+1] = vp.y;
        vertexBuffer[3*vid+2] = vp.z;
    }
    for (int faceid = 0; faceid < ccm_FaceCount(halfedge_mesh); faceid++) {
        glm::vec3 barycenter = FaceBarycenter(halfedge_mesh, faceid);
        int idx = halfedge_mesh->vertexCount+faceid;
        vertexBuffer[3*idx+0] = barycenter.x;
        vertexBuffer[3*idx+1] = barycenter.y;
        vertexBuffer[3*idx+2] = barycenter.z;
    }


    unsigned int* indexBuffer = (unsigned int*) rtcSetNewGeometryBuffer(geom,
                                                                        RTC_BUFFER_TYPE_INDEX,
                                                                        0,
                                                                        RTC_FORMAT_UINT3,
                                                                        3*sizeof(unsigned int),
                                                                        ccm_HalfedgeCount(halfedge_mesh));

    for (int halfedgeID = 0; halfedgeID < ccm_HalfedgeCount(halfedge_mesh); halfedgeID++) {
        indexBuffer[3*halfedgeID+0] = (unsigned int)(halfedge_mesh->vertexCount + ccm_HalfedgeFaceID(halfedge_mesh, halfedgeID));
        indexBuffer[3*halfedgeID+1] = (unsigned int)ccm_HalfedgeVertexID(halfedge_mesh, halfedgeID);
        indexBuffer[3*halfedgeID+2] = (unsigned int)ccm_HalfedgeVertexID(halfedge_mesh, ccm_HalfedgeNextID(halfedge_mesh, halfedgeID));
    }


    std::vector<glm::vec3> halfedgeNormals(ccm_HalfedgeCount(halfedge_mesh));
    for (int halfedgeID = 0; halfedgeID < ccm_HalfedgeCount(halfedge_mesh); halfedgeID++) {
        glm::vec3 vertices[3];
        for (int k = 0; k < 3; k++) {
            unsigned int vid = indexBuffer[3*halfedgeID+k];
            vertices[k] = glm::vec3(vertexBuffer[3*vid+0], vertexBuffer[3*vid+1], vertexBuffer[3*vid+2]);
        }

        glm::vec3 n = glm::normalize(glm::cross(vertices[1]-vertices[0], vertices[2]-vertices[0]));
        halfedgeNormals[halfedgeID] = n;
    }

    rtcSetGeometryVertexAttributeCount(geom, 1);
    float* normalBuffer = (float*) rtcSetNewGeometryBuffer(geom,
                                                           RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE,
                                                           0,
                                                           RTC_FORMAT_FLOAT3,
                                                           3*sizeof(float),
                                                           ccm_VertexCount(halfedge_mesh) + ccm_FaceCount(halfedge_mesh));
    for (int vid = 0; vid < ccm_VertexCount(halfedge_mesh); vid++) {
        glm::vec3 n = ComputeVertexNormal(halfedge_mesh, vid, halfedgeNormals);
        normalBuffer[3*vid+0] = n.x;
        normalBuffer[3*vid+1] = n.y;
        normalBuffer[3*vid+2] = n.z;
    }
    for (int faceid = 0; faceid < ccm_FaceCount(halfedge_mesh); faceid++) {
        glm::vec3 n = ComputeFacePointNormal(halfedge_mesh, faceid, halfedgeNormals);
        int idx = halfedge_mesh->vertexCount + faceid;
        normalBuffer[3*idx+0] = n.x;
        normalBuffer[3*idx+1] = n.y;
        normalBuffer[3*idx+2] = n.z;
    }

    rtcCommitGeometry(geom);
    rtcAttachGeometry(scene, geom);
    rtcReleaseGeometry(geom);
    rtcCommitScene(scene);

    uint8_t* texels = new uint8_t[width*height*4];
    for (int i = 0; i < width*height; i++) {
        texels[4*i+0] = 0x00;
        texels[4*i+1] = 0xff;
        texels[4*i+2] = 0x00;
        texels[4*i+3] = 0xff;
    }
    for (int quadID = 0; quadID < halfedge_mesh->edgeCount; quadID++) {
        int halfedge1 = ccm_EdgeToHalfedgeID(halfedge_mesh, quadID);
        int halfedge2 = ccm_HalfedgeTwinID(halfedge_mesh, halfedge1);
        int halfedge_max = (halfedge1 > halfedge2) ? halfedge1 : halfedge2;
        int halfedge_min = (halfedge1 < halfedge2) ? halfedge1 : halfedge2;

        GenerateTexture(halfedge_mesh, halfedge_max, texels, width, height, true, scene, geom);
        if (halfedge_min >= 0) {
            GenerateTexture(halfedge_mesh, halfedge_min, texels, width, height, false, scene, geom);
        } else {
            // TODO
        }

        Htex::QuadInfo quadInfo{Htex::Res(log2_res, log2_res), quadID};
        if (!writer->writeQuad(quadID, quadInfo, texels)) {
            writer->close(err);
            fprintf(stderr, "Failed to write quad %d: %s\n", quadID, err.c_str());
            exit(1);
        }

        if ((quadID + 1) % 100 == 0) {
            printf("%i / %i\n", quadID+1, halfedge_mesh->edgeCount);
        }
    }

    delete[] texels;


    if (!writer->close(err)) {
        fprintf(stderr, "Failed to write Htex file: %s\n", err.c_str());
        exit(1);
    }
}
