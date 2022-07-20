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

void OrthonormalBasis(glm::vec3 n, glm::vec3& b1, glm::vec3& b2)
{
    if (fabs(n.y) > 0.9) {
        b1 = glm::normalize(glm::cross(n, glm::vec3(1,0,0)));
    } else {
        b1 = glm::normalize(glm::cross(n, glm::vec3(0,1,0)));
    }

    b2 = glm::normalize(glm::cross(n, b1));

    float eps = 1e-5;
    assert(fabs(glm::dot(b1, n)) < eps);
    assert(fabs(glm::dot(b2, n)) < eps);
    assert(fabs(glm::dot(b1, b2)) < eps);
    assert(fabs(glm::length(n) - 1) < eps);
    assert(fabs(glm::length(b1) -1) < eps);
    assert(fabs(glm::length(b2) - 1) < eps);
}

glm::vec3 RandomCosineWeightedHemisphere(glm::vec3 n)
{
    glm::vec3 b1, b2;
    OrthonormalBasis(n, b1, b2);

    float theta = 2.f*(float)M_PI*(float)drand48();
    float r = sqrtf((float)drand48());

    float x = r*cosf(theta);
    float y = r*sinf(theta);
    float z = sqrtf(1-r*r);

    return x*b1 + y*b2 + z*n;
}

glm::vec3 RandomUnitVectorHemisphere(glm::vec3 n)
{
    glm::vec3 res;
    do {
        res.x = (float)drand48()*2-1;
        res.y = (float)drand48()*2-1;
        res.z = (float)drand48()*2-1;
    //} while (glm::dot(res, res) > 1 || abs(glm::dot(glm::normalize(res), n)) < 5e-2);
    } while (glm::dot(res, res) > 1);

    if (glm::dot(res, n) < 0) {
        res *= -1;
    }

    return glm::normalize(res);
}

float ComputeAO(glm::vec3 p, glm::vec3 n, int sample_count, RTCScene scene)
{
    double sum = 0;

    RTCIntersectContext context;
    rtcInitIntersectContext(&context);
    context.flags = RTC_INTERSECT_CONTEXT_FLAG_COHERENT;


    std::vector<glm::vec3> directions(sample_count);
    for (int s = 0; s < sample_count; s++) {
        //directions[s] = RandomUnitVectorHemisphere(n);
        directions[s] = RandomCosineWeightedHemisphere(n);
    }

    for (int s = 0; s < sample_count; s++) {
        RTCRayHit rayhit;
        glm::vec3 dir = directions[s];

        rayhit.ray.org_x = p.x;
        rayhit.ray.org_y = p.y;
        rayhit.ray.org_z = p.z;
        rayhit.ray.dir_x = dir.x;
        rayhit.ray.dir_y = dir.y;
        rayhit.ray.dir_z = dir.z;
        rayhit.ray.tnear = 1e-2;

        rayhit.ray.tfar = INFINITY;
        rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;

        rtcOccluded1(scene, &context, &rayhit.ray);
        //rtcIntersect1(scene, &context, &rayhit);


        if (rayhit.ray.tfar != -INFINITY) {
            //if (rayhit.hit.geomID == RTC_INVALID_GEOMETRY_ID) {
            sum += 1;
            /*double pdf = 0.5f / M_PI;
            double val = glm::dot(dir, n) / (float)M_PI;
            sum += val / pdf;*/
        } else {
#if 0
            printf("%f %f %f\n", dir.x,dir.y, dir.z);
                    fprintf(stderr, "intersection: %f %f %f / %f %f %f\n", rayhit.ray.dir_x, rayhit.ray.dir_y, rayhit.ray.dir_z,
                            n.x, n.y, n.z);
                    //fprintf(stderr, "t: %f\n", rayhit.ray.tfar);
                    fprintf(stderr, "%f\n\n", glm::dot(dir, n));
#endif
        }
    }

    float ao = (float)(glm::clamp(sum / (double)sample_count, 0., 1.));

    return ao;
}

void GenerateTexture(cc_Mesh* mesh, int halfedgeID, uint8_t* quad_texels, int width, int height, bool is_max_halfedge, RTCScene scene, RTCGeometry geom, glm::vec3* halfedgeNormals)
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
#if 1
            rtcInterpolate0(geom, halfedgeID, u, v, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, N, 3);
            glm::vec3 n = glm::normalize(glm::vec3(N[0], N[1], N[2]));
#else
            glm::vec3 n = halfedgeNormals[halfedgeID];
#endif

            glm::vec3 p(P[0], P[1], P[2]);
            float ao = ComputeAO(p, n, 160, scene);

            quad_texels[4*px+0] = (uint8_t)(ao*255.f);
            quad_texels[4*px+1] = (uint8_t)(ao*255.f);
            quad_texels[4*px+2] = (uint8_t)(ao*255.f);
            quad_texels[4*px+3] = 0xff;
        }
    }
}

void embree_error(void* userPtr, RTCError code, const char* str)
{
    fprintf(stderr,  "Embree error: %s\n", str);
    exit(1);
}

bool StrEndsWith(const char* str, const char* substr)
{
    unsigned int len_substr = strlen(substr);
    unsigned int len_str = strlen(str);

    if (len_str < len_substr) return false;

    for (unsigned int i = 0; i < len_substr; i++) {
        if (str[len_str-1-i] != substr[len_substr-1-i]) {
            return false;
        }
    }

    return true;
}
int main(int argc, char** argv) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <input mesh (.ccm/.htx)> <log2 resolution> <output path>\nYou can use the obj_to_ccm program to generate a .ccm from a .obj file\n", argv[0]);
        return 1;
    }

    srand48(123456789);

    const char* ccmPath = argv[1];
    int log2_res = atoi(argv[2]);
    if (log2_res <= 0 || log2_res > 255) {
        fprintf(stderr, "Invalid resolution\n");
        return 1;
    }
    const char* outputPath = argv[3];

    cc_Mesh* halfedge_mesh;
    if (StrEndsWith(ccmPath, ".ccm")) {
        halfedge_mesh = ccm_Load(ccmPath);
        if (!halfedge_mesh) {
            fprintf(stderr, "Failed to load .ccm\n");
            return 1;
        }
    } else if (StrEndsWith(ccmPath, ".htx")) {
        Htex::String err;
        HtexTexture* inputHtex = HtexTexture::open(ccmPath, err);
        if (inputHtex) {
            halfedge_mesh = inputHtex->getHalfedgeMesh();
        } else {
            fprintf(stderr, "Failed to load Htex file: %s\n", err.c_str());
            return 1;
        }
    } else {
        fprintf(stderr, "%s: unknown file format\n", ccmPath);
        return 1;
    }

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


    int done = 0;

    // ground truth: 1/sqrt(2) (approx. 0.707) (without cosine term)
    //               1/2 (with cosine term)
    //printf("%f\n", ComputeAO(glm::vec3(0,-0.5,0), glm::vec3(0,1,0), 1e7, scene));
    //return 0;

#pragma omp parallel default(none) shared(halfedge_mesh, width, height, scene, geom, writer, log2_res, stderr, done, halfedgeNormals) private(err)
    {
        std::vector<uint8_t> texels(width*height*4);

#pragma omp for
        for (int quadID = 0; quadID < halfedge_mesh->edgeCount; quadID++) {
            int halfedge1 = ccm_EdgeToHalfedgeID(halfedge_mesh, quadID);
            int halfedge2 = ccm_HalfedgeTwinID(halfedge_mesh, halfedge1);
            int halfedge_max = (halfedge1 > halfedge2) ? halfedge1 : halfedge2;
            int halfedge_min = (halfedge1 < halfedge2) ? halfedge1 : halfedge2;

            GenerateTexture(halfedge_mesh, halfedge_max, texels.data(), width, height, true, scene, geom, halfedgeNormals.data());
            if (halfedge_min >= 0) {
                GenerateTexture(halfedge_mesh, halfedge_min, texels.data(), width, height, false, scene, geom, halfedgeNormals.data());
            } else {
                // TODO
            }

            Htex::QuadInfo quadInfo{Htex::Res(log2_res, log2_res), quadID};

#pragma omp critical
            if (!writer->writeQuad(quadID, quadInfo, texels.data())) {
                writer->close(err);
                fprintf(stderr, "Failed to write quad %d: %s\n", quadID, err.c_str());
                exit(1);
            }

#pragma omp atomic
            done++;

            if ((done+1) % 100 == 0) {
                printf("%i / %i\n", done+1, halfedge_mesh->edgeCount);
            }
        }
    }


    if (!writer->close(err)) {
        fprintf(stderr, "Failed to write Htex file: %s\n", err.c_str());
        exit(1);
    }
}
