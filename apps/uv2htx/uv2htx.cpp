#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "Htexture.h"
#include "glm/glm.hpp"
#include "CatmullClark.h"

struct Texture {
    unsigned char* texels;
    int width;
    int height;
};

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

Color sampleTexture(glm::vec2 uv, const Texture& texture)
{
    glm::vec2 unnorm = uv * glm::vec2(texture.width, texture.height) - glm::vec2(0.5f);
    int x = (int)unnorm.x;
    int y = (int)unnorm.y;

    float tx = unnorm.x - (float)x;
    float ty = unnorm.y - (float)y;

    int x_left = max(x, 0);
    int y_bot = max(y, 0);
    int x_right = min(x+1, texture.width-1);
    int y_top = min(y+1, texture.height-1);

    int px_botleft = y_bot*texture.width+x_left;
    int px_botright = y_bot*texture.width+x_right;
    int px_topleft = y_top*texture.width+x_left;
    int px_topright = y_top*texture.width+x_right;

    glm::vec4 botLeft(texture.texels[4*px_botleft+0], texture.texels[4*px_botleft+1], texture.texels[4*px_botleft+2], texture.texels[4*px_botleft+3]);
    glm::vec4 botRight(texture.texels[4*px_botright+0], texture.texels[4*px_botright+1], texture.texels[4*px_botright+2], texture.texels[4*px_botright+3]);
    glm::vec4 topLeft(texture.texels[4*px_topleft+0], texture.texels[4*px_topleft+1], texture.texels[4*px_topleft+2], texture.texels[4*px_topleft+3]);
    glm::vec4 topRight(texture.texels[4*px_topright+0], texture.texels[4*px_topright+1], texture.texels[4*px_topright+2], texture.texels[4*px_topright+3]);

    glm::vec4 result = botLeft*(1-tx)*(1-ty) + botRight*tx*(1-ty) + topLeft*(1-tx)*ty + topRight*tx*ty;

    return {
            (uint8_t)result.x,
            (uint8_t)result.y,
            (uint8_t)result.z,
            (uint8_t)result.w
    };
}

glm::vec2 VertexUvToVec2(cc_VertexUv uv)
{
    return glm::vec2(uv.u, uv.v);
}

glm::vec2 FacePointUV(const cc_Mesh* mesh, int faceID)
{
    int startHalfedge = ccm_FaceToHalfedgeID(mesh, faceID);
    int n = 0;
    glm::vec2 uv(0);

    int halfedgeID = startHalfedge;
    do {
        cc_VertexUv vertexUV = ccm_HalfedgeVertexUv(mesh, halfedgeID);
        uv += glm::vec2(vertexUV.u, vertexUV.v);
        n++;

        halfedgeID = ccm_HalfedgeNextID(mesh, halfedgeID);
    } while (halfedgeID != startHalfedge);

    uv /= n;
    return uv;
}

glm::vec2 BarycentricInterpolation(float u, float v, glm::vec2 v0, glm::vec2 v1, glm::vec2 v2)
{
    return (1-u-v)*v0 + u*v1 + v*v2;
}

void generateTexture(const cc_Mesh* mesh, int quadID, uint8_t* texels, int width, int height, const Texture& input_texture) {
    int halfedge1 = ccm_EdgeToHalfedgeID(mesh, quadID);
    int halfedge2 = ccm_HalfedgeTwinID(mesh, halfedge1);

    int halfedgeMax = halfedge1 > halfedge2 ? halfedge1 : halfedge2;
    int halfedgeMin = halfedge1 < halfedge2 ? halfedge1 : halfedge2;

    glm::vec2 max_uv0 = FacePointUV(mesh, ccm_HalfedgeFaceID(mesh, halfedgeMax));
    glm::vec2 max_uv1 = VertexUvToVec2(ccm_HalfedgeVertexUv(mesh, halfedgeMax));
    glm::vec2 max_uv2 = VertexUvToVec2(ccm_HalfedgeVertexUv(mesh, ccm_HalfedgeNextID(mesh, halfedgeMax)));

    glm::vec2 min_uv0(-1);
    glm::vec2 min_uv1(-1);
    glm::vec2 min_uv2(-1);
    if (halfedgeMin >= 0) {
        min_uv0 = FacePointUV(mesh, ccm_HalfedgeFaceID(mesh, halfedgeMin));
        min_uv1 = VertexUvToVec2(ccm_HalfedgeVertexUv(mesh, halfedgeMin));
        min_uv2 = VertexUvToVec2(ccm_HalfedgeVertexUv(mesh, ccm_HalfedgeNextID(mesh, halfedgeMin)));
    }

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            float u = ((float)j + .5f) / (float)width;
            float v = ((float)i + .5f) / (float)height;

            int px = i*width+j;

            glm::vec2 input_texcoords(0);
            if (u+v <= 1.f) {
                input_texcoords = BarycentricInterpolation(u, v, max_uv0, max_uv1, max_uv2);
            } else {
                if (halfedgeMin >= 0) {
                    input_texcoords = BarycentricInterpolation(1-u, 1-v, min_uv0, min_uv1, min_uv2);
                } else {
                    input_texcoords = BarycentricInterpolation(1-u, 1-v, max_uv0, max_uv2, max_uv1);
                }
            }

            Color color = sampleTexture(input_texcoords, input_texture);

            texels[4*px+0] = color.r;
            texels[4*px+1] = color.g;
            texels[4*px+2] = color.b;
            texels[4*px+3] = color.a;
        }
    }
}


int main(int argc, char** argv) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <ccm file> <texture file> <log2 resolution> <output path>\n", argv[0]);
        return 1;
    }

    const char* ccmPath = argv[1];
    const char* inputTexturePath = argv[2];
    int log2_res = atoi(argv[3]);
    if (log2_res <= 0 || log2_res > 255) {
        fprintf(stderr, "Invalid resolution\n");
        return 1;
    }
    const char* outputPath = argv[4];

    cc_Mesh* halfedge_mesh = ccm_Load(ccmPath);

    Texture input_texture{};
    {
        int width, height, channels;
        stbi_set_flip_vertically_on_load(true);
        unsigned char* input_texels = stbi_load(inputTexturePath, &width, &height, &channels, 4);
        input_texture.texels = input_texels;
        input_texture.width = width;
        input_texture.height = height;
    }

    Htex::String err;
    HtexWriter* writer = HtexWriter::open(outputPath, halfedge_mesh, Htex::mt_quad, Htex::dt_uint8, 4, 3, err);
    if (!writer) {
        fprintf(stderr, "Failed to create HtexWriter: %s\n", err.c_str());
        exit(1);
    }

    const int width = 1 << log2_res;
    const int height = 1 << log2_res;

    uint8_t* texels = new uint8_t[width*height*4];
    for (int quadID = 0; quadID < halfedge_mesh->edgeCount; quadID++) {
        generateTexture(halfedge_mesh, quadID, texels, width, height, input_texture);
        Htex::QuadInfo faceInfo{Htex::Res(log2_res, log2_res), quadID};
        if (!writer->writeQuad(quadID, faceInfo, texels)) {
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
