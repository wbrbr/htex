#include <cstdlib>
#include <cstdio>
#include <vector>

#include "glad/glad.h"
#include "GLFW/glfw3.h"
#include "imgui.h"
#include "imgui_impl.h"
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"

#include "Htexture.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define DJ_OPENGL_IMPLEMENTATION
#include "dj_opengl.h"


////////////////////////////////////////////////////////////////////////////////
// Tweakable Macros
//
////////////////////////////////////////////////////////////////////////////////
#define FPS (60)
#define FRAME_COUNT (FPS * 4)

#define VIEWER_DEFAULT_WIDTH  (1920)
#define VIEWER_DEFAULT_HEIGHT (1080)

#define LOG(fmt, ...) fprintf(stdout, fmt "\n", ##__VA_ARGS__); fflush(stdout);

#ifndef PATH_TO_SRC_DIRECTORY
#   define PATH_TO_SRC_DIRECTORY "./"
#endif

#define PATH_TO_SHADER_DIRECTORY PATH_TO_SRC_DIRECTORY "shaders/"

#define BUFFER_OFFSET(i) ((char *)NULL + (i))
#define BUFFER_SIZE(x)    ((int)(sizeof(x)/sizeof(x[0])))

#define NUM_LOG2_RESOLUTIONS 12

////////////////////////////////////////////////////////////////////////////////
// Global Variables
//
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
// Window Manager
enum { AA_NONE, AA_MSAA2, AA_MSAA4, AA_MSAA8, AA_MSAA16 };
struct Window {
    GLFWwindow* handle;
    const char *name;
    int32_t width, height;
    uint32_t frameID, maxFrameID;
    bool showHud;
    struct {
        int32_t major, minor;
    } glversion;
    struct {
        uint32_t isRunning, frameID, captureID;
    } recorder;
    struct {
        int32_t aa;
        struct { int32_t fixed; } msaa;
    } framebuffer;
} g_window = {
    NULL,
    "Htex",
    VIEWER_DEFAULT_WIDTH, VIEWER_DEFAULT_HEIGHT,
    0u, ~0u,
    true,
    {4, 5},
    {0, 0, 0},
    {
        AA_NONE,
        {false}
    }
};

#define CAMERA_FOVY 55.f
// -----------------------------------------------------------------------------
// Camera Controller
struct CameraController {
    float fovy, zNear, zFar;    // perspective settings
    struct {
        float theta, phi;
    } angles;                   // axis
    float radius;               // 3D position
    float center_x, center_y;
    struct {
        float x, y, factor;
    } frameZoom;
} g_camera =
{
    CAMERA_FOVY, 0.01f, 1024.f,
    {0.5, 0},
    4.2f,
    0.f, 0.f,
    {0.f,0.f,0.f}
};

enum FilteringType {
    FILTERING_NONE,
    FILTERING_BORDERLESS,
};

enum MinMagFilter {
    MINMAG_FILTER_TRILINEAR,
    MINMAG_FILTER_LINEAR,
    MINMAG_FILTER_NEAREST,
};

enum TextureType {
    TEXTURETYPE_COLOR,
    TEXTURETYPE_DISPLACEMENT,
    TEXTURETYPE_NORMAL,
    TEXTURETYPE_AO,

    TEXTURETYPE_COUNT
};

enum ShadingMode {
    SHADING_DIFFUSE_AO,
    SHADING_DIFFUSE,
    SHADING_BASECOLOR,
    SHADING_NORMAL,
    SHADING_DISPLACEMENT,
    SHADING_AO
};

// -----------------------------------------------------------------------------
struct HtexController {
    cc_Mesh *mesh;
    int32_t tessFactor;
    std::vector<int> quadLog2Resolutions;
    enum FilteringType filter;
    MinMagFilter minmagFilter;
    struct {bool wireframe, debugTexture, displace;} flags;
    uint64_t textureFlags;
} g_htex = {
        NULL,
        0,
        {},
        FILTERING_BORDERLESS,
        MINMAG_FILTER_TRILINEAR,
        {false, false, false},
        0
};


// -----------------------------------------------------------------------------
// OpenGL resources
enum { FRAMEBUFFER_SCENE, FRAMEBUFFER_COUNT };
enum { STREAM_TRANSFORM, STREAM_COUNT };
enum { CLOCK_RENDER, CLOCK_COUNT };
enum { VERTEXARRAY_EMPTY, VERTEXARRAY_COUNT };
enum {
    BUFFER_VERTEX_TO_HALFEDGE,
    BUFFER_EDGE_TO_HALFEDGE,
    BUFFER_FACE_TO_HALFEDGE,
    BUFFER_HALFEDGES,
    BUFFER_VERTEX_POINTS,
    BUFFER_VERTEX_UVS,
    BUFFER_COUNTERS,
    BUFFER_HTEX_TEXTURE_HANDLES,
    BUFFER_HTEX_QUAD_LOG2_RESOLUTIONS,
    BUFFER_HTEX_ALPHA_TEXTURE_HANDLES,
    BUFFER_HALFEDGE_NORMALS,
    BUFFER_BARYCENTERS,
    BUFFER_FRUSTUM_PLANES,
    BUFFER_OUTPUT_VERTICES,
    BUFFER_EDGE_IMAGE_HANDLES,

    BUFFER_COUNT
};
enum {
    TEXTURE_SCENE_COLOR_BUFFER,
    TEXTURE_SCENE_DEPTH_BUFFER,
    TEXTURE_VIEWER_COLOR_BUFFER,
    TEXTURE_BLACK,

    TEXTURE_COUNT
};
enum {
    PROGRAM_VIEWER,
    PROGRAM_MAIN,
    PROGRAM_COMPUTE_HALFEDGE_NORMALS,

    PROGRAM_COUNT
};
enum {
    UNIFORM_VIEWER_SCENE_FRAMEBUFFER_SAMPLER,
    UNIFORM_VIEWER_EXPOSURE,
    UNIFORM_VIEWER_GAMMA,
    UNIFORM_VIEWER_VIEWPORT,

    UNIFORM_MVP,
    UNIFORM_TESS_FACTOR,
    UNIFORM_SCREEN_RESOLUTION,
    UNIFORM_ALPHA_TEXTURES,
    UNIFORM_DISPLACEMENT_SCALE,
    UNIFORM_DISPLACEMENT_BIAS,
    UNIFORM_MODELVIEW,
    UNIFORM_MODEL,
    UNIFORM_LOD_FACTOR,
    UNIFORM_ENABLE_ADAPTIVE_TESS,

    UNIFORM_COUNT
};

enum {
    SAMPLER_EDGE_TEXTURES,

    SAMPLER_COUNT,
};
struct OpenGLManager {
    GLuint programs[PROGRAM_COUNT];
    GLuint framebuffers[FRAMEBUFFER_COUNT];
    GLuint textures[TEXTURE_COUNT];
    GLuint vertexArrays[VERTEXARRAY_COUNT];
    GLuint buffers[BUFFER_COUNT];
    GLint uniforms[UNIFORM_COUNT];
    GLuint samplers[SAMPLER_COUNT];
    std::vector<GLuint> edgeTextures;
    std::vector<uint64_t> edgeTextureHandles;
    std::vector<uint64_t> alphaTextureHandles;

    djg_buffer *streams[STREAM_COUNT];
    djg_clock *clocks[CLOCK_COUNT];
} g_gl = {};

struct MiscState {
    float displacementScale;
    float displacementBias;
    float anisotropy;
    float edgeLength;
    ShadingMode shadingMode;
    bool disableHtex;
    bool enableAdaptiveTess;
} g_state = { 1, 0.5f, 1, 7, SHADING_BASECOLOR, false, false };

////////////////////////////////////////////////////////////////////////////////
// Utility functions
//
////////////////////////////////////////////////////////////////////////////////

float Radians(float degrees)
{
    return degrees * M_PI / 180.f;
}

static void APIENTRY
DebugOutputLogger(
    GLenum source,
    GLenum type,
    GLuint,
    GLenum severity,
    GLsizei,
    const GLchar* message,
    const GLvoid*
) {
    char srcstr[32], typestr[32];

    switch(source) {
        case GL_DEBUG_SOURCE_API: strcpy(srcstr, "OpenGL"); break;
        case GL_DEBUG_SOURCE_WINDOW_SYSTEM: strcpy(srcstr, "Windows"); break;
        case GL_DEBUG_SOURCE_SHADER_COMPILER: strcpy(srcstr, "Shader Compiler"); break;
        case GL_DEBUG_SOURCE_THIRD_PARTY: strcpy(srcstr, "Third Party"); break;
        case GL_DEBUG_SOURCE_APPLICATION: strcpy(srcstr, "Application"); break;
        case GL_DEBUG_SOURCE_OTHER: strcpy(srcstr, "Other"); break;
        default: strcpy(srcstr, "???"); break;
    };

    switch(type) {
        case GL_DEBUG_TYPE_ERROR: strcpy(typestr, "Error"); break;
        case GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR: strcpy(typestr, "Deprecated Behavior"); break;
        case GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR: strcpy(typestr, "Undefined Behavior"); break;
        case GL_DEBUG_TYPE_PORTABILITY: strcpy(typestr, "Portability"); break;
        case GL_DEBUG_TYPE_PERFORMANCE: strcpy(typestr, "Performance"); break;
        case GL_DEBUG_TYPE_OTHER: strcpy(typestr, "Message"); break;
        default: strcpy(typestr, "???"); break;
    }

    if(severity == GL_DEBUG_SEVERITY_HIGH) {
        LOG("djg_error: %s %s\n"                \
                "-- Begin -- GL_debug_output\n" \
                "%s\n"                              \
                "-- End -- GL_debug_output\n",
                srcstr, typestr, message);
    } else if(severity == GL_DEBUG_SEVERITY_MEDIUM) {
        LOG("djg_warn: %s %s\n"                 \
                "-- Begin -- GL_debug_output\n" \
                "%s\n"                              \
                "-- End -- GL_debug_output\n",
                srcstr, typestr, message);
    }
}

void InitDebugOutput(void)
{
    glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
    glDebugMessageCallback(&DebugOutputLogger, NULL);
}


////////////////////////////////////////////////////////////////////////////////
// Program Configuration
//
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
// set viewer program uniforms
void ConfigureViewerProgram()
{
    glProgramUniform1i(g_gl.programs[PROGRAM_VIEWER],
                       g_gl.uniforms[UNIFORM_VIEWER_SCENE_FRAMEBUFFER_SAMPLER],
                       TEXTURE_SCENE_COLOR_BUFFER);
}

// -----------------------------------------------------------------------------
void ConfigureMainProgram()
{
    glProgramUniform2f(g_gl.programs[PROGRAM_MAIN],
                       g_gl.uniforms[UNIFORM_SCREEN_RESOLUTION],
                       g_window.width, g_window.height);
    glProgramUniform1f(g_gl.programs[PROGRAM_MAIN],
                       g_gl.uniforms[UNIFORM_TESS_FACTOR],
                       1u << g_htex.tessFactor);
}


// -----------------------------------------------------------------------------
/**
 * Load the Viewer Program
 *
 * This program is responsible for blitting the scene framebuffer to
 * the back framebuffer, while applying gamma correction and tone mapping to
 * the rendering.
 */
bool LoadViewerProgram()
{
    djg_program *djp = djgp_create();
    GLuint *program = &g_gl.programs[PROGRAM_VIEWER];

    LOG("Loading {Viewer-Program}");
    if (g_window.framebuffer.aa >= AA_MSAA2 && g_window.framebuffer.aa <= AA_MSAA16)
        djgp_push_string(djp, "#define MSAA_FACTOR %i\n", 1 << g_window.framebuffer.aa);

    djgp_push_file(djp, PATH_TO_SHADER_DIRECTORY "Viewer.glsl");
    LOG("%s", PATH_TO_SHADER_DIRECTORY "Viewer.glsl");

    if (!djgp_to_gl(djp, 450, false, true, program)) {
        LOG("=> Failure <=\n");
        djgp_release(djp);

        return false;
    }
    djgp_release(djp);

    g_gl.uniforms[UNIFORM_VIEWER_SCENE_FRAMEBUFFER_SAMPLER] =
        glGetUniformLocation(*program, "u_FramebufferSampler");

    ConfigureViewerProgram();

    return (glGetError() == GL_NO_ERROR);
}


// -----------------------------------------------------------------------------
void djp_setup_halfedge(djg_program* djp)
{
    djgp_push_string(djp,
                     "#define CC_BUFFER_BINDING_CAGE_VERTEX_TO_HALFEDGE %i\n",
                     BUFFER_VERTEX_TO_HALFEDGE);
    djgp_push_string(djp,
                     "#define CC_BUFFER_BINDING_CAGE_EDGE_TO_HALFEDGE %i\n",
                     BUFFER_EDGE_TO_HALFEDGE);
    djgp_push_string(djp,
                     "#define CC_BUFFER_BINDING_CAGE_FACE_TO_HALFEDGE %i\n",
                     BUFFER_FACE_TO_HALFEDGE);
    djgp_push_string(djp,
                     "#define CC_BUFFER_BINDING_CAGE_HALFEDGE %i\n",
                     BUFFER_HALFEDGES);
    djgp_push_string(djp,
                     "#define CC_BUFFER_BINDING_CAGE_VERTEX_POINT %i\n",
                     BUFFER_VERTEX_POINTS);
    djgp_push_string(djp,
                     "#define CC_BUFFER_BINDING_CAGE_COUNTERS %i\n",
                     BUFFER_COUNTERS);
    djgp_push_file(djp, PATH_TO_SHADER_DIRECTORY "Halfedge.glsl");
}

bool LoadMainProgram()
{
    djg_program *djp = djgp_create();
    GLuint *glp = &g_gl.programs[PROGRAM_MAIN];

    LOG("Loading {Main-Program}");

    djgp_push_string(djp, "#extension GL_NV_gpu_shader5 : enable\n");
    djgp_push_string(djp, "#extension GL_ARB_bindless_texture : enable\n");
    djp_setup_halfedge(djp);

    djgp_push_string(djp, "#define HTEX_NUM_LOG2_RESOLUTIONS %i\n", NUM_LOG2_RESOLUTIONS);
    djgp_push_string(djp, "#define HTEX_BUFFER_BINDING_TEXTURE_HANDLES %i\n", BUFFER_HTEX_TEXTURE_HANDLES);
    djgp_push_string(djp, "#define HTEX_BUFFER_BINDING_QUAD_LOG2_RESOLUTIONS %i\n", BUFFER_HTEX_QUAD_LOG2_RESOLUTIONS);
    djgp_push_string(djp, "#define HTEX_BUFFER_BINDING_ALPHA_TEXTURE_HANDLES %i\n", BUFFER_HTEX_ALPHA_TEXTURE_HANDLES);
    djgp_push_file(djp, PATH_TO_ROOT_DIRECTORY "Htex.glsl");

    if (g_htex.flags.wireframe) {
        djgp_push_string(djp, "#define FLAG_WIRE 1\n");
    }
    if (g_htex.flags.debugTexture) {
        djgp_push_string(djp, "#define FLAG_DEBUG_TEXTURE 1\n");
    }
    if (g_htex.flags.displace) {
        djgp_push_string(djp, "#define FLAG_DISPLACE 1\n");
    }
    djgp_push_string(djp, "#define FLAG_FILTERING %i\n", g_htex.filter);
    djgp_push_string(djp, "#define FLAG_DISABLE_HTEX %i\n", g_state.disableHtex);

    djgp_push_string(djp, "#define FILTERING_NONE %i\n", FILTERING_NONE);
    djgp_push_string(djp, "#define FILTERING_BORDERLESS %i\n", FILTERING_BORDERLESS);

    djgp_push_string(djp, "#define TEXTURETYPE_COLOR %i\n", TEXTURETYPE_COLOR);
    djgp_push_string(djp, "#define TEXTURETYPE_DISPLACEMENT %i\n", TEXTURETYPE_DISPLACEMENT);
    djgp_push_string(djp, "#define TEXTURETYPE_NORMAL %i\n", TEXTURETYPE_NORMAL);
    djgp_push_string(djp, "#define TEXTURETYPE_AO %i\n", TEXTURETYPE_AO);


    djgp_push_string(djp, "#define SHADING_DIFFUSE_AO %i\n", SHADING_DIFFUSE_AO);
    djgp_push_string(djp, "#define SHADING_DIFFUSE %i\n", SHADING_DIFFUSE);
    djgp_push_string(djp, "#define SHADING_NORMAL %i\n", SHADING_NORMAL);
    djgp_push_string(djp, "#define SHADING_BASECOLOR %i\n", SHADING_BASECOLOR);
    djgp_push_string(djp, "#define SHADING_DISPLACEMENT %i\n", SHADING_DISPLACEMENT);
    djgp_push_string(djp, "#define SHADING_AO %i\n", SHADING_AO);
    djgp_push_string(djp, "#define SHADING_MODE %i\n", g_state.shadingMode);


    djgp_push_string(djp, "#define BUFFER_BINDING_HALFEDGE_NORMALS %i\n", BUFFER_HALFEDGE_NORMALS);
    djgp_push_string(djp, "#define BUFFER_BINDING_FRUSTUM_PLANES %i\n", BUFFER_FRUSTUM_PLANES);
    djgp_push_string(djp, "#define CC_BUFFER_BINDING_CAGE_VERTEX_UVS %i\n", BUFFER_VERTEX_UVS);
    djgp_push_string(djp, "#define BUFFER_BINDING_OUTPUT_VERTICES %i\n", BUFFER_OUTPUT_VERTICES);

    djgp_push_file(djp, PATH_TO_SHADER_DIRECTORY "Main.glsl");

    if (!djgp_to_gl(djp, 450, false, true, glp)) {
        djgp_release(djp);

        return false;
    }
    djgp_release(djp);

    g_gl.uniforms[UNIFORM_MVP] =
        glGetUniformLocation(*glp, "u_ModelViewProjection");
    g_gl.uniforms[UNIFORM_SCREEN_RESOLUTION] = glGetUniformLocation(*glp, "u_ScreenResolution");
    g_gl.uniforms[UNIFORM_TESS_FACTOR] = glGetUniformLocation(*glp, "u_TessFactor");
    g_gl.uniforms[UNIFORM_ALPHA_TEXTURES] = glGetUniformLocation(*glp, "u_AlphaTextures");
    g_gl.uniforms[UNIFORM_DISPLACEMENT_SCALE] = glGetUniformLocation(*glp, "u_DisplacementScale");
    g_gl.uniforms[UNIFORM_DISPLACEMENT_BIAS] = glGetUniformLocation(*glp, "u_DisplacementBias");
    g_gl.uniforms[UNIFORM_MODELVIEW] = glGetUniformLocation(*glp, "u_ModelView");
    g_gl.uniforms[UNIFORM_MODEL] = glGetUniformLocation(*glp, "u_Model");
    g_gl.uniforms[UNIFORM_LOD_FACTOR] = glGetUniformLocation(*glp, "u_LodFactor");
    g_gl.uniforms[UNIFORM_ENABLE_ADAPTIVE_TESS] = glGetUniformLocation(*glp, "u_EnableAdaptiveTess");

    ConfigureMainProgram();

    return (glGetError() == GL_NO_ERROR);
}


bool LoadHalfedgeNormalsComputeProgram()
{
    djg_program *djp = djgp_create();
    GLuint *glp = &g_gl.programs[PROGRAM_COMPUTE_HALFEDGE_NORMALS];

    LOG("Loading {ComputeHalfedgeNormals-Program}");

    djp_setup_halfedge(djp);
    djgp_push_string(djp, "#define BUFFER_BINDING_HALFEDGE_NORMALS %i\n", BUFFER_HALFEDGE_NORMALS);

    djgp_push_file(djp, PATH_TO_SHADER_DIRECTORY "ComputeHalfedgeNormals.glsl");

    if (!djgp_to_gl(djp, 450, false, true, glp)) {
        djgp_release(djp);

        return false;
    }
    djgp_release(djp);

    return (glGetError() == GL_NO_ERROR);
}

/**
 * Load All Programs
 *
 */
bool LoadPrograms()
{
    bool success = true;

    if (success) success = LoadViewerProgram();
    if (success) success = LoadMainProgram();
    if (success) success = LoadHalfedgeNormalsComputeProgram();

    return success;
}


////////////////////////////////////////////////////////////////////////////////
// Buffer Loading
//
////////////////////////////////////////////////////////////////////////////////
bool LoadBuffer(
        int32_t bufferID,
        GLsizeiptr bufferByteSize,
        const void *data,
        GLbitfield flags
) {
    GLuint* buffer = &g_gl.buffers[bufferID];
    if (glIsBuffer(*buffer))
        glDeleteBuffers(1, buffer);
    glGenBuffers(1, buffer);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, *buffer);
    glBufferStorage(GL_SHADER_STORAGE_BUFFER,
                    bufferByteSize,
                    data,
                    flags);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, bufferID, *buffer);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    return glGetError() == GL_NO_ERROR;
}

bool LoadCageVertexToHalfedgeIDsBuffer(const cc_Mesh *cage)
{
    return LoadBuffer(BUFFER_VERTEX_TO_HALFEDGE,
                      ccm_VertexCount(cage) * sizeof(int32_t),
                      cage->vertexToHalfedgeIDs,
                      0);
}

bool LoadCageEdgeToHalfedgeIDsBuffer(const cc_Mesh *cage)
{
    return LoadBuffer(BUFFER_EDGE_TO_HALFEDGE,
                      ccm_EdgeCount(cage) * sizeof(int32_t),
                      cage->edgeToHalfedgeIDs,
                      0);
}

bool LoadCageFaceToHalfedgeIDsBuffer(const cc_Mesh *cage)
{
    return LoadBuffer(BUFFER_FACE_TO_HALFEDGE,
                      ccm_FaceCount(cage) * sizeof(int32_t),
                      cage->faceToHalfedgeIDs,
                      0);
}

bool LoadCageHalfedgeBuffer(const cc_Mesh *cage)
{
    return LoadBuffer(BUFFER_HALFEDGES,
                      sizeof(cc_Halfedge) * ccm_HalfedgeCount(cage),
                      cage->halfedges,
                      0);
}

bool LoadCageVertexPointBuffer(const cc_Mesh *cage)
{
    return LoadBuffer(BUFFER_VERTEX_POINTS,
                      sizeof(cc_VertexPoint) * ccm_VertexCount(cage),
                      cage->vertexPoints,
                      0);
}

bool LoadCageVertexUvBuffer(const cc_Mesh *cage)
{
    if (ccm_UvCount(cage) != 0){
        return LoadBuffer(BUFFER_VERTEX_UVS,
                          sizeof(cc_VertexUv) * ccm_UvCount(cage),
                          cage->uvs,
                          0);
    }
    return true;
}

bool LoadCageCounterBuffer(const cc_Mesh *cage)
{
    const struct {
        int32_t vertexCount;
        int32_t halfedgeCount;
        int32_t edgeCount;
        int32_t faceCount;
        int32_t uvCount;
    } counters = {
        ccm_VertexCount(cage),
        ccm_HalfedgeCount(cage),
        ccm_EdgeCount(cage),
        ccm_FaceCount(cage),
        ccm_UvCount(cage)
    };

    return LoadBuffer(BUFFER_COUNTERS,
                      sizeof(counters),
                      &counters,
                      0);
}


bool LoadBuffers()
{
    bool success = true;

    if (success) success = LoadCageVertexToHalfedgeIDsBuffer(g_htex.mesh);
    if (success) success = LoadCageEdgeToHalfedgeIDsBuffer(g_htex.mesh);
    if (success) success = LoadCageFaceToHalfedgeIDsBuffer(g_htex.mesh);
    if (success) success = LoadCageHalfedgeBuffer(g_htex.mesh);
    if (success) success = LoadCageVertexPointBuffer(g_htex.mesh);
    if (success) success = LoadCageVertexUvBuffer(g_htex.mesh);
    if (success) success = LoadCageCounterBuffer(g_htex.mesh);
    if (success) success = LoadBuffer(BUFFER_HALFEDGE_NORMALS, g_htex.mesh->halfedgeCount * 4 * sizeof(float), nullptr,
                                      0);
    if (success) success = LoadBuffer(BUFFER_HTEX_QUAD_LOG2_RESOLUTIONS,
                                      g_htex.quadLog2Resolutions.size() * sizeof(g_htex.quadLog2Resolutions[0]),
                                      g_htex.quadLog2Resolutions.data(), 0);
    if (success) success = LoadBuffer(BUFFER_BARYCENTERS, g_htex.mesh->faceCount * 4 * sizeof(float), nullptr, 0);
    if (success) success = LoadBuffer(BUFFER_OUTPUT_VERTICES, g_htex.mesh->halfedgeCount * 3 * 64 * 4 * sizeof(float),
                                      nullptr, GL_MAP_READ_BIT);

    return success;
}

////////////////////////////////////////////////////////////////////////////////
// Texture Loading
//
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
/**
 * Load the Scene Framebuffer Textures
 *
 * Depending on the scene framebuffer AA mode, this function load 2 or
 * 3 textures. In FSAA mode, two RGBA32F and one DEPTH24_STENCIL8 textures
 * are created. In other modes, one RGBA32F and one DEPTH24_STENCIL8 textures
 * are created.
 */
bool LoadSceneFramebufferTexture()
{
    if (glIsTexture(g_gl.textures[TEXTURE_SCENE_COLOR_BUFFER]))
        glDeleteTextures(1, &g_gl.textures[TEXTURE_SCENE_COLOR_BUFFER]);
    if (glIsTexture(g_gl.textures[TEXTURE_SCENE_DEPTH_BUFFER]))
        glDeleteTextures(1, &g_gl.textures[TEXTURE_SCENE_DEPTH_BUFFER]);
    glGenTextures(1, &g_gl.textures[TEXTURE_SCENE_DEPTH_BUFFER]);
    glGenTextures(1, &g_gl.textures[TEXTURE_SCENE_COLOR_BUFFER]);

    switch (g_window.framebuffer.aa) {
        case AA_NONE:
            LOG("Loading {Scene-Z-Framebuffer-Texture}");
            glActiveTexture(GL_TEXTURE0 + TEXTURE_SCENE_DEPTH_BUFFER);
            glBindTexture(GL_TEXTURE_2D, g_gl.textures[TEXTURE_SCENE_DEPTH_BUFFER]);
            glTexStorage2D(GL_TEXTURE_2D,
                           1,
                           GL_DEPTH24_STENCIL8,
                           g_window.width,
                           g_window.height);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

            LOG("Loading {Scene-RGBA-Framebuffer-Texture}");
            glActiveTexture(GL_TEXTURE0 + TEXTURE_SCENE_COLOR_BUFFER);
            glBindTexture(GL_TEXTURE_2D, g_gl.textures[TEXTURE_SCENE_COLOR_BUFFER]);
            glTexStorage2D(GL_TEXTURE_2D,
                           1,
                           //GL_RGBA32F,
                           GL_RGBA8,
                           g_window.width,
                           g_window.height);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
            break;
        case AA_MSAA2:
        case AA_MSAA4:
        case AA_MSAA8:
        case AA_MSAA16: {
            int32_t samples = 1u << g_window.framebuffer.aa;
            int32_t maxSamples;

            glGetIntegerv(GL_MAX_INTEGER_SAMPLES, &maxSamples);
            if (samples > maxSamples) {
                LOG("note: MSAA is %ix", maxSamples);
                samples = maxSamples;
            }
            LOG("Loading {Scene-MSAA-Z-Framebuffer-Texture}");
            glActiveTexture(GL_TEXTURE0 + TEXTURE_SCENE_DEPTH_BUFFER);
            glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, g_gl.textures[TEXTURE_SCENE_DEPTH_BUFFER]);
            glTexStorage2DMultisample(GL_TEXTURE_2D_MULTISAMPLE,
                                      samples,
                                      GL_DEPTH24_STENCIL8,
                                      g_window.width,
                                      g_window.height,
                                      g_window.framebuffer.msaa.fixed);

            LOG("Loading {Scene-MSAA-RGBA-Framebuffer-Texture}");
            glActiveTexture(GL_TEXTURE0 + TEXTURE_SCENE_COLOR_BUFFER);
            glBindTexture(GL_TEXTURE_2D_MULTISAMPLE,
                          g_gl.textures[TEXTURE_SCENE_COLOR_BUFFER]);
            glTexStorage2DMultisample(GL_TEXTURE_2D_MULTISAMPLE,
                                      samples,
                                      //GL_RGBA32F,
                                      GL_RGBA8,
                                      g_window.width,
                                      g_window.height,
                                      g_window.framebuffer.msaa.fixed);
        } break;
    }
    glActiveTexture(GL_TEXTURE0);

    return (glGetError() == GL_NO_ERROR);
}


// -----------------------------------------------------------------------------
/**
 * Load the Viewer Framebuffer Texture
 *
 * This loads an RGBA8 texture used as a color buffer for the back
 * framebuffer.
 */
bool LoadViewerFramebufferTexture()
{
    LOG("Loading {Viewer-Framebuffer-Texture}");
    if (glIsTexture(g_gl.textures[TEXTURE_VIEWER_COLOR_BUFFER]))
        glDeleteTextures(1, &g_gl.textures[TEXTURE_VIEWER_COLOR_BUFFER]);
    glGenTextures(1, &g_gl.textures[TEXTURE_VIEWER_COLOR_BUFFER]);

    glActiveTexture(GL_TEXTURE0 + TEXTURE_VIEWER_COLOR_BUFFER);
    glBindTexture(GL_TEXTURE_2D, g_gl.textures[TEXTURE_VIEWER_COLOR_BUFFER]);
    glTexStorage2D(GL_TEXTURE_2D,
                   1,
                   GL_RGBA8,
                   g_window.width,
                   g_window.height);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    glActiveTexture(GL_TEXTURE0);

    return (glGetError() == GL_NO_ERROR);
}

bool LoadBlackTexture()
{
    LOG("Loading {Black-Texture}");
    if (glIsTexture(g_gl.textures[TEXTURE_BLACK]))
        glDeleteTextures(1, &g_gl.textures[TEXTURE_BLACK]);

    unsigned int& texture = g_gl.textures[TEXTURE_BLACK];
    glCreateTextures(GL_TEXTURE_2D, 1, &texture);
    glTextureStorage2D(texture, 1, GL_RGBA8, 1, 1);
    uint8_t texel[4] = { 0, 0, 0, 0};
    glTextureSubImage2D(texture,
                        0,
                        0,0,
                        1,1,
                        GL_RGBA,
                        GL_UNSIGNED_BYTE,
                        texel);

    return (glGetError() == GL_NO_ERROR);
}



// -----------------------------------------------------------------------------
/**
 * Load All Textures
 */
bool LoadTextures()
{
    bool success = true;

    if (success) success = LoadSceneFramebufferTexture();
    if (success) success = LoadViewerFramebufferTexture();
    if (success) success = LoadBlackTexture();

    return success;
}


////////////////////////////////////////////////////////////////////////////////
// Vertex Array Loading
//
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
/**
 * Load an Empty Vertex Array
 *
 * This will be used to draw procedural geometry, e.g., a fullscreen quad.
 */
bool LoadEmptyVertexArray()
{
    LOG("Loading {Empty-VertexArray}");
    if (glIsVertexArray(g_gl.vertexArrays[VERTEXARRAY_EMPTY]))
        glDeleteVertexArrays(1, &g_gl.vertexArrays[VERTEXARRAY_EMPTY]);

    glGenVertexArrays(1, &g_gl.vertexArrays[VERTEXARRAY_EMPTY]);
    glBindVertexArray(g_gl.vertexArrays[VERTEXARRAY_EMPTY]);
    glBindVertexArray(0);

    return (glGetError() == GL_NO_ERROR);
}


// -----------------------------------------------------------------------------
/**
 * Load All Vertex Arrays
 *
 */
bool LoadVertexArrays()
{
    bool success = true;

    if (success) success = LoadEmptyVertexArray();

    return success;
}


////////////////////////////////////////////////////////////////////////////////
// Framebuffer Loading
//
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
/**
 * Load the Scene Framebuffer
 *
 * This framebuffer is used to draw the 3D scene.
 * A single framebuffer is created, holding a color and Z buffer.
 * The scene writes directly to it.
 */
bool LoadSceneFramebuffer()
{
    LOG("Loading {Scene-Framebuffer}");
    if (glIsFramebuffer(g_gl.framebuffers[FRAMEBUFFER_SCENE]))
        glDeleteFramebuffers(1, &g_gl.framebuffers[FRAMEBUFFER_SCENE]);

    glGenFramebuffers(1, &g_gl.framebuffers[FRAMEBUFFER_SCENE]);
    glBindFramebuffer(GL_FRAMEBUFFER, g_gl.framebuffers[FRAMEBUFFER_SCENE]);

    if (g_window.framebuffer.aa >= AA_MSAA2 && g_window.framebuffer.aa <= AA_MSAA16) {
        glFramebufferTexture2D(GL_FRAMEBUFFER,
                               GL_COLOR_ATTACHMENT0,
                               GL_TEXTURE_2D_MULTISAMPLE,
                               g_gl.textures[TEXTURE_SCENE_COLOR_BUFFER],
                               0);
        glFramebufferTexture2D(GL_FRAMEBUFFER,
                               GL_DEPTH_STENCIL_ATTACHMENT,
                               GL_TEXTURE_2D_MULTISAMPLE,
                               g_gl.textures[TEXTURE_SCENE_DEPTH_BUFFER],
                               0);
    } else {
        glFramebufferTexture2D(GL_FRAMEBUFFER,
                               GL_COLOR_ATTACHMENT0,
                               GL_TEXTURE_2D,
                               g_gl.textures[TEXTURE_SCENE_COLOR_BUFFER],
                               0);
        glFramebufferTexture2D(GL_FRAMEBUFFER,
                               GL_DEPTH_STENCIL_ATTACHMENT,
                               GL_TEXTURE_2D,
                               g_gl.textures[TEXTURE_SCENE_DEPTH_BUFFER],
                               0);
    }

    glDrawBuffer(GL_COLOR_ATTACHMENT0);
    if (GL_FRAMEBUFFER_COMPLETE != glCheckFramebufferStatus(GL_FRAMEBUFFER)) {
        LOG("=> Failure <=");

        return false;
    }

    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    return (glGetError() == GL_NO_ERROR);
}

// -----------------------------------------------------------------------------
/**
 * Load All Framebuffers
 *
 */
bool LoadFramebuffers()
{
    bool success = true;

    if (success) success = LoadSceneFramebuffer();

    return success;
}


////////////////////////////////////////////////////////////////////////////////
// OpenGL Resource Loading
//
////////////////////////////////////////////////////////////////////////////////

struct Args {
    const char* colorFile;
    const char* displacementFile;
    const char* normalFile;
    const char* aoFile;
    const char* ccmFile;
};

bool ParseCommandLine(int argc, char **argv, Args& args)
{
    args.colorFile = nullptr;
    args.displacementFile = nullptr;
    args.normalFile = nullptr;
    args.aoFile = nullptr;
    args.ccmFile = nullptr;

    for (int i = 1; i < argc; ++i) {
        if (!strcmp("--frame-limit", argv[i])) {
            g_window.maxFrameID = atoi(argv[++i]);
            LOG("Note: frame limit set to %i\n", g_window.maxFrameID);
        } else if (!strcmp("--color", argv[i])) {
            args.colorFile = argv[++i];
        } else if (!strcmp("--disp", argv[i])) {
            args.displacementFile = argv[++i];
        } else if (!strcmp("--normal", argv[i])) {
            args.normalFile = argv[++i];
        } else if (!strcmp("--ao", argv[i])) {
            args.aoFile = argv[++i];
        } else if (!strcmp("--ccm", argv[i])) {
            args.ccmFile = argv[++i];
        }
    }

    return args.colorFile || args.displacementFile;
}


//TODO: support non-RGBA8 Htex textures

bool LoadHalfedgeTextureHandles()
{
    g_gl.edgeTextureHandles.resize(g_gl.edgeTextures.size());
    for (unsigned int i = 0; i < g_gl.edgeTextureHandles.size(); i++) {
        if (g_gl.edgeTextures[i] > 0) {
            g_gl.edgeTextureHandles[i] = glGetTextureSamplerHandleARB(g_gl.edgeTextures[i], g_gl.samplers[SAMPLER_EDGE_TEXTURES]);
        } else {
            g_gl.edgeTextureHandles[i] = glGetTextureSamplerHandleARB(g_gl.textures[TEXTURE_BLACK], g_gl.samplers[SAMPLER_EDGE_TEXTURES]);
        }
        if (!glIsTextureHandleResidentARB(g_gl.edgeTextureHandles[i])) {
            glMakeTextureHandleResidentARB(g_gl.edgeTextureHandles[i]);
        }
    }

    if (glGetError() != GL_NO_ERROR) return false;

    return LoadBuffer(BUFFER_HTEX_TEXTURE_HANDLES, g_gl.edgeTextureHandles.size() * sizeof(g_gl.edgeTextureHandles[0]),
                      g_gl.edgeTextureHandles.data(), 0);
}

bool CreateHalfedgeSampler()
{
    GLuint* sampler = &g_gl.samplers[SAMPLER_EDGE_TEXTURES];
    if (glIsSampler(*sampler)) {
        glDeleteSamplers(1, sampler);
    }
    glGenSamplers(1, sampler);
    int min_filter;
    int mag_filter;
    switch(g_htex.minmagFilter) {
        case MINMAG_FILTER_TRILINEAR:
            min_filter = GL_LINEAR_MIPMAP_LINEAR;
            mag_filter = GL_LINEAR;
            break;

        case MINMAG_FILTER_LINEAR:
            mag_filter = GL_LINEAR;
            min_filter = GL_LINEAR;
            break;

        case MINMAG_FILTER_NEAREST:
            mag_filter = GL_NEAREST;
            min_filter = GL_NEAREST;
            break;
    }
    int wrap;
    switch(g_htex.filter) {
        case FILTERING_NONE:
            wrap = GL_CLAMP_TO_EDGE;
            break;

        case FILTERING_BORDERLESS:
            wrap = GL_CLAMP_TO_BORDER;
            break;

        default:
            abort();
    }

    glSamplerParameteri(*sampler, GL_TEXTURE_MIN_FILTER, min_filter);
    glSamplerParameteri(*sampler, GL_TEXTURE_MAG_FILTER, mag_filter);
    glSamplerParameteri(*sampler, GL_TEXTURE_WRAP_S, wrap);
    glSamplerParameteri(*sampler, GL_TEXTURE_WRAP_T, wrap);
    glSamplerParameterf(*sampler, GL_TEXTURE_MAX_ANISOTROPY, g_state.anisotropy);

    return (glGetError() == GL_NO_ERROR);
}

bool PreprocessCorners(const std::vector<unsigned int>& edgeTextures, int textureType)
{
    djg_program *djp_preprocessCorners = djgp_create();
    djgp_push_string(djp_preprocessCorners, "#extension GL_NV_gpu_shader5 : enable\n");
    djgp_push_string(djp_preprocessCorners, "#extension GL_ARB_bindless_texture : enable\n");
    djp_setup_halfedge(djp_preprocessCorners);
    djgp_push_string(djp_preprocessCorners,
                     "#define BUFFER_BINDING_EDGE_IMAGE_HANDLES %i\n",
                     BUFFER_EDGE_IMAGE_HANDLES);
    djgp_push_string(djp_preprocessCorners,
                    "#define BUFFER_BINDING_HTEX_QUAD_LOG2_RESOLUTIONS %i\n",
                    BUFFER_HTEX_QUAD_LOG2_RESOLUTIONS);
    djgp_push_file(djp_preprocessCorners, PATH_TO_SHADER_DIRECTORY "PreprocessCorners.glsl");

    GLuint program_preprocessCorners;
    LOG("Loading {PreprocessCorners-Program}");
    LOG("%s", PATH_TO_SHADER_DIRECTORY "PreprocessCorners.glsl");
    if (!djgp_to_gl(djp_preprocessCorners, 450, false, true, &program_preprocessCorners)) {
        LOG("=> Failure <=\n");
        djgp_release(djp_preprocessCorners);

        return false;
    }
    djgp_release(djp_preprocessCorners);

    djg_program *djp_preprocessCornersFacePoints = djgp_create();
    djgp_push_string(djp_preprocessCornersFacePoints, "#extension GL_NV_gpu_shader5 : enable\n");
    djgp_push_string(djp_preprocessCornersFacePoints, "#extension GL_ARB_bindless_texture : enable\n");
    djp_setup_halfedge(djp_preprocessCornersFacePoints);
    djgp_push_string(djp_preprocessCornersFacePoints,
                     "#define BUFFER_BINDING_EDGE_IMAGE_HANDLES %i\n",
                     BUFFER_EDGE_IMAGE_HANDLES);
    djgp_push_string(djp_preprocessCornersFacePoints,
                    "#define BUFFER_BINDING_HTEX_QUAD_LOG2_RESOLUTIONS %i\n",
                    BUFFER_HTEX_QUAD_LOG2_RESOLUTIONS);
    djgp_push_file(djp_preprocessCornersFacePoints, PATH_TO_SHADER_DIRECTORY "PreprocessCornersFacePoints.glsl");

    GLuint program_preprocessCornersFacePoints;
    LOG("Loading {PreprocessCornersFacePoints-Program}");
    LOG("%s", PATH_TO_SHADER_DIRECTORY "PreprocessCornersFacePoints.glsl");
    if (!djgp_to_gl(djp_preprocessCornersFacePoints, 450, false, true, &program_preprocessCornersFacePoints)) {
        LOG("=> Failure <=\n");
        djgp_release(djp_preprocessCornersFacePoints);

        return false;
    }
    djgp_release(djp_preprocessCornersFacePoints);


    std::vector<uint64_t> edgeImageHandles;
    edgeImageHandles.resize(g_htex.mesh->edgeCount);

    int offset = textureType * g_htex.mesh->edgeCount;
    for (int edgeID = 0; edgeID < g_htex.mesh->edgeCount; edgeID++) {
        uint64_t handle = glGetImageHandleARB(edgeTextures[offset+edgeID], 0, false, 0, GL_RGBA8);
        edgeImageHandles[edgeID] = handle;
        if (!glIsImageHandleResidentARB(handle)) {
            glMakeImageHandleResidentARB(handle, GL_READ_WRITE);
        }
    }

    LoadBuffer(BUFFER_EDGE_IMAGE_HANDLES, edgeImageHandles.size() * sizeof(edgeImageHandles[0]),
               edgeImageHandles.data(), 0);

    glUseProgram(program_preprocessCorners);
    glDispatchCompute(g_htex.mesh->vertexCount, 1, 1);

    glUseProgram(program_preprocessCornersFacePoints);
    glDispatchCompute(g_htex.mesh->faceCount, 1, 1);

    glMemoryBarrier(GL_TEXTURE_FETCH_BARRIER_BIT);

    for (uint64_t handle : edgeImageHandles) {
        if (glIsImageHandleResidentARB(handle)) {
            glMakeImageHandleNonResidentARB(handle);
        }
    }

    return (glGetError() == GL_NO_ERROR);
}

bool LoadHtex(const char *pathToFile, int textureType)
{
    Htex::String error;
    HtexTexture *htex = HtexTexture::open(pathToFile, error);
    if (!htex) {
        fprintf(stderr, "Failed to load htex file: %s\n", error.c_str());
        return false;
    }

    const float* vertexPositions = NULL;
    int vertexCount;

    if (!g_htex.mesh) {
        g_htex.mesh = htex->getHalfedgeMesh();
    }

    int offset = textureType * g_htex.mesh->edgeCount;
    if ((int)g_gl.edgeTextures.size()>offset) {
        glDeleteTextures(g_gl.edgeTextures.size()-offset, &g_gl.edgeTextures[offset]);
    }
    g_gl.edgeTextures.resize((textureType+1) * g_htex.mesh->edgeCount);

    glCreateTextures(GL_TEXTURE_2D, g_htex.mesh->edgeCount, &g_gl.edgeTextures[offset]);

    size_t max_res = 1 << NUM_LOG2_RESOLUTIONS;
    g_htex.quadLog2Resolutions.resize(2 * g_htex.mesh->edgeCount);

    uint8_t* buf = new uint8_t[4*max_res*max_res];
    for (int edgeID = 0; edgeID < g_htex.mesh->edgeCount; edgeID++) {
        htex->getData(edgeID, buf, 0);

        Htex::Res res = htex->getQuadInfo(edgeID).res;
        g_htex.quadLog2Resolutions[2 * edgeID + 0] = res.ulog2;
        g_htex.quadLog2Resolutions[2 * edgeID + 1] = res.vlog2;

        int mipmapLevels = 1 + std::max(res.ulog2, res.vlog2);
        int width = res.u();
        int height = res.v();

        GLuint texture = g_gl.edgeTextures[offset+edgeID];
        glTextureStorage2D(texture,
                           mipmapLevels,
                           GL_RGBA8,
                           width,
                           height);
        glTextureSubImage2D(texture,
                            0,
                            0,0,
                            width,
                            height,
                            GL_RGBA,
                            GL_UNSIGNED_BYTE,
                            buf);
        glGenerateTextureMipmap(texture);
        glTextureParameteri(texture, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTextureParameteri(texture, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTextureParameteri(texture, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
        glTextureParameteri(texture, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
    }

    delete[] buf;

    return (glGetError() == GL_NO_ERROR);
}

bool CreateAlphaTextures()
{
    g_gl.alphaTextureHandles.resize(NUM_LOG2_RESOLUTIONS * NUM_LOG2_RESOLUTIONS);
    int res_v = 1;
    for (int i = 0; i < NUM_LOG2_RESOLUTIONS; i++) {
        int res_u = 1;
        for (int j = 0; j < NUM_LOG2_RESOLUTIONS; j++) {
            GLuint texture;
            glGenTextures(1, &texture);
            int idx = i * NUM_LOG2_RESOLUTIONS + j;
            glActiveTexture(GL_TEXTURE0 + TEXTURE_COUNT + idx);
            glBindTexture(GL_TEXTURE_2D, texture);

            std::vector<uint8_t> texels;
            texels.resize(res_u*res_v*4, 0xff);
            glTexStorage2D(GL_TEXTURE_2D,
                           std::max(i, j)+1,
                           GL_RGBA8,
                           res_u,
                           res_v);
            glTexSubImage2D(GL_TEXTURE_2D,
                            0,
                            0, 0,
                            res_u, res_v,
                            GL_RGBA,
                            GL_UNSIGNED_BYTE,
                            texels.data());
            glGenerateMipmap(GL_TEXTURE_2D);

            g_gl.alphaTextureHandles[idx] = glGetTextureSamplerHandleARB(texture, g_gl.samplers[SAMPLER_EDGE_TEXTURES]);
            if (!glIsTextureHandleResidentARB(g_gl.alphaTextureHandles[idx])) {
                glMakeTextureHandleResidentARB(g_gl.alphaTextureHandles[idx]);
            }

            res_u *= 2;
        }
        res_v *= 2;
    }

    LoadBuffer(BUFFER_HTEX_ALPHA_TEXTURE_HANDLES,
               g_gl.alphaTextureHandles.size() * sizeof(g_gl.alphaTextureHandles[0]),
               g_gl.alphaTextureHandles.data(),
               0);
    return (glGetError() == GL_NO_ERROR);
}

bool Load(int argc, char **argv)
{
    Args args{};
    bool success = ParseCommandLine(argc, argv, args);

    if (!success) {
        fprintf(stderr, "Usage: %s [--color <color htex>] [--disp <displacement htex>]\n", argv[0]);
    }

    if (args.colorFile) {
        g_state.shadingMode = SHADING_BASECOLOR;
    } else if (args.displacementFile) {
        g_state.shadingMode = SHADING_DISPLACEMENT;
    }

    if (!args.displacementFile) g_state.displacementBias = 0;

    InitDebugOutput();
    if (success) success = LoadTextures();
    if (success && args.colorFile) success = LoadHtex(args.colorFile, TEXTURETYPE_COLOR);
    if (success && args.displacementFile) success = LoadHtex(args.displacementFile, TEXTURETYPE_DISPLACEMENT);

    if (success) success = LoadBuffers();
    if (success) success = LoadFramebuffers();
    if (success) success = LoadPrograms();
    if (success) success = LoadVertexArrays();
    if (success) success = CreateHalfedgeSampler();
    if (success) success = CreateAlphaTextures();
    if (success && args.displacementFile) success = PreprocessCorners(g_gl.edgeTextures, TEXTURETYPE_DISPLACEMENT);
    if (success) success = LoadHalfedgeTextureHandles();

    for (int i = 0; i < CLOCK_COUNT; ++i) {
        g_gl.clocks[i] = djgc_create();
    }

    return success;
}

void Release()
{
    int32_t i;

    for (i = 0; i < STREAM_COUNT; ++i)
        if (g_gl.streams[i])
            djgb_release(g_gl.streams[i]);
    for (i = 0; i < PROGRAM_COUNT; ++i)
        if (glIsProgram(g_gl.programs[i]))
            glDeleteProgram(g_gl.programs[i]);
    for (i = 0; i < TEXTURE_COUNT; ++i)
        if (glIsTexture(g_gl.textures[i]))
            glDeleteTextures(1, &g_gl.textures[i]);
    for (i = 0; i < BUFFER_COUNT; ++i)
        if (glIsBuffer(g_gl.buffers[i]))
            glDeleteBuffers(1, &g_gl.buffers[i]);
    for (i = 0; i < FRAMEBUFFER_COUNT; ++i)
        if (glIsFramebuffer(g_gl.framebuffers[i]))
            glDeleteFramebuffers(1, &g_gl.framebuffers[i]);
    for (i = 0; i < VERTEXARRAY_COUNT; ++i)
        if (glIsVertexArray(g_gl.vertexArrays[i]))
            glDeleteVertexArrays(1, &g_gl.vertexArrays[i]);
    for (i = 0; i < CLOCK_COUNT; ++i)
        if (g_gl.clocks[i])
            djgc_release(g_gl.clocks[i]);

    ccm_Release(g_htex.mesh);
}


// -----------------------------------------------------------------------------
/**
 * Render Scene
 *
 * This procedure renders the scene to the back buffer.
 */
glm::mat4 ProjectionMatrix()
{
    const glm::mat4 P = glm::perspective(Radians(g_camera.fovy),
                                         (float)g_window.width / g_window.height,
                                         g_camera.zNear,
                                         g_camera.zFar);
    return P;
}

glm::mat4 ModelMatrix()
{
    return glm::mat4(1);
}

glm::mat4 ViewMatrix()
{
    const glm::mat4 theta = glm::rotate(glm::mat4(1.0f), Radians(g_camera.angles.theta * 180.0f - 90.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    const glm::mat4 phi = glm::rotate(glm::mat4(1.0f), Radians(g_camera.angles.phi * 360.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    const glm::mat4 translation = glm::translate(glm::mat4(1.0f), glm::vec3(g_camera.center_x, g_camera.center_y, g_camera.radius));
    const glm::mat4 viewInv = phi * theta * translation;
    const glm::mat4 V = glm::inverse(viewInv);

    return V;
}

/**
 * Extract Frustum Planes from MVP Matrix
 *
 * Based on "Fast Extraction of Viewing Frustum Planes from the World-
 * View-Projection Matrix", by Gil Gribb and Klaus Hartmann.
 * This procedure computes the planes of the frustum and normalizes
 * them.
 */
void LoadFrustum(glm::mat4 mvp, glm::vec4* planes)
{
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 2; ++j) {
            planes[i*2+j].x = mvp[0][3] + (j == 0 ? mvp[0][i] : -mvp[0][i]);
            planes[i*2+j].y = mvp[1][3] + (j == 0 ? mvp[1][i] : -mvp[1][i]);
            planes[i*2+j].z = mvp[2][3] + (j == 0 ? mvp[2][i] : -mvp[2][i]);
            planes[i*2+j].w = mvp[3][3] + (j == 0 ? mvp[3][i] : -mvp[3][i]);
            planes[i*2+j] *= length(glm::vec3(planes[i * 2 + j]));
        }
    }
}

float ComputeLodFactor()
{
    const float zoomFactor = exp2(-g_camera.frameZoom.factor);
    float tmp = (float)g_window.height
                /(2.0f * tan(glm::radians(g_camera.fovy) / 2.0f))
                / g_state.edgeLength
                / zoomFactor;
    return tmp;
}

void RenderScene()
{
    const float zoomFactor = exp2(-g_camera.frameZoom.factor);
    const float x = g_camera.frameZoom.x;
    const float y = g_camera.frameZoom.y;

    glm::mat4 zoom = glm::ortho(x - zoomFactor,
                                x + zoomFactor,
                                y - zoomFactor,
                                y + zoomFactor,
                                1.f,
                                -1.f);

    //glPatchParameteri(GL_PATCH_VERTICES, 1);
    const glm::mat4 model = ModelMatrix();
    const glm::mat4 view = ViewMatrix();
    const glm::mat4 proj =  zoom *  ProjectionMatrix();
    const glm::mat4 modelView = view * model;
    const glm::mat4 mvp = proj * modelView;
    glm::vec4 planes[6];
    LoadFrustum(mvp, planes);
    LoadBuffer(BUFFER_FRUSTUM_PLANES, sizeof(planes), planes, 0);

    djgc_start(g_gl.clocks[CLOCK_RENDER]);

    glUseProgram(g_gl.programs[PROGRAM_COMPUTE_HALFEDGE_NORMALS]);
    glDispatchCompute(g_htex.mesh->halfedgeCount, 1, 1);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    //glDisable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    glUseProgram(g_gl.programs[PROGRAM_MAIN]);
    glUniformMatrix4fv(g_gl.uniforms[UNIFORM_MVP], 1, 0, &mvp[0][0]);
    glUniformMatrix4fv(g_gl.uniforms[UNIFORM_MODELVIEW], 1, 0, &modelView[0][0]);
    glUniformMatrix4fv(g_gl.uniforms[UNIFORM_MODEL], 1, 0, &model[0][0]);
    glUniform1f(g_gl.uniforms[UNIFORM_DISPLACEMENT_SCALE], g_state.displacementScale);
    glUniform1f(g_gl.uniforms[UNIFORM_DISPLACEMENT_BIAS], g_state.displacementBias);
    glUniform1f(g_gl.uniforms[UNIFORM_LOD_FACTOR], ComputeLodFactor());
    glUniform1i(g_gl.uniforms[UNIFORM_ENABLE_ADAPTIVE_TESS], g_state.enableAdaptiveTess);
    glBindVertexArray(g_gl.vertexArrays[VERTEXARRAY_EMPTY]);
    glDrawArrays(GL_PATCHES, 0, 3 * g_htex.mesh->halfedgeCount);
    glBindVertexArray(0);
    glUseProgram(0);
    glDisable(GL_DEPTH_TEST);

    djgc_stop(g_gl.clocks[CLOCK_RENDER]);

    // enable with uniform tessellation to check for cracks
#if 0
    int tessFactor = 1 << g_htex.tessFactor;
    float* outputVertices = (float*)glMapNamedBufferRange(g_gl.buffers[BUFFER_OUTPUT_VERTICES], 0, 3*g_htex.mesh->halfedgeCount*tessFactor*4*sizeof(float), GL_MAP_READ_BIT);

    for (int halfedgeID = 0; halfedgeID < g_htex.mesh->halfedgeCount; halfedgeID++) {
        int nextID = ccm_HalfedgeNextID(g_htex.mesh, halfedgeID);
        int twinID = ccm_HalfedgeTwinID(g_htex.mesh, halfedgeID);


        // check for prev/next equality
        for (int i = 0; i < tessFactor; i++) {
            int idx = (3*halfedgeID+2)*tessFactor;
            glm::vec4 v(outputVertices[4*(idx+i)+0],
                        outputVertices[4*(idx+i)+1],
                        outputVertices[4*(idx+i)+2],
                        outputVertices[4*(idx+i)+3]);

            int idx_adj = (3*nextID+0)*tessFactor;
            glm::vec4 v_adj(
                    outputVertices[4*(idx_adj+tessFactor-i)+0],
                    outputVertices[4*(idx_adj+tessFactor-i)+1],
                    outputVertices[4*(idx_adj+tessFactor-i)+2],
                    outputVertices[4*(idx_adj+tessFactor-i)+3]);

            if (v != v_adj) {
                printf("Crack found! (prev/next): %d\n", i);
                printf("%f %f %f / %f %f %f (%g)\n", v.x,v.y,v.z, v_adj.x, v_adj.y, v_adj.z, glm::length(v-v_adj));
            }
        }

        // check for twin equality
        if (twinID >= 0){
            for (int i = 0; i < tessFactor; i++) {
                int idx = (3*halfedgeID+1)*tessFactor;
                glm::vec4 v(outputVertices[4*(idx+i)+0],
                            outputVertices[4*(idx+i)+1],
                            outputVertices[4*(idx+i)+2],
                            outputVertices[4*(idx+i)+3]);

                int idx_adj = (3*twinID+1)*tessFactor;

                glm::vec4 v_adj(
                        outputVertices[4*(idx_adj+tessFactor-i)+0],
                        outputVertices[4*(idx_adj+tessFactor-i)+1],
                        outputVertices[4*(idx_adj+tessFactor-i)+2],
                        outputVertices[4*(idx_adj+tessFactor-i)+3]);

                if (v != v_adj) {
                    printf("Crack found! (twin): %d\n", i);
                    printf("%f %f %f / %f %f %f (%g)\n", v.x,v.y,v.z, v_adj.x, v_adj.y, v_adj.z, glm::length(v_adj-v));
                }
            }
        }
    }
    glUnmapNamedBuffer(g_gl.buffers[BUFFER_OUTPUT_VERTICES]);
#endif
}

// -----------------------------------------------------------------------------
void PrintLargeNumber(const char *label, int32_t value)
{
    const int32_t M = value / 1000000;
    const int32_t K = (value - M * 1000000) / 1000;
    const int32_t r = value - M * 1000000 - K * 1000;

    if (value >= 1000000) {
        ImGui::Text("%s: %i,%03i,%03i", label, M, K, r);
    } else if (value >= 1000) {
        ImGui::Text("%s: %i,%03i", label, K, r);
    } else {
        ImGui::Text("%s: %i", label, value);
    }
}

void RenderViewer()
{
    // render framebuffer
    glUseProgram(g_gl.programs[PROGRAM_VIEWER]);
    glBindVertexArray(g_gl.vertexArrays[VERTEXARRAY_EMPTY]);
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

    // draw HUD
    if (g_window.showHud) {
        // ImGui
        glUseProgram(0);
        // Viewer Widgets
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // Performance Widget
        ImGui::SetNextWindowPos(ImVec2(g_window.width - 310, 10), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(300, 460), ImGuiCond_FirstUseEver);
        ImGui::Begin("Performances");
        {
            double cpuDt, gpuDt;

            djgc_ticks(g_gl.clocks[CLOCK_RENDER], &cpuDt, &gpuDt);
            ImGui::Text("GPU : %s", glGetString(GL_RENDERER));
            ImGui::Text("FPS : %.1f (%.3f ms)", 1.f / gpuDt, gpuDt * 1e3);
            djgc_ticks(g_gl.clocks[CLOCK_RENDER], &cpuDt, &gpuDt);
        }
        ImGui::End();

        // Camera Widgets
        ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(250, 120), ImGuiCond_FirstUseEver);
        ImGui::Begin("Camera");
        {
            const char* eAA[] = {
                "None",
                "MSAAx2",
                "MSAAx4",
                "MSAAx8",
                "MSAAx16"
            };

            if (ImGui::Combo("AA", &g_window.framebuffer.aa, &eAA[0], BUFFER_SIZE(eAA))) {
                LoadSceneFramebufferTexture();
                LoadSceneFramebuffer();
                LoadViewerProgram();
            }

            ImGui::SliderFloat("FOVY", &g_camera.fovy, 1.0f, 179.0f);
            if (ImGui::SliderFloat("zNear", &g_camera.zNear, 0.001f, 100.f)) {
                if (g_camera.zNear >= g_camera.zFar)
                    g_camera.zNear = g_camera.zFar - 0.01f;
            }
            if (ImGui::SliderFloat("zFar", &g_camera.zFar, 1.f, 1500.f)) {
                if (g_camera.zFar <= g_camera.zNear)
                    g_camera.zFar = g_camera.zNear + 0.01f;
            }

            ImGui::DragFloat("radius", &g_camera.radius, 0.1f, 0.f, 100.f);
            ImGui::DragFloat("center X", &g_camera.center_x, 0.1f, -10.f, 10.f);
            ImGui::DragFloat("center Y", &g_camera.center_y, 0.1f, -10.f, 10.f);
            ImGui::DragFloat("zoom factor", &g_camera.frameZoom.factor, 0.1f, 0.f, 10.f);
            ImGui::DragFloat("zoom x ", &g_camera.frameZoom.x, 0.01f, -10.f, 10.f);
            ImGui::DragFloat("zoom y", &g_camera.frameZoom.y, 0.01f, -10.f, 10.f);
        }
        ImGui::End();

        // Mesh Parameters
        ImGui::SetNextWindowPos(ImVec2(270, 10), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(380, 235), ImGuiCond_FirstUseEver);
        ImGui::Begin("Htex Settings");
        {
            if (ImGui::Checkbox("Wire", &g_htex.flags.wireframe)) {
                LoadPrograms();
            }
            ImGui::SameLine();
            if (ImGui::Checkbox("Displace", &g_htex.flags.displace)) {
                LoadPrograms();
            }

            if (g_htex.flags.displace) {
                ImGui::DragFloat("Displacement scale", &g_state.displacementScale, 0.1f, 0.f, 10.f);
                ImGui::DragFloat("Displacement bias", &g_state.displacementBias, 0.1f, 0, 1);
            }

            const char* filter_items[] = {
                    "None",
                    "Borderless",
            };
            if (ImGui::Combo("Filter type", (int*)&g_htex.filter, filter_items, sizeof(filter_items) / sizeof(filter_items[0]))) {
                LoadPrograms();
                CreateHalfedgeSampler();
                LoadHalfedgeTextureHandles();
                CreateAlphaTextures();
            }

            const char* intra_filter_items[] = {
                    "Trilinear",
                    "Linear",
                    "Nearest"
            };
            if (ImGui::Combo("Min/mag filter", (int*)&g_htex.minmagFilter, intra_filter_items, 3)) {
                CreateHalfedgeSampler();
                LoadHalfedgeTextureHandles();
                CreateAlphaTextures();
            }

            if (ImGui::DragFloat("Anisotropy", &g_state.anisotropy, 0.1f, 1.f, 16.f)) {
                CreateHalfedgeSampler();
                LoadHalfedgeTextureHandles();
            }

            ImGui::Checkbox("Adaptive tessellation", &g_state.enableAdaptiveTess);

            if (g_state.enableAdaptiveTess) {
                ImGui::SliderFloat("Edge length", &g_state.edgeLength, 1, 20.f);
            } else {
                if (ImGui::SliderInt("Tessellation factor", &g_htex.tessFactor, 0, 6)) {
                    ConfigureMainProgram();
                }
            }


            const char* shading_mode_items[] = {
                    "Diffuse + AO",
                    "Diffuse",
                    "BaseColor",
                    "Normal",
                    "Displacement",
                    "AO"
            };

            if (ImGui::Combo("Shading", (int*)&g_state.shadingMode, shading_mode_items, sizeof(shading_mode_items)/sizeof(shading_mode_items[0]))) {
                LoadPrograms();
            }

            if (ImGui::Checkbox("Disable Htex", &g_state.disableHtex)) {
                LoadPrograms();
            }
        }
        ImGui::End();

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    }

    // screen recording
    if (g_window.recorder.isRunning) {
        char name[64];

        glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
        sprintf(name, "capture_%02i_%09i",
                g_window.recorder.captureID,
                g_window.recorder.frameID);
        djgt_save_glcolorbuffer_bmp(GL_BACK, GL_RGB, name);
        ++g_window.recorder.frameID;
    }
}


// -----------------------------------------------------------------------------
/**
 * Blit the Composited Framebuffer to the Window Backbuffer
 *
 * Final drawing step: the composited framebuffer is blitted to the
 * OpenGL window backbuffer
 */
void Render()
{
    glBindFramebuffer(GL_FRAMEBUFFER, g_gl.framebuffers[FRAMEBUFFER_SCENE]);
    glViewport(0, 0, g_window.width, g_window.height);
    glClearColor(1,1,1,1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    RenderScene();

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glViewport(0, 0, g_window.width, g_window.height);
    glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT);
    RenderViewer();
}


////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
void
KeyboardCallback(
    GLFWwindow*,
    int key, int, int action, int
) {
    ImGuiIO& io = ImGui::GetIO();
    if (io.WantCaptureKeyboard)
        return;

    if (action == GLFW_PRESS) {
        switch (key) {
        case GLFW_KEY_ESCAPE:
            g_window.showHud = !g_window.showHud;
            break;
        case GLFW_KEY_C:
            if (g_window.recorder.isRunning) {
                g_window.recorder.frameID = 0;
                ++g_window.recorder.captureID;
            }
            g_window.recorder.isRunning = !g_window.recorder.isRunning;
            break;
        case GLFW_KEY_R:
            LoadPrograms();
            break;
        case GLFW_KEY_T: {
            char name[64], path[1024];
            static int frameID = 0;

            glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
            sprintf(name, "screenshot_%02i", frameID);
            sprintf(path, "./%s", name);
            djgt_save_glcolorbuffer_bmp(GL_FRONT, GL_RGBA, path);
            ++frameID;
            break;
        }
        case GLFW_KEY_F:
            printf("radius: %f\ntheta: %f\nphi: %f\ncenter_x: %f\ncenter_y: %f\n",
                   g_camera.radius,
                   g_camera.angles.theta,
                   g_camera.angles.phi,
                   g_camera.center_x,
                   g_camera.center_y);
            break;
        default: break;
        }
    }
}

void MouseButtonCallback(GLFWwindow*, int, int, int)
{
    ImGuiIO& io = ImGui::GetIO();

    if (io.WantCaptureMouse)
        return;
}

void MouseMotionCallback(GLFWwindow* window, double x, double y)
{
    static float x0 = 0, y0 = 0;
    const float dx = x - x0, dy = y - y0;

    ImGuiIO& io = ImGui::GetIO();
    if (io.WantCaptureMouse)
        return;

    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
        g_camera.angles.theta-= dy * 1e-3f;
        g_camera.angles.phi-= dx * 1e-3f;

        if (g_camera.angles.theta < 0.0f)
            g_camera.angles.theta = 0.0f;
        if (g_camera.angles.theta > 1.0f)
            g_camera.angles.theta = 1.0f;

    } else if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS) {
        g_camera.radius-= dx * 5e-3f;

        if (g_camera.radius < 1.0f)
            g_camera.radius = 1.0f;
    } else if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_PRESS) {
        g_camera.center_x -= dx * 1e-3f * g_camera.radius;
        g_camera.center_y += dy * 1e-3f * g_camera.radius;
    }

    x0 = x;
    y0 = y;
}

void MouseScrollCallback(GLFWwindow* window, double xoffset, double yoffset)
{
    ImGuiIO& io = ImGui::GetIO();
    ImGui_ImplGlfw_ScrollCallback(window, xoffset, yoffset);

    if (io.WantCaptureMouse)
        return;
}

void ResizeCallback(GLFWwindow* window, int width, int height)
{
    (void)window;

    g_window.width = width;
    g_window.height = height;

    LoadSceneFramebufferTexture();
    LoadSceneFramebuffer();
    ConfigureMainProgram();
}

bool CheckExtensions(const std::vector<const char*>& extensions)
{
    int numExtensions;
    glGetIntegerv(GL_NUM_EXTENSIONS, &numExtensions);
    for (const char* requestedExtension : extensions) {
        bool found = false;
        for (int i = 0; i < numExtensions; i++) {
            const char* ext = reinterpret_cast<const char*>(glGetStringi(GL_EXTENSIONS, i));
            if (strcmp(requestedExtension, ext) == 0) {
                found = true;
                break;
            }
        }
        if (!found) {
            fprintf(stderr, "Unsupported extension: %s\n", requestedExtension);
            return false;
        }
    }

    return true;
}

int main(int argc, char **argv)
{
    bool isVisible = true;

    for (int i = 0; i < argc; ++i) {
        if (!strcmp("--hidden", argv[i])) {
            isVisible = false;
        }
    }

    LOG("Loading {OpenGL Window}");
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, g_window.glversion.major);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, g_window.glversion.minor);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_VISIBLE, isVisible);
#ifndef NDEBUG
    glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GL_TRUE);
#endif
    g_window.handle = glfwCreateWindow(g_window.width,
                                       g_window.height,
                                       g_window.name,
                                       NULL, NULL);
    if (g_window.handle == NULL) {
        LOG("=> Failure <=");
        glfwTerminate();

        return -1;
    }
    glfwMakeContextCurrent(g_window.handle);
    glfwSetKeyCallback(g_window.handle, &KeyboardCallback);
    glfwSetCursorPosCallback(g_window.handle, &MouseMotionCallback);
    glfwSetMouseButtonCallback(g_window.handle, &MouseButtonCallback);
    glfwSetScrollCallback(g_window.handle, &MouseScrollCallback);
    glfwSetWindowSizeCallback(g_window.handle, &ResizeCallback);

    // load OpenGL functions
    LOG("Loading {OpenGL Functions}");
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        LOG("=> Failure <=");
        glfwTerminate();

        return -1;
    }

    // initialize ImGUI
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(g_window.handle, false);
    ImGui_ImplOpenGL3_Init("#version 450");

    if (!CheckExtensions({"GL_NV_gpu_shader5", "GL_ARB_bindless_texture"})) {
        ImGui_ImplGlfw_Shutdown();
        ImGui::DestroyContext();
        glfwTerminate();
        return -1;
    }

    // initialize
    LOG("Loading {Demo}");
    if (!Load(argc, argv)) {
        LOG("=> Failure <=");
        ImGui_ImplGlfw_Shutdown();
        ImGui::DestroyContext();
        glfwTerminate();

        return -1;
    }

    // main loop
    while (!glfwWindowShouldClose(g_window.handle) && g_window.frameID != g_window.maxFrameID) {
        glfwPollEvents();
        Render();
        glfwSwapBuffers(g_window.handle);

        if (g_window.maxFrameID != ~0u)
            ++g_window.frameID;
    }

    // cleanup
    Release();

    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwTerminate();

    return 0;
}
