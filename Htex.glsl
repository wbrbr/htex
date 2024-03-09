layout(std430, binding = HTEX_BUFFER_BINDING_TEXTURE_HANDLES)
readonly buffer HtexTextureHandles {
    uint64_t u_HtexTextureHandles[];
};

layout(std430, binding = HTEX_BUFFER_BINDING_ALPHA_TEXTURE_HANDLES)
readonly buffer HtexAlphaTextureHandles {
    uint64_t u_HtexAlphaTextureHandles[];
};

layout(std430, binding = HTEX_BUFFER_BINDING_QUAD_LOG2_RESOLUTIONS)
readonly buffer HtexQuadLog2Resolutions {
    ivec2 u_QuadLog2Resolutions[];
};

vec4 SampleQuad(int quadID, vec2 xy, int channel)
{
    uint64_t texture_handle = u_HtexTextureHandles[channel*ccm_EdgeCount()+quadID];
    #ifdef GL_AMD_gpu_shader_int64
	sampler2D quadTexture = sampler2D(unpackUint2x32(texture_handle));
    #else
	sampler2D quadTexture = sampler2D(texture_handle);
    #endif
    return texture(quadTexture, xy);
}

float SampleAlpha(int quadID, vec2 xy)
{
    ivec2 res = u_QuadLog2Resolutions[quadID];
    uint64_t alpha_handle = u_HtexAlphaTextureHandles[res.y*HTEX_NUM_LOG2_RESOLUTIONS+res.x];
    #ifdef GL_AMD_gpu_shader_int64
	sampler2D alphaTexture = sampler2D(unpackUint2x32(alpha_handle));
    #else
	sampler2D alphaTexture = sampler2D(alpha_handle);
    #endif
    return texture(alphaTexture, xy).r;
}

vec2 TriangleToQuadUV(int halfedgeID, vec2 uv)
{
    int twinID = ccm_HalfedgeTwinID(halfedgeID);

    if (halfedgeID > twinID) {
        return uv;
    } else {
        return 1-uv;
    }
}

vec4 Htexture(int halfedgeID, vec2 uv, int channel)
{
    int nextID = ccm_HalfedgeNextID(halfedgeID);
    int prevID = ccm_HalfedgePrevID(halfedgeID);

    vec4 c = vec4(0);
    float alpha = 0;

    int quadID = ccm_HalfedgeEdgeID(halfedgeID);
    vec2 xy = TriangleToQuadUV(halfedgeID, uv);
    c += SampleQuad(quadID, xy, channel);
    alpha += SampleAlpha(quadID, xy);

    quadID = ccm_HalfedgeEdgeID(nextID);
    xy = TriangleToQuadUV(nextID, vec2(uv.y, -uv.x));
    c += SampleQuad(quadID, xy, channel);
    alpha += SampleAlpha(quadID, xy);

    quadID = ccm_HalfedgeEdgeID(prevID);
    xy = TriangleToQuadUV(prevID, vec2(-uv.y, uv.x));
    c += SampleQuad(quadID, xy, channel);
    alpha += SampleAlpha(quadID, xy);

    return c / alpha;
}
