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
    sampler2D quadTexture = sampler2D(u_HtexTextureHandles[channel*ccm_EdgeCount()+quadID]);
    return texture(quadTexture, xy);
}

float SampleAlpha(int quadID, vec2 xy)
{
    ivec2 res = u_QuadLog2Resolutions[quadID];
    sampler2D alphaTexture = sampler2D(u_HtexAlphaTextureHandles[res.y*HTEX_NUM_LOG2_RESOLUTIONS+res.x]);
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

/// quadIDs contains the ID of the quad for the current, next and prev halfedges
vec4 Htexture_sample(ivec3 quadIDs, vec2 xy, vec2 xyNext, vec2 xyPrev, int channel)
{
    vec4 c = vec4(0);
    float alpha = 0;

    c += SampleQuad(quadIDs.x, xy, channel);
    alpha += SampleAlpha(quadIDs.x, xy);

    c += SampleQuad(quadIDs.y, xyNext, channel);
    alpha += SampleAlpha(quadIDs.y, xy);

    c += SampleQuad(quadIDs.z, xyPrev, channel);
    alpha += SampleAlpha(quadIDs.z, xy);

    return c / alpha;
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
