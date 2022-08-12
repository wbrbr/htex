layout(local_size_x = 1, local_size_y = 1, local_size_z = 1) in;

uniform vec2 u_MousePosition;
uniform vec2 u_ScreenResolution;
uniform mat4 u_ModelViewProjection;

layout(std430, binding = BUFFER_BINDING_HTEX_IMAGE_HANDLES)
readonly buffer HtexImageHandles {
    uint64_t u_HtexImageHandles[];
};

#ifdef COMPUTE_SHADER
vec3 computeBarycenter(int faceID)
{
    int startHalfedge = ccm_FaceToHalfedgeID(faceID);
    vec3 barycenter = vec3(0);
    int numVertices = 0;
    {
        int currentHEdge = startHalfedge;
        do {
            vec3 v = ccm_HalfedgeVertexPoint(currentHEdge);
            barycenter += v;
            numVertices += 1;
            currentHEdge = ccm_HalfedgeNextID(currentHEdge);
        } while (currentHEdge != startHalfedge);
    }
    barycenter /= numVertices;

    return barycenter;
}

void main() {
    int quadID = int(gl_GlobalInvocationID.x);
    ivec2 ij = ivec2(gl_GlobalInvocationID.yz);

    int offset = TEXTURETYPE_COLOR * ccm_EdgeCount();
    layout(rgba8) image2D img = layout(rgba8) image2D(u_HtexImageHandles[offset+quadID]);

    int halfedgeID = ccm_EdgeToHalfedgeID(quadID);
    int twinID = ccm_HalfedgeTwinID(halfedgeID);
    int halfedgeMax = max(halfedgeID, twinID);
    int halfedgeMin = min(halfedgeID, twinID);

    ivec2 size = imageSize(img);

    ivec2 xy;
    if (ij.x + ij.y > size.x) {
        halfedgeID = halfedgeMin;
        xy = size - 1 - ij;
    } else {
        halfedgeID = halfedgeMax;
        xy = ij;
    }
    int nextID = ccm_HalfedgeNextID(halfedgeID);

    vec2 texel_uv = (vec2(xy) + 0.5) / vec2(size);

    vec3 v0 = computeBarycenter(ccm_HalfedgeFaceID(halfedgeID));
    vec3 v1 = ccm_HalfedgeVertexPoint(halfedgeID);
    vec3 v2 = ccm_HalfedgeVertexPoint(nextID);

    //imageStore(img, ij, vec4(0,1,0,1));

    vec4 texel_pos = vec4((1-texel_uv.x-texel_uv.y)*v0 + texel_uv.x*v1 + texel_uv.y*v2, 1);
    vec4 texel_screenspace = u_ModelViewProjection * texel_pos;
    texel_screenspace /= texel_screenspace.w;

    texel_screenspace.xy = texel_screenspace.xy*0.5+0.5;
    //imageStore(img, ij, vec4(1,0,0,1));

    //imageStore(img, ij, vec4(texel_uv, 0, 1));
    //imageStore(img, ij, texel_pos);
    //imageStore(img, ij, vec4(texel_screenspace.xy, 0, 1));
    //imageStore(img, ij, vec4(u_MousePosition, 0, 1));
    vec2 aspect = vec2(1, 9.0/16.0);
    float d = distance(texel_screenspace.xy*aspect, u_MousePosition*aspect);
    float sigma = 0.03;
    float coef = exp(-d*d/(sigma*sigma));

    vec4 oldTexel = imageLoad(img, ij);
    vec4 newTexel = coef * vec4(1,0,0,1) + (1-coef) * oldTexel;
    //imageStore(img, ij, vec4(coef,coef,coef, 1));
    imageStore(img, ij, newTexel);
}
    #endif