
uniform mat4 u_ModelViewProjection = mat4(1.0f);
uniform mat4 u_ModelView;
uniform mat4 u_Model;

uniform float u_TessFactor;
uniform float u_Displacement;
uniform float u_LodFactor;
uniform float u_DisplacementBias;
uniform float u_DisplacementScale;

uniform bool u_EnableAdaptiveTess;


layout(std430, binding = BUFFER_BINDING_PTEX_HALFEDGE_NORMALS)
readonly buffer HalfedgeNormals {
    vec4 u_HalfedgeNormals[];
};

layout(std430, binding = BUFFER_BINDING_INPUT_TEXTURE_HANDLES)
readonly buffer InputTextureHandles {
    uint64_t u_InputTextureHandles[];
};

layout(std430, binding = CC_BUFFER_BINDING_CAGE_VERTEX_UVS)
readonly buffer ccm_VertexUvBuffer {
    float ccmu_VertexUvs[];
};

layout(std430, binding = BUFFER_BINDING_OUTPUT_VERTICES)
writeonly buffer ccm_OutputVertices {
    vec4 u_OutputVertices[];
};

vec2 ccm_HalfedgeVertexUv(int halfedgeID)
{
    int uvID = ccmu_Halfedges[halfedgeID].uvID;
    float u = ccmu_VertexUvs[2*uvID+0];
    float v = ccmu_VertexUvs[2*uvID+1];

    return vec2(u,v);
}

vec2 computeBarycenterUV(int faceID)
{
    int halfedgeID = ccm_FaceToHalfedgeID(faceID);
    int currentHalfedge = halfedgeID;
    vec2 sum = vec2(0);
    int n = 0;
    do {
        sum += ccm_HalfedgeVertexUv(currentHalfedge);
        n++;
        currentHalfedge = ccm_HalfedgeNextID(currentHalfedge);
    } while(currentHalfedge != halfedgeID);
    return sum / n;
}


vec2 interpolateVertexUv(vec2 halfedgeUV, int halfedgeID)
{
    vec2 uv0 = computeBarycenterUV(ccm_HalfedgeFaceID(halfedgeID));
    vec2 uv1 = ccm_HalfedgeVertexUv(halfedgeID);
    vec2 uv2 = ccm_HalfedgeVertexUv(ccm_HalfedgeNextID(halfedgeID));

    float u = halfedgeUV.x;
    float v = halfedgeUV.y;
    float w = 1-u-v;

    return w*uv0 + u*uv1 + v*uv2;
}

vec4 SampleUVMapped(int halfedgeID, vec2 halfedgeUV, int textureType)
{
    vec2 texCoords = interpolateVertexUv(halfedgeUV, halfedgeID);
    sampler2D sampler = sampler2D(u_InputTextureHandles[textureType]);
    return texture(sampler, texCoords);
}

/*******************************************************************************
* TriangleLevelOfDetail -- Computes the LoD assocaited to a triangle
*
* This function is used to garantee a user-specific pixel edge length in
* screen space. The reference edge length is that of the longest edge of the
* input triangle.In practice, we compute the LoD as:
*      LoD = 2 * log2(EdgePixelLength / TargetPixelLength)
* where the factor 2 is because the number of segments doubles every 2
* subdivision level.
*/
float TriangleLevelOfDetail_Perspective(vec3 v0, vec3 v2)
{
#if 1 //  human-readable version
    vec3 edgeCenter = (v0 + v2); // division by 2 was moved to u_LodFactor
    vec3 edgeVector = (v2 - v0);
    float distanceToEdgeSqr = dot(edgeCenter, edgeCenter);
    float edgeLengthSqr = dot(edgeVector, edgeVector);

    float distanceToEdge = sqrt(distanceToEdgeSqr);
    float edgeLength = sqrt(edgeLengthSqr);

    float tmp = edgeLength * u_LodFactor / distanceToEdge;
    return exp2(floor(log2(max(1, tmp))));
#else // optimized version
    float sqrMagSum = dot(v0, v0) + dot(v2, v2);
    float twoDotAC = 2.0f * dot(v0, v2);
    float distanceToEdgeSqr = sqrMagSum + twoDotAC;
    float edgeLengthSqr     = sqrMagSum - twoDotAC;

    return u_LodFactor + log2(edgeLengthSqr / distanceToEdgeSqr);
#endif
}

#ifdef VERTEX_SHADER
void main() {}
#endif

#ifdef TESS_CONTROL_SHADER
layout (vertices = 1) out;
patch out Triangle {
    vec4 vertexPoints[3];
    vec4 vertexNormals[3];
    int halfedgeID;
} o_Triangle;


vec3 computeVertexNormal(int vertexID)
{
    int startHalfedge = ccm_VertexToHalfedgeID(vertexID);
    vec3 n = vec3(0);
    int currentHEdge = startHalfedge;
    do {
        n += u_HalfedgeNormals[currentHEdge].xyz;
        currentHEdge = ccm_HalfedgeTwinID(currentHEdge);
        if (currentHEdge < 0) break;
        currentHEdge = ccm_HalfedgeNextID(currentHEdge);
    } while (currentHEdge != startHalfedge);

    // boundary, we do a backwards traversal
    if (currentHEdge != startHalfedge) {
        currentHEdge = ccm_HalfedgeTwinID(ccm_HalfedgePrevID(startHalfedge));
        while (currentHEdge >= 0) {
            n += u_HalfedgeNormals[currentHEdge].xyz;
            currentHEdge = ccm_HalfedgeTwinID(ccm_HalfedgePrevID(currentHEdge));
        }
    }

    return normalize(n);
}

vec3 computeFacePointNormal(int faceID)
{
    int startHalfedge = ccm_FaceToHalfedgeID(faceID);
    vec3 n = vec3(0);
    int currentHEdge = startHalfedge;
    do {
        n += computeVertexNormal(ccm_HalfedgeVertexID(currentHEdge)).xyz;
        currentHEdge = ccm_HalfedgeNextID(currentHEdge);
    } while(currentHEdge != startHalfedge);
    return normalize(n);
}

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

void main()
{
    int halfedgeID = int(gl_PrimitiveID);

    o_Triangle.halfedgeID = halfedgeID;

    int nextID = ccm_HalfedgeNextID(halfedgeID);

    o_Triangle.vertexNormals[0] = vec4(computeFacePointNormal(ccm_HalfedgeFaceID(halfedgeID)), 0);
    o_Triangle.vertexNormals[1] = vec4(computeVertexNormal(ccm_HalfedgeVertexID(halfedgeID)), 0);
    o_Triangle.vertexNormals[2] = vec4(computeVertexNormal(ccm_HalfedgeVertexID(nextID)), 0);

    vec3 v0 = computeBarycenter(ccm_HalfedgeFaceID(halfedgeID));
    vec3 v1 = ccm_HalfedgeVertexPoint(halfedgeID);
    vec3 v2 = ccm_HalfedgeVertexPoint(nextID);

    int faceID = ccm_HalfedgeFaceID(halfedgeID);

    o_Triangle.vertexPoints[0] = vec4(v0, 1);
    o_Triangle.vertexPoints[1] = vec4(v1, 1);
    o_Triangle.vertexPoints[2] = vec4(v2, 1);

    if (u_EnableAdaptiveTess) {
        vec4 viewSpaceVertices[3]; for (int i = 0; i < 3; i++) {
            viewSpaceVertices[i] = u_ModelView * o_Triangle.vertexPoints[i];
        }

        gl_TessLevelOuter[0] = TriangleLevelOfDetail_Perspective(viewSpaceVertices[2].xyz, viewSpaceVertices[0].xyz);
        gl_TessLevelOuter[1] = TriangleLevelOfDetail_Perspective(viewSpaceVertices[0].xyz, viewSpaceVertices[1].xyz);
        gl_TessLevelOuter[2] = TriangleLevelOfDetail_Perspective(viewSpaceVertices[1].xyz, viewSpaceVertices[2].xyz);
    } else {
        gl_TessLevelOuter[0] =
        gl_TessLevelOuter[1] =
        gl_TessLevelOuter[2] = u_TessFactor;
    }

    gl_TessLevelInner[0] = min(min(gl_TessLevelOuter[0], gl_TessLevelOuter[1]), gl_TessLevelOuter[2]);
}
#endif

#ifdef TESS_EVALUATION_SHADER
layout (triangles, ccw, equal_spacing) in;
//layout (triangles, ccw, fractional_even_spacing) in;
patch in Triangle {
    vec4 vertexPoints[3];
    vec4 vertexNormals[3];
    int halfedgeID;
} i_Patch;

layout(location = 0) flat out int o_HalfedgeID;
layout(location = 1) out vec2 o_VertexUv;
layout(location = 2) out vec3 o_VertexPosition;

void main()
{
    vec4 vertices[3] = i_Patch.vertexPoints;
    vec4 normals[3] = i_Patch.vertexNormals;

    float u = gl_TessCoord.x;
    float v = gl_TessCoord.y;
    float w = gl_TessCoord.z;

    #if 0
    vec4 vertexPoint = orderIndependentMAD(vertices, w, u, v);
    vec3 vertexNormal = orderIndependentMAD(normals, w, u, v).xyz;
    #else
    precise vec4 vertexPoint = w*vertices[0] + u*vertices[1] + v*vertices[2];
    precise vec3 vertexNormal = (normals[0]*w + normals[1]*u + normals[2]*v).xyz;
    #endif
    vertexNormal = normalize(vertexNormal);

    float displacement = 0;
    #if FLAG_DISPLACE

    displacement = Htexture(i_Patch.halfedgeID, gl_TessCoord.xy, TEXTURETYPE_DISPLACEMENT).r;
    vertexPoint.xyz = vertexPoint.xyz + vertexNormal * (displacement-u_DisplacementBias) * u_DisplacementScale;
    #endif

    if (u == 0) {
        // edge 2-0
        int wi = int(round(w*u_TessFactor));
        int idx = (3*i_Patch.halfedgeID+2) * int(u_TessFactor)+ wi;
        //u_OutputVertices[idx] = vertexPoint;
        u_OutputVertices[idx] = vec4(displacement);
    } else if (v == 0) {
        // edge 0-1
        int ui = int(round(u*u_TessFactor));
        int idx = (3*i_Patch.halfedgeID+0) * int(u_TessFactor) + ui;
        //u_OutputVertices[idx] = vertexPoint;
        u_OutputVertices[idx] = vec4(displacement);
    } else if (w == 0) {
        // edge 1-2
        int vi = int(round(v*u_TessFactor));
        int idx = (3*i_Patch.halfedgeID+1) * int(u_TessFactor) + vi;
        //u_OutputVertices[idx] = vertexPoint;
        u_OutputVertices[idx] = vec4(displacement);
    }

    gl_Position = u_ModelViewProjection * vertexPoint;
    o_VertexUv = gl_TessCoord.xy;
    o_HalfedgeID = i_Patch.halfedgeID;
    o_VertexPosition = vertexPoint.xyz;
}
#endif

#ifdef GEOMETRY_SHADER
layout (triangles) in;
layout (triangle_strip, max_vertices = 3) out;
layout (location = 0) flat in int i_HalfedgeID[];
layout(location = 1) in vec2 i_VertexUV[];
layout(location = 2) in vec3 i_VertexPosition[];
layout (location = 0) flat out int o_HalfedgeID;
layout(location = 1) flat out vec2 o_VertexUV;
layout (location = 2) noperspective out vec3 o_Distance;
layout(location = 3) out vec3 o_VertexNormal;

uniform vec2 u_ScreenResolution;

void main()
{
    vec2 p0 = u_ScreenResolution * gl_in[0].gl_Position.xy / gl_in[0].gl_Position.w;
    vec2 p1 = u_ScreenResolution * gl_in[1].gl_Position.xy / gl_in[1].gl_Position.w;
    vec2 p2 = u_ScreenResolution * gl_in[2].gl_Position.xy / gl_in[2].gl_Position.w;
    vec2 v[3] = vec2[3](p2 - p1, p2 - p0, p1 - p0);
    float area = abs(v[1].x * v[2].y - v[1].y * v[2].x);

    vec3 n = normalize(cross(i_VertexPosition[1]-i_VertexPosition[0], i_VertexPosition[2] - i_VertexPosition[0]));

    for (int i = 0; i < 3; ++i) {
        o_HalfedgeID = i_HalfedgeID[i];
        o_VertexUV = i_VertexUV[i];
        o_VertexNormal = n;
        o_Distance = vec3(0);
        o_Distance[i] = area * inversesqrt(dot(v[i],v[i]));
        gl_Position = gl_in[i].gl_Position;
        EmitVertex();
    }

    EndPrimitive();
}
#endif

#ifdef FRAGMENT_SHADER
layout (location = 0) flat in int i_HalfedgeID;
layout(location = 1) in vec2 i_FragmentUV;
layout (location = 2) noperspective in vec3 i_Distance;
layout(location = 3) in vec3 i_FragmentNormal;
layout (location = 0) out vec4 o_FragmentColor;

void main()
{
#if SHADING_MODE == SHADING_DIFFUSE_AO
    vec3 objectSpaceNormal = PtexTexture(i_FragmentUV, i_HalfedgeID, TEXTURETYPE_NORMAL).xyz * 2 -1;
    vec3 worldSpaceNormal = normalize(vec3(u_Model * vec4(objectSpaceNormal, 0)));

    vec3 L = normalize(vec3(0.3, 1.0, 0.3));
    float val = clamp(dot(worldSpaceNormal, L), 0, 1) * 0.5;

    float ao = clamp(PtexTexture(i_FragmentUV, i_HalfedgeID, TEXTURETYPE_AO).r, 0, 1);
    vec3 ambient = vec3(0.5) * ao;
    vec4 color = vec4(vec3(ambient + val), 1);
#elif SHADING_MODE == SHADING_DIFFUSE
    vec3 objectSpaceNormal = PtexTexture(i_FragmentUV, i_HalfedgeID, TEXTURETYPE_NORMAL).xyz * 2 -1;
    vec3 worldSpaceNormal = normalize(vec3(u_Model * vec4(objectSpaceNormal, 0)));
    vec3 L = normalize(vec3(0.4, 1.0, 1.0));

    float val = clamp(dot(worldSpaceNormal, L), 0, 1);
    float ambient = 1;
    ambient = 0;

    /* ambient = 0;
    val = dot(worldSpaceNormal, L) * 0.5 + 0.5; */
    vec4 color = vec4(vec3(val*0.8 + ambient*0.15), 1);
#elif SHADING_MODE == SHADING_NORMAL
    vec4 color = PtexTexture(i_FragmentUV, i_HalfedgeID, TEXTURETYPE_NORMAL);
#elif SHADING_MODE == SHADING_BASECOLOR
    vec4 color = Htexture(i_HalfedgeID, i_FragmentUV, TEXTURETYPE_COLOR);
#elif SHADING_MODE == SHADING_DISPLACEMENT
    vec4 color = Htexture(i_HalfedgeID, i_FragmentUV, TEXTURETYPE_DISPLACEMENT);
#elif SHADING_MODE == SHADING_AO
    vec4 color = clamp(PtexTexture(i_FragmentUV, i_HalfedgeID, TEXTURETYPE_AO), 0, 1);
#endif
//color = vec4(0.5);

    vec2 texCoords = interpolateVertexUv(i_FragmentUV, i_HalfedgeID);
    //color = vec4(texCoords, 0, 1);

    //color = vec4(abs(i_FragmentNormal), 1);

#if FLAG_WIRE
    const float wireScale = 1.0; // scale of the wire in pixel
    vec4 wireColor = vec4(0.0, 0.0, 0.0, 1.0);
    vec3 distanceSquared = i_Distance * i_Distance;
    float nearestDistance = min(min(distanceSquared.x, distanceSquared.y), distanceSquared.z);
    float blendFactor = exp2(-nearestDistance / wireScale);

    color = mix(color, wireColor, blendFactor);
#endif

    o_FragmentColor = color;
}
#endif
