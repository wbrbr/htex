
struct cc_Halfedge {
    int twinID;
    int nextID;
    int prevID;
    int faceID;
    int edgeID;
    int vertexID;
    int uvID;
};


// -----------------------------------------------------------------------------
// Buffers
#ifndef CC_BUFFER_BINDING_CAGE_VERTEX_TO_HALFEDGE
#   error Unspecified Buffer Binding
#endif
#ifndef CC_BUFFER_BINDING_CAGE_EDGE_TO_HALFEDGE
#   error Unspecified Buffer Binding
#endif
#ifndef CC_BUFFER_BINDING_CAGE_FACE_TO_HALFEDGE
#   error Unspecified Buffer Binding
#endif
#ifndef CC_BUFFER_BINDING_CAGE_HALFEDGE
#   error User must specify the binding of the cage halfedge buffer
#endif
#ifndef CC_BUFFER_BINDING_CAGE_VERTEX_POINT
#   error User must specify the binding of the cage vertex buffer
#endif
#ifndef CC_BUFFER_BINDING_CAGE_COUNTERS
#   error User must specify the binding of the cage counter
#endif

layout(std430, binding = CC_BUFFER_BINDING_CAGE_VERTEX_TO_HALFEDGE)
readonly buffer ccm_HalfedgeToVertexBuffer {
    int ccmu_VertexToHalfedgeIDs[];
};

layout(std430, binding = CC_BUFFER_BINDING_CAGE_EDGE_TO_HALFEDGE)
readonly buffer ccm_EdgeToHalfedgeBuffer {
    int ccmu_EdgeToHalfedgeIDs[];
};

layout(std430, binding = CC_BUFFER_BINDING_CAGE_FACE_TO_HALFEDGE)
readonly buffer ccm_FaceToHalfedgeBuffer {
    int ccmu_FaceToHalfedgeIDs[];
};

layout(std430, binding = CC_BUFFER_BINDING_CAGE_HALFEDGE)
readonly buffer ccm_HalfedgeBuffer {
    cc_Halfedge ccmu_Halfedges[];
};

layout(std430, binding = CC_BUFFER_BINDING_CAGE_VERTEX_POINT)
readonly buffer ccm_VertexPointBuffer {
    float ccmu_VertexPoints[];
};

layout(std430, binding = CC_BUFFER_BINDING_CAGE_COUNTERS)
readonly buffer ccm_CounterBuffer {
    int ccmu_VertexCount;
    int ccmu_HalfedgeCount;
    int ccmu_EdgeCount;
    int ccmu_FaceCount;
    int ccmu_UvCount;
};

// -----------------------------------------------------------------------------

// mesh queries
int ccm_FaceCount();
int ccm_EdgeCount();
int ccm_HalfedgeCount();
int ccm_CreaseCount();
int ccm_VertexCount();
int ccm_UvCount();

// counts at a given Catmull-Clark subdivision depth
int ccm_HalfedgeCountAtDepth(int depth);
int ccm_FaceCountAtDepth     (int depth);
int ccm_FaceCountAtDepth_Fast(int depth);
int ccm_EdgeCountAtDepth     (int depth);
int ccm_EdgeCountAtDepth_Fast(int depth);
int ccm_VertexCountAtDepth     (int depth);
int ccm_VertexCountAtDepth_Fast(int depth);

// data-access (O(1))
int ccm_HalfedgeTwinID(int halfedgeID);
int ccm_HalfedgePrevID(int halfedgeID);
int ccm_HalfedgeNextID(int halfedgeID);
int ccm_HalfedgeFaceID(int halfedgeID);
int ccm_HalfedgeEdgeID(int halfedgeID);
int ccm_HalfedgeVertexID(int halfedgeID);
int ccm_HalfedgeUvID(int halfedgeID);
float ccm_HalfedgeSharpnnes(int halfedgeID);
vec3 ccm_HalfedgeVertexPoint(int halfedgeID);
vec2 ccm_HalfedgeVertexUv(int halfedgeID);
int ccm_CreaseNextID(int edgeID);
int ccm_CreasePrevID(int edgeID);
float ccm_CreaseSharpness(int edgeID);
vec3 ccm_VertexPoint(int vertexID);
int ccm_HalfedgeNextID_Quad(int halfedgeID);
int ccm_HalfedgePrevID_Quad(int halfedgeID);
int ccm_HalfedgeFaceID_Quad(int halfedgeID);

// (vertex, edge, face) -> halfedge mappings (O(1))
int ccm_VertexToHalfedgeID(int vertexID);
int ccm_EdgeToHalfedgeID(int edgeID);
int ccm_FaceToHalfedgeID(int faceID);
int ccm_FaceToHalfedgeID_Quad(int faceID);

// halfedge remappings (O(1))
int ccm_NextVertexHalfedgeID(int halfedgeID);
int ccm_PrevVertexHalfedgeID(int halfedgeID);

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------


/*******************************************************************************
 * UV Encoding / Decoding routines
 *
 */
vec2 cc__DecodeUv(int uvEncoded)
{
    const uint tmp = uint(uvEncoded);
    const vec2 uv = vec2(
        ((tmp >>  0) & 0xFFFFu) / 65535.0f,
        ((tmp >> 16) & 0xFFFFu) / 65535.0f
    );

    return uv;
}

int cc__EncodeUv(vec2 uv)
{
    const uint u = uint(round(uv[0] * 65535.0f));
    const uint v = uint(round(uv[1] * 65535.0f));
    const uint tmp = ((u & 0xFFFFu) | ((v & 0xFFFFu) << 16));

    return int(tmp);
}


/*******************************************************************************
 * FaceCount -- Returns the number of faces
 *
 */
int ccm_FaceCount()
{
    return ccmu_FaceCount;
}


/*******************************************************************************
 * EdgeCount -- Returns the number of edges
 *
 */
int ccm_EdgeCount()
{
    return ccmu_EdgeCount;
}


/*******************************************************************************
 * HalfedgeCount -- Returns the number of half edges
 *
 */
int ccm_HalfedgeCount()
{
    return ccmu_HalfedgeCount;
}


/*******************************************************************************
 * CreaseCount -- Returns the number of creases
 *
 */
int ccm_CreaseCount()
{
    return ccm_EdgeCount();
}


/*******************************************************************************
 * VertexCount -- Returns the number of vertices
 *
 */
int ccm_VertexCount()
{
    return ccmu_VertexCount;
}


/*******************************************************************************
 * UvCount -- Returns the number of uvs
 *
 */
int ccm_UvCount()
{
    return ccmu_UvCount;
}


/*******************************************************************************
 * FaceCountAtDepth -- Returns the number of faces at a given subdivision depth
 *
 * The number of faces follows the rule
 *          F^{d+1} = H^d
 * Therefore, the number of half edges at a given subdivision depth d>= 0 is
 *          F^d = 4^{d - 1} H^0,
 * where H0 denotes the number of half-edges of the control cage.
 *
 */
int ccm_FaceCountAtDepth_Fast(int depth)
{
    const int H0 = ccm_HalfedgeCount();

    return (H0 << (2 * (depth - 1)));
}

int ccm_FaceCountAtDepth(int depth)
{
    if (depth == 0) {
        return ccm_FaceCount();
    } else {
        return ccm_FaceCountAtDepth_Fast(depth);
    }
}


/*******************************************************************************
 * EdgeCountAtDepth -- Returns the number of edges at a given subdivision depth
 *
 * The number of edges follows the rule
 *          E^{d+1} = 2 E^d + H^d
 * Therefore, the number of edges at a given subdivision depth d>= 0 is
 *          E^d = 2^{d - 1} (2 E^0 + (2^d - 1) H^0),
 * where H0 and E0 respectively denote the number of half-edges and edges
 * of the control cage.
 *
 */
int ccm_EdgeCountAtDepth_Fast(int depth)
{
    const int E0 = ccm_EdgeCount();
    const int H0 = ccm_HalfedgeCount();
    const int tmp = ~(0xFFFFFFFF << depth); // (2^d - 1)

    return ((E0 << 1) + (tmp * H0)) << (depth - 1);
}

int ccm_EdgeCountAtDepth(int depth)
{
    if (depth == 0) {
        return ccm_EdgeCount();
    } else {
        return ccm_EdgeCountAtDepth_Fast(depth);
    }
}


/*******************************************************************************
 * HalfedgeCountAtDepth -- Returns the number of half edges at a given subd depth
 *
 * The number of half edges is multiplied by 4 at each subdivision step.
 * Therefore, the number of half edges at a given subdivision depth d>= 0 is
 *          4^d H0,
 * where H0 denotes the number of half-edges of the control cage.
 *
 */
int ccm_HalfedgeCountAtDepth(int depth)
{
    const int H0 = ccm_HalfedgeCount();

    return H0 << (depth << 1);
}


/*******************************************************************************
 * CreaseCountAtDepth -- Returns the number of creases at a given subd depth
 *
 * The number of creases is multiplied by 2 at each subdivision step.
 * Therefore, the number of halfedges at a given subdivision depth d>= 0 is
 *          2^d C0,
 * where C0 denotes the number of creases of the control cage.
 *
 */
int ccm_CreaseCountAtDepth(int depth)
{
    const int C0 = ccm_CreaseCount();

    return C0 << depth;
}


/*******************************************************************************
 * VertexCountAtDepth -- Returns the number of vertices at a given subd depth
 *
 * The number of vertices follows the rule
 *          V^{d+1} = V^d + E^d + F^d
 * For a quad mesh, the number of vertices at a given subdivision depth d>= 0 is
 *          V^d = V0 + (2^{d-1} - 1)E0 + (2^{d-1} - 1)^2F0,
 * where:
 * - V0 denotes the number of vertices of the control cage
 * - E0 denotes the number of edges of the control cage
 * - F0 denotes the number of faces of the control cage
 * Note that since the input mesh may contain non-quad faces, we compute
 * the first subdivision step by hand and then apply the formula.
 *
 */
int ccm_VertexCountAtDepth_Fast(int depth)
{
    const int V0 = ccm_VertexCount();
    const int F0 = ccm_FaceCount();
    const int E0 = ccm_EdgeCount();
    const int H0 = ccm_HalfedgeCount();
    const int F1 = H0;
    const int E1 = 2 * E0 + H0;
    const int V1 = V0 + E0 + F0;
    const int tmp =  ~(0xFFFFFFFF << (depth - 1)); // 2^{d-1} - 1

    return V1 + tmp * (E1 + tmp * F1);
}

int ccm_VertexCountAtDepth(int depth)
{
    if (depth == 0) {
        return ccm_VertexCount();
    } else {
        return ccm_VertexCountAtDepth_Fast(depth);
    }
}


/*******************************************************************************
 * Halfedge data accessors
 *
 */
cc_Halfedge ccm__Halfedge(int halfedgeID)
{
    return ccmu_Halfedges[halfedgeID];
}

int ccm_HalfedgeTwinID(int halfedgeID)
{
    return ccm__Halfedge(halfedgeID).twinID;
}

int ccm_HalfedgeNextID(int halfedgeID)
{
    return ccm__Halfedge(halfedgeID).nextID;
}

int ccm_HalfedgePrevID(int halfedgeID)
{
    return ccm__Halfedge(halfedgeID).prevID;
}

int ccm_HalfedgeVertexID(int halfedgeID)
{
    return ccm__Halfedge(halfedgeID).vertexID;
}

int ccm_HalfedgeUvID(int halfedgeID)
{
    return ccm__Halfedge(halfedgeID).uvID;
}

int ccm_HalfedgeEdgeID(int halfedgeID)
{
    return ccm__Halfedge(halfedgeID).edgeID;
}

int ccm_HalfedgeFaceID(int halfedgeID)
{
    return ccm__Halfedge(halfedgeID).faceID;
}

//float ccm_HalfedgeSharpness(int halfedgeID)
//{
//    return ccm_CreaseSharpness(ccm_HalfedgeEdgeID(halfedgeID));
//}

vec3 ccm_HalfedgeVertexPoint(int halfedgeID)
{
    return ccm_VertexPoint(ccm_HalfedgeVertexID(halfedgeID));
}

int ccm_HalfedgeFaceID_Quad(int halfedgeID)
{
    return halfedgeID >> 2;
}


/*******************************************************************************
 * Halfedge Iteration (Quad-only special case)
 *
 */
int ccm__ScrollFaceHalfedgeID_Quad(int halfedgeID, int dir)
{
    const int base = 3;
    const int localID = (halfedgeID & base) + dir;

    return (halfedgeID & ~base) | (localID & base);
}

int ccm_HalfedgeNextID_Quad(int halfedgeID)
{
    return ccm__ScrollFaceHalfedgeID_Quad(halfedgeID, +1);
}

int ccm_HalfedgePrevID_Quad(int halfedgeID)
{
    return ccm__ScrollFaceHalfedgeID_Quad(halfedgeID, -1);
}


/*******************************************************************************
 * Vertex queries
 *
 */
vec3 ccm_VertexPoint(int vertexID)
{
#define vertexPoints ccmu_VertexPoints
    const float x = vertexPoints[3 * vertexID + 0];
    const float y = vertexPoints[3 * vertexID + 1];
    const float z = vertexPoints[3 * vertexID + 2];
#undef vertexPoints

    return vec3(x, y, z);
}


/*******************************************************************************
 * VertexToHalfedgeID -- Returns a half edge ID that carries a given vertex
 *
 */
int ccm_VertexToHalfedgeID(int vertexID)
{
    return ccmu_VertexToHalfedgeIDs[vertexID];
}


/*******************************************************************************
 * EdgeToHalfedgeID -- Returns a halfedge associated with a given edge
 *
 */
int ccm_EdgeToHalfedgeID(int edgeID)
{
    return ccmu_EdgeToHalfedgeIDs[edgeID];
}


/*******************************************************************************
 * FaceToHalfedgeID -- Returns a halfedge associated with a given face
 *
 */
int ccm_FaceToHalfedgeID(int faceID)
{
    return ccmu_FaceToHalfedgeIDs[faceID];
}

int ccm_FaceToHalfedgeID_Quad(int faceID)
{
    return faceID << 2;
}


/*******************************************************************************
 * Vertex Halfedge Iteration
 *
 */
int ccm_NextVertexHalfedgeID(int halfedgeID)
{
    const int twinID = ccm_HalfedgeTwinID(halfedgeID);

    return twinID >= 0 ? ccm_HalfedgeNextID(twinID) : -1;
}

int ccm_PrevVertexHalfedgeID(int halfedgeID)
{
    const int prevID = ccm_HalfedgePrevID(halfedgeID);

    return ccm_HalfedgeTwinID(prevID);
}
