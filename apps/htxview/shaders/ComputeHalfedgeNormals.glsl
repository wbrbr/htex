
layout(local_size_x = 1, local_size_y = 1, local_size_z = 1) in;

layout(std430, binding = BUFFER_BINDING_PTEX_HALFEDGE_NORMALS)
writeonly buffer HalfedgeNormals {
    vec4 u_HalfedgeNormals[];
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

void main()
{
    int halfedgeID = int(gl_GlobalInvocationID.x);

    vec3 v0 = computeBarycenter(ccm_HalfedgeFaceID(halfedgeID));
    vec3 v1 = ccm_HalfedgeVertexPoint(halfedgeID);
    vec3 v2 = ccm_HalfedgeVertexPoint(ccm_HalfedgeNextID(halfedgeID));

    vec3 n = normalize(cross(v1-v0, v2-v0));

    u_HalfedgeNormals[halfedgeID] = vec4(n, 0);
}
#endif
