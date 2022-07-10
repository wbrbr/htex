
#ifdef COMPUTE_SHADER
layout(local_size_x = 1, local_size_y = 1, local_size_z = 1) in;

layout(std430, binding=BUFFER_BINDING_EDGE_IMAGE_HANDLES)
readonly buffer EdgeImageHandles {
    uint64_t u_EdgeImageHandles[];
};

vec4 SampleVertexTexel(int halfedgeID)
{
    int twinID = ccm_HalfedgeTwinID(halfedgeID);
    ivec2 uv = (halfedgeID > twinID) ? ivec2(1,0) : ivec2(0,1);

    int edgeID = ccm_HalfedgeEdgeID(halfedgeID);

    layout(rgba8) image2D edgeImage = layout(rgba8) image2D(u_EdgeImageHandles[edgeID]);
    ivec2 size = imageSize(edgeImage).xy;

    return imageLoad(edgeImage, uv * (size-1));
}

void StoreVertexTexel(int halfedgeID, vec4 color)
{
    int twinID = ccm_HalfedgeTwinID(halfedgeID);
    ivec2 uv = (halfedgeID > twinID) ? ivec2(1,0) : ivec2(0,1);

    int edgeID = ccm_HalfedgeEdgeID(halfedgeID);
    layout(rgba8) image2D edgeImage = layout(rgba8) image2D(u_EdgeImageHandles[edgeID]);
    ivec2 size = imageSize(edgeImage).xy;

    imageStore(edgeImage, uv*(size-1), color);
}

void main()
{
    int vertexID = int(gl_GlobalInvocationID.x);
    if (vertexID >= ccm_VertexCount()) {
        return;
    }

    int startHalfedge = ccm_VertexToHalfedgeID(vertexID);
    int currentHalfedge = startHalfedge;

    vec4 color = vec4(0);
    {
        int n = 0;
        do {
            color += SampleVertexTexel(currentHalfedge);
            n += 1;

            currentHalfedge = ccm_HalfedgeTwinID(currentHalfedge);
            if (currentHalfedge < 0) break;
            currentHalfedge = ccm_HalfedgeNextID(currentHalfedge);
        } while (currentHalfedge != startHalfedge);

        if (currentHalfedge < 0) {
            int prevHalfedge = ccm_HalfedgePrevID(startHalfedge);
            currentHalfedge = ccm_HalfedgeTwinID(prevHalfedge);

            while (currentHalfedge >= 0) {
                color += SampleVertexTexel(currentHalfedge);
                n += 1;

                prevHalfedge = ccm_HalfedgePrevID(currentHalfedge);
                currentHalfedge = ccm_HalfedgeTwinID(prevHalfedge);
            }

            int edgeID = ccm_HalfedgeEdgeID(prevHalfedge);
            layout(rgba8) image2D edgeImage = layout(rgba8) image2D(u_EdgeImageHandles[edgeID]);
            ivec2 size = imageSize(edgeImage).xy;
            color += imageLoad(edgeImage, ivec2(0,1)*(size-1));
            n += 1;
        }

        color /= n;
    }

    //color = vec4(0,1,0,1);
    currentHalfedge = startHalfedge;

    do {
        StoreVertexTexel(currentHalfedge, color);

        currentHalfedge = ccm_HalfedgeTwinID(currentHalfedge);
        if (currentHalfedge < 0) break;
        currentHalfedge = ccm_HalfedgeNextID(currentHalfedge);
    } while(currentHalfedge != startHalfedge);

    if (currentHalfedge < 0) {
        int prevHalfedge = ccm_HalfedgePrevID(startHalfedge);
        currentHalfedge = ccm_HalfedgeTwinID(prevHalfedge);
        while (currentHalfedge >= 0) {
            StoreVertexTexel(currentHalfedge, color);

            prevHalfedge = ccm_HalfedgePrevID(currentHalfedge);
            currentHalfedge = ccm_HalfedgeTwinID(prevHalfedge);
        }

        int edgeID = ccm_HalfedgeEdgeID(prevHalfedge);
        layout(rgba8) image2D edgeImage = layout(rgba8) image2D(u_EdgeImageHandles[edgeID]);
        ivec2 size = imageSize(edgeImage).xy;
        imageStore(edgeImage, ivec2(0,1)*(size-1), color);
    }
}
#endif