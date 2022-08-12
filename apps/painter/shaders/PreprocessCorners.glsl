
#ifdef COMPUTE_SHADER
layout(local_size_x = 1, local_size_y = 1, local_size_z = 1) in;

layout(std430, binding=BUFFER_BINDING_EDGE_IMAGE_HANDLES)
readonly buffer EdgeImageHandles {
    uint64_t u_EdgeImageHandles[];
};

layout(std430, binding=BUFFER_BINDING_HTEX_QUAD_LOG2_RESOLUTIONS)
readonly buffer QuadLog2Resolutions {
    ivec2 u_QuadLog2Resolutions[];
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

void StoreVertexTexels(int halfedgeID, vec4 color, int min_res)
{
    int twinID = ccm_HalfedgeTwinID(halfedgeID);

    int edgeID = ccm_HalfedgeEdgeID(halfedgeID);
    layout(rgba8) image2D edgeImage = layout(rgba8) image2D(u_EdgeImageHandles[edgeID]);
    ivec2 log_res = u_QuadLog2Resolutions[edgeID];
    ivec2 size = ivec2(1 << log_res.x, 1 << log_res.y);

    ivec2 log_corner_size = log_res - ivec2(min_res, min_res);
    ivec2 corner_size = ivec2(1) << log_corner_size;

    for (int j = 0; j < corner_size.y; j++) {
        for (int i = 0; i < corner_size.x; i++) {
            ivec2 tex_coords = ivec2(size.x-1-i, j);
            if (halfedgeID < twinID) tex_coords = size - ivec2(1) - tex_coords;

            imageStore(edgeImage, tex_coords, color);
        }
    }
}

void main()
{
    int vertexID = int(gl_GlobalInvocationID.x);
    if (vertexID >= ccm_VertexCount()) {
        return;
    }

    int startHalfedge = ccm_VertexToHalfedgeID(vertexID);
    int currentHalfedge = startHalfedge;

    int min_res = 1000;

    vec4 color = vec4(0);
    {
        int n = 0;
        do {
            color += SampleVertexTexel(currentHalfedge);
            n += 1;

            ivec2 res = u_QuadLog2Resolutions[ccm_HalfedgeEdgeID(currentHalfedge)];
            min_res = min(min_res, min(res.x, res.y));

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

                ivec2 res = u_QuadLog2Resolutions[ccm_HalfedgeEdgeID(currentHalfedge)];
                min_res = min(min_res, min(res.x, res.y));

                prevHalfedge = ccm_HalfedgePrevID(currentHalfedge);
                currentHalfedge = ccm_HalfedgeTwinID(prevHalfedge);
            }

            int edgeID = ccm_HalfedgeEdgeID(prevHalfedge);
            layout(rgba8) image2D edgeImage = layout(rgba8) image2D(u_EdgeImageHandles[edgeID]);
            ivec2 size = imageSize(edgeImage).xy;
            color += imageLoad(edgeImage, ivec2(0,1)*(size-1));
            n += 1;

            ivec2 res = u_QuadLog2Resolutions[edgeID];
            min_res = min(min_res, min(res.x, res.y));
        }

        color /= n;
    }

    //color = vec4(0,1,0,1);
    currentHalfedge = startHalfedge;

    do {
        StoreVertexTexels(currentHalfedge, color, min_res);

        currentHalfedge = ccm_HalfedgeTwinID(currentHalfedge);
        if (currentHalfedge < 0) break;
        currentHalfedge = ccm_HalfedgeNextID(currentHalfedge);
    } while(currentHalfedge != startHalfedge);

    if (currentHalfedge < 0) {
        int prevHalfedge = ccm_HalfedgePrevID(startHalfedge);
        currentHalfedge = ccm_HalfedgeTwinID(prevHalfedge);
        while (currentHalfedge >= 0) {
            StoreVertexTexels(currentHalfedge, color, min_res);

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
