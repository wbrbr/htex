
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

vec4 SampleFacePoint(int halfedgeID)
{
    int edgeID = ccm_HalfedgeEdgeID(halfedgeID);
    int twinID = ccm_HalfedgeTwinID(halfedgeID);
    ivec2 uv = (halfedgeID > twinID) ? ivec2(0,0) : ivec2(1,1);

    layout(rgba8)  image2D edgeImage = layout(rgba8) image2D(u_EdgeImageHandles[edgeID]);
    ivec2 size = imageSize(edgeImage).xy;
    return imageLoad(edgeImage, uv*(size-1));
}

void StoreFacePointTexels(int halfedgeID, vec4 color, int min_res)
{
    int twinID = ccm_HalfedgeTwinID(halfedgeID);

    int edgeID = ccm_HalfedgeEdgeID(halfedgeID);
    layout(rgba8) image2D edgeImage = layout(rgba8) image2D(u_EdgeImageHandles[edgeID]);

    ivec2 log_res = u_QuadLog2Resolutions[edgeID];
    ivec2 size = ivec2(1 << log_res.x, 1 << log_res.y);

    ivec2 log_corner_size = log_res - ivec2(min_res, min_res);
    ivec2 corner_size = ivec2(1) << log_corner_size;

    /*if (ivec2(min_res) != log_res) {
        color = vec4(0,1,0,1);
    }*/

    for (int j = 0; j < corner_size.y; j++) {
        for (int i = 0; i < corner_size.x; i++) {
            ivec2 tex_coords = ivec2(i, j);
            if (halfedgeID < twinID) tex_coords = size - ivec2(1) - tex_coords;

            imageStore(edgeImage, tex_coords, color);
        }
    }

    //ivec2 uv = (halfedgeID > twinID) ? ivec2(0,0) : ivec2(1,1);
    //imageStore(edgeImage, uv*(size-1), color);
}

void main()
{
    int faceID = int(gl_GlobalInvocationID.x);
    if (faceID >= ccm_FaceCount()) {
        return;
    }

    int startHalfedge = ccm_FaceToHalfedgeID(faceID);
    int currentHalfedge = startHalfedge;

    int min_res = 1000;
    vec4 color = vec4(0);
    {
        int n = 0;
        layout(rgba8) image2D halfedgeImage;
        do {
            color += SampleFacePoint(currentHalfedge);
            n += 1;

            ivec2 log2_res = u_QuadLog2Resolutions[ccm_HalfedgeEdgeID(currentHalfedge)];
            min_res = min(min_res, min(log2_res.x, log2_res.y));

            currentHalfedge = ccm_HalfedgeNextID(currentHalfedge);
        } while (currentHalfedge != startHalfedge);

        color /= n;
    }

    currentHalfedge = startHalfedge;
    do {
        StoreFacePointTexels(currentHalfedge, color, min_res);

        currentHalfedge = ccm_HalfedgeNextID(currentHalfedge);
    } while(currentHalfedge != startHalfedge);
}
    #endif
