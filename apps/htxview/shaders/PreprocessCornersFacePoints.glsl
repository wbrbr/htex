
#ifdef COMPUTE_SHADER
layout(local_size_x = 1, local_size_y = 1, local_size_z = 1) in;

layout(std430, binding=BUFFER_BINDING_EDGE_IMAGE_HANDLES)
readonly buffer EdgeImageHandles {
    uint64_t u_EdgeImageHandles[];
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

void StoreFacePoint(int halfedgeID, vec4 color)
{
    int edgeID = ccm_HalfedgeEdgeID(halfedgeID);
    int twinID = ccm_HalfedgeTwinID(halfedgeID);
    ivec2 uv = (halfedgeID > twinID) ? ivec2(0,0) : ivec2(1,1);

    layout(rgba8)  image2D edgeImage = layout(rgba8) image2D(u_EdgeImageHandles[edgeID]);
    ivec2 size = imageSize(edgeImage).xy;
    imageStore(edgeImage, uv*(size-1), color);
}

void main()
{
    int faceID = int(gl_GlobalInvocationID.x);
    if (faceID >= ccm_FaceCount()) {
        return;
    }

    int startHalfedge = ccm_FaceToHalfedgeID(faceID);
    int currentHalfedge = startHalfedge;

    vec4 color = vec4(0);
    {
        int n = 0;
        layout(rgba8) image2D halfedgeImage;
        do {
            color += SampleFacePoint(currentHalfedge);
            n += 1;

            currentHalfedge = ccm_HalfedgeNextID(currentHalfedge);
        } while (currentHalfedge != startHalfedge);

        color /= n;
    }

    currentHalfedge = startHalfedge;
    do {
        StoreFacePoint(currentHalfedge, color);

        currentHalfedge = ccm_HalfedgeNextID(currentHalfedge);
    } while(currentHalfedge != startHalfedge);
}
    #endif
