#include <cstdio>
#include "Ptexture.h"

int main(int argc,char** argv) {
    if (argc != 2) {
        fprintf(stderr, "invalid usage\n");
        return 1;
    }

    cc_Mesh* mesh = ccm_Load(argv[1]);
    if (!mesh) {
        fprintf(stderr, "failed to load ccm");
        return 1;
    }

    return 0;
}