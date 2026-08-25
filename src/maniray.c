#include <lis.h>

#include "maniray/maniray.h"
#include "maniray/utils/misc.h"

int mr_initialize(int *argc, char ***argv) {
    return lis_initialize(argc, argv) == LIS_SUCCESS ? MR_SUCCESS : MR_FAILURE;
}

int mr_finalize() {
    return lis_finalize() == LIS_SUCCESS ? MR_SUCCESS : MR_FAILURE;
}