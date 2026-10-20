#ifndef EDLIB_H
#define EDLIB_H

#ifdef __cplusplus
extern "C" {
#endif

#define EDLIB_STATUS_OK 0
#define EDLIB_STATUS_ERROR 1

typedef enum {
    EDLIB_MODE_NW,
    EDLIB_MODE_SHW,
    EDLIB_MODE_HW
} EdlibAlignMode;

typedef enum {
    EDLIB_TASK_DISTANCE,
    EDLIB_TASK_LOC,
    EDLIB_TASK_PATH
} EdlibAlignTask;

typedef enum {
    EDLIB_CIGAR_STANDARD,
    EDLIB_CIGAR_EXTENDED
} EdlibCigarFormat;

#define EDLIB_EDOP_MATCH 0
#define EDLIB_EDOP_INSERT 1
#define EDLIB_EDOP_DELETE 2
#define EDLIB_EDOP_MISMATCH 3

typedef struct {
    char first;
    char second;
} EdlibEqualityPair;

typedef struct {
    int k;
    EdlibAlignMode mode;
    EdlibAlignTask task;
    EdlibEqualityPair* additionalEqualities;
    int additionalEqualitiesLength;
} EdlibAlignConfig;

typedef struct {
    int status;
    int editDistance;
    int* endLocations;
    int* startLocations;
    int numLocations;
    unsigned char* alignment;
    int alignmentLength;
    int alphabetLength;
} EdlibAlignResult;

EdlibAlignConfig edlibNewAlignConfig(
    int k,
    EdlibAlignMode mode,
    EdlibAlignTask task,
    EdlibEqualityPair* additionalEqualities,
    int additionalEqualitiesLength
);

EdlibAlignConfig edlibDefaultAlignConfig(void);

void edlibFreeAlignResult(EdlibAlignResult result);

EdlibAlignResult edlibAlign(
    const char* query,
    int queryLength,
    const char* target,
    int targetLength,
    const EdlibAlignConfig config
);

char* edlibAlignmentToCigar(
    const unsigned char* alignment,
    int alignmentLength,
    EdlibCigarFormat cigarFormat
);

#ifdef __cplusplus
}
#endif

#endif
