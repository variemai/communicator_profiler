#include <mpi.h>
#include <functional>
#include <cstdlib>
#include <unistd.h>
#include <cstring>
#include <unordered_map>

#define NAMELEN 32
#define NUM_OF_PRIMS 6 // Simulate 6 primitives
#define NUM_BUCKETS 8
#define KEY_START_VALUE 100

int keyval[2]; // keyval[0]  contains metadata
                   // keyval[1]  contains profiling data

const int64_t buckets[NUM_BUCKETS-1] = {128,1024,8192,65536,262144,1048576,33554432};
// Metadata associated with an MPI communicator
typedef struct profiler_metadata {
    char name[NAMELEN];
    int size;
    int comms;
    char id;
} prof_metadata;

// Structure to store profiling data for a primitive and bucket
typedef struct PrimBucketInfo {
    double time;
    int num_messages;
    uint64_t volume;
}primBucketInfo;

typedef struct comm_profiler {
    std::unordered_map<int, primBucketInfo> map;
} comm_profiler;

// Helper function to create unique key
int getPrimBucketKey(int prim, int bucketIndex) {
    std::hash<int> hash_int;
    size_t combinedHash = hash_int(prim) ^ (hash_int(bucketIndex) << 1);
    int keyRangeSize = NUM_OF_PRIMS * NUM_BUCKETS;
    return combinedHash % keyRangeSize;
}

void init_comm_prof(MPI_Comm comm, prof_metadata *metadata) {
    comm_profiler prof;
    PMPI_Comm_set_attr(comm, keyval, &prof);
}

void insertOrUpdatePrimBucketInfo(std::unordered_map<int, primBucketInfo>& map,
                                  int key,
                                  const primBucketInfo& newInfo) {

    // Check if the key exists in the map
    auto it = map.find(key);

    if (it == map.end()) {
        // Key not found, insert the new pair
        map[key] = newInfo;
    } else {
        // Key found, update the existing value
        it->second.time += newInfo.time;
        it->second.num_messages += newInfo.num_messages;
        it->second.volume += newInfo.volume;
    }
}

int main(int argc, char *argv[])
{
    int rank, size;
    int i;
    MPI_Init(&argc,&argv);

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    prof_metadata metadata;
    strncpy(metadata.name, "WORLD", NAMELEN);
    metadata.size = size;
    metadata.comms = 1;
    metadata.id = 'W';

    for (i=0; i<2; i++) {
        // create key for the object
        MPI_Comm_create_keyval(MPI_COMM_DUP_FN,MPI_COMM_NULL_DELETE_FN,
                               &keyval[i],NULL);
    }
    for (i=0; i<2; i++) {
        printf("keyval[%d] = %d\n",i,keyval[i]);
    }
    PMPI_Comm_set_attr(MPI_COMM_WORLD,keyval[0],&metadata);

    prof_metadata *met;
    int flag;
    PMPI_Comm_get_attr(MPI_COMM_WORLD,keyval[0],&met,&flag);
    if (flag){
        printf("Rank %d: name = %s, size = %d, comms = %d, id = %c\n",
               rank,met->name,met->size,met->comms,met->id);
    }
    else {
        printf("Attribute not found\n");
    }
    // std::unordered_map<int, primBucketInfo> primBucketMap;
    comm_profiler prof;
    PMPI_Comm_set_attr(MPI_COMM_WORLD, keyval[1], &prof);

    primBucketInfo profData;
    profData.time = 1.0;
    profData.num_messages = 10;
    profData.volume = 1000;

    comm_profiler *comm_prof;
    PMPI_Comm_get_attr(MPI_COMM_WORLD, keyval[1], &comm_prof, &flag);
    insertOrUpdatePrimBucketInfo(comm_prof->map, getPrimBucketKey(0, 1), profData);
    // PMPI_Comm_get_attr(MPI_COMM_WORLD, keyval[1], &map, &flag);
    if (flag) {
        printf("Rank %d: Found the map\n", rank);
    } else {
        printf("Rank %d: Map not found\n", rank);
    }
    comm_prof->map[getPrimBucketKey(1, 2)] = profData;
    printf("Rank %d: Inserted data into map\n", rank);
    // get inserted data
    primBucketInfo data = comm_prof->map[getPrimBucketKey(1, 2)];
    printf("Rank %d: Retrieved data from map: time = %f, num_messages = %d, volume = %llu\n", rank, data.time, data.num_messages, data.volume);

    primBucketInfo profData2;
    profData2.time = 2.0;
    profData2.num_messages = 20;
    profData2.volume = 2000;

    insertOrUpdatePrimBucketInfo(comm_prof->map, getPrimBucketKey(1, 2), profData2);
    data = comm_prof->map[getPrimBucketKey(1, 2)];
    printf("Rank %d: Updated data in map: time = %f, num_messages = %d, volume = %llu\n", rank, data.time, data.num_messages, data.volume);

    MPI_Comm newcomm;
    MPI_Comm_split(MPI_COMM_WORLD, rank, rank, &newcomm);

    prof_metadata metadata2;
    strncpy(metadata2.name, "SPLIT", NAMELEN);
    metadata2.size = size;
    metadata2.comms = 1;
    metadata2.id = 'T';

    prof_metadata *met2;
    PMPI_Comm_set_attr(newcomm,keyval[0],&metadata2);
    PMPI_Comm_get_attr(newcomm,keyval[0],&met2,&flag);
    if (flag){
        printf("Rank %d: name = %s, size = %d, comms = %d, id = %c\n",
               rank,met2->name,met2->size,met2->comms,met2->id);
    }
    else {
        printf("Attribute not found\n");
    }
    PMPI_Comm_set_attr(newcomm, keyval[1], &prof);

    PMPI_Comm_get_attr(newcomm, keyval[1], &comm_prof, &flag);
    comm_prof->map[getPrimBucketKey(2, 3)] = profData;
    data = comm_prof->map[getPrimBucketKey(2, 3)];
    printf("Rank %d: Retrieved data from map: time = %f, num_messages = %d, volume = %llu\n", rank, data.time, data.num_messages, data.volume);



    MPI_Finalize();
    return 0;
}
