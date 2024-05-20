#include <mpi.h>
#include <functional>
#include <cstdlib>
#include <unistd.h>
#include <cstring>
#include <unordered_map>
#include <cstdint>
#include <cstdio>

#define NAMELEN 32
#define NUM_OF_PRIMS 6 // Simulate 6 primitives
#define NUM_BUCKETS 8
#define KEY_START_VALUE 100

int local_cid = 0;
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

void init_comm_prof(MPI_Comm comm, char c) {
    prof_metadata *metadata;
    comm_profiler *prof;
    metadata = new prof_metadata();
    prof = new comm_profiler();
    PMPI_Comm_size(comm, &metadata->size);
    metadata->comms = local_cid++;
    metadata->id = c;
    PMPI_Comm_set_attr(comm, keyval[0], metadata);
    PMPI_Comm_set_attr(comm, keyval[1], prof);
}

void insertOrUpdatePrimBucketInfo(std::unordered_map<int, primBucketInfo>& map,
                                  int key, double time, uint64_t volume) {

    // Create a new primBucketInfo object
    primBucketInfo newInfo;
    newInfo.time = time;
    newInfo.num_messages = 1;
    newInfo.volume = volume;

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


// Define a function to free the prof_metadata object
void delete_comm_prof_metadata(MPI_Comm comm, int keyval, void* attr_val, void* extra_state) {
  delete static_cast<prof_metadata*>(attr_val);
}

void delete_comm_profiler(MPI_Comm comm, int keyval, void* attr_val, void* extra_state) {
  delete static_cast<comm_profiler*>(attr_val);
  std::unordered_map<int, primBucketInfo*> *map = (std::unordered_map<int, primBucketInfo*> *)attr_val;
  for (auto it = map->begin(); it != map->end(); ++it) {
    delete it->second;
  }
  delete map;
}

int main(int argc, char *argv[])
{
    int rank, size;
    int i;
    MPI_Init(&argc,&argv);

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);


    for (i=0; i<2; i++) {
        // create key for the object
        MPI_Comm_create_keyval(MPI_COMM_DUP_FN,MPI_COMM_NULL_DELETE_FN,
                               &keyval[i],NULL);
    }
    for (i=0; i<2; i++) {
        printf("keyval[%d] = %d\n",i,keyval[i]);
    }
    init_comm_prof(MPI_COMM_WORLD, 'W');

    prof_metadata *met;
    int flag;
    PMPI_Comm_get_attr(MPI_COMM_WORLD,keyval[0],&met,&flag);
    if (flag){
        printf("Rank %d: nsize = %d, comms = %d, id = %c\n",
               rank,met->size,met->comms,met->id);
    }
    else {
        printf("Attribute not found\n");
    }
    comm_profiler *comm_prof;
    PMPI_Comm_get_attr(MPI_COMM_WORLD, keyval[1], &comm_prof, &flag);
    if (flag) {
        printf("Rank %d: Found the map\n", rank);
    } else {
        printf("Rank %d: Map not found\n", rank);
    }
    printf("Rank %d: Inserting data into map\n", rank);
    // insert data into map
    insertOrUpdatePrimBucketInfo(comm_prof->map, getPrimBucketKey(1, 2), 1.0,  1000);
    //comm_prof->map[getPrimBucketKey(1, 2)] = profData;
    // printf("Rank %d: Inserted data into map\n", rank);
    // // get inserted data
    insertOrUpdatePrimBucketInfo(comm_prof->map, getPrimBucketKey(1, 2), 2.0,  1000);
    insertOrUpdatePrimBucketInfo(comm_prof->map, getPrimBucketKey(0, 0), 5.0, 100);

    // primBucketInfo profData2;
    // profData2.time = 2.0;
    // profData2.num_messages = 20;
    // profData2.volume = 2000;

    // insertOrUpdatePrimBucketInfo(comm_prof->map, getPrimBucketKey(1, 2), profData2);
    primBucketInfo data;
    data = comm_prof->map[getPrimBucketKey(1, 2)];
    printf("Rank %d: Updated data in map: time = %f, num_messages = %d, volume = %lu\n", rank, data.time, data.num_messages, data.volume);
    data = comm_prof->map[getPrimBucketKey(0, 0)];
    printf("Rank %d: Updated data in map: time = %f, num_messages = %d, volume = %lu\n", rank, data.time, data.num_messages, data.volume);

    MPI_Comm newcomm;
    MPI_Comm_split(MPI_COMM_WORLD, rank, rank, &newcomm);
    init_comm_prof(newcomm, 'S');
    PMPI_Comm_get_attr(newcomm,keyval[0],&met,&flag);
    if (flag){
        printf("Rank %d: nsize = %d, comms = %d, id = %c\n",
               rank,met->size,met->comms,met->id);
    }
    else {
        printf("Attribute not found\n");
    }
    PMPI_Comm_get_attr(newcomm, keyval[1], &comm_prof, &flag);
    if (flag) {
        printf("Rank %d: Found the map\n", rank);
    } else {
        printf("Rank %d: Map not found\n", rank);
    }
    insertOrUpdatePrimBucketInfo(comm_prof->map, getPrimBucketKey(1, 2), 0.5,  10);
    insertOrUpdatePrimBucketInfo(comm_prof->map, getPrimBucketKey(1, 2), 5.0, 100);
    data = comm_prof->map[getPrimBucketKey(1, 2)];
    printf("Rank %d: Updated data in map: time = %f, num_messages = %d, volume = %lu\n", rank, data.time, data.num_messages, data.volume);

    // prof_metadata metadata2;
    // strncpy(metadata2.name, "SPLIT", NAMELEN);
    // metadata2.size = size;
    // metadata2.comms = 1;
    // metadata2.id = 'T';

    // prof_metadata *met2;
    // PMPI_Comm_set_attr(newcomm,keyval[0],&metadata2);
    // PMPI_Comm_get_attr(newcomm,keyval[0],&met2,&flag);
    // if (flag){
    //     printf("Rank %d: name = %s, size = %d, comms = %d, id = %c\n",
    //            rank,met2->name,met2->size,met2->comms,met2->id);
    // }
    // else {
    //     printf("Attribute not found\n");
    // }
    // PMPI_Comm_set_attr(newcomm, keyval[1], &prof);

    // PMPI_Comm_get_attr(newcomm, keyval[1], &comm_prof, &flag);
    // comm_prof->map[getPrimBucketKey(2, 3)] = profData;
    // data = comm_prof->map[getPrimBucketKey(2, 3)];
    // printf("Rank %d: Retrieved data from map: time = %f, num_messages = %d, volume = %lu\n", rank, data.time, data.num_messages, data.volume);



    MPI_Finalize();
    return 0;
}
