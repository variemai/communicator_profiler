#include <mpi.h>
#include <functional>
#include <cstdlib>
#include <unistd.h>
#include <cstring>
#include <unordered_map>
#include <cstdint>
#include <cstdio>
#include <vector>
#include <iostream>

#define NAMELEN 32
#define NUM_OF_PRIMS 6 // Simulate 6 primitives
#define NUM_BUCKETS 4
#define KEY_START_VALUE 100

int local_cid = 0;
int keyval[2]; // keyval[0]  contains metadata
               // keyval[1]  contains profiling data
int64_t buckets[NUM_BUCKETS-1] = {128,1024,8192};
const char prim_names[][NAMELEN] = {"Send", "Recv", "Allreduce", "Bcast", "Alltoall", "Reduce"};
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


typedef struct comm_data{
    int prim;
    int bucketIndex;
    double time;
    int num_messages;
    uint64_t volume;
} comm_data;

typedef struct comm_all{
    char name[NAMELEN];
    int size;
    std::vector<comm_data> data;
} comm_all;

// Helper function to create unique key
// int getPrimBucketKey(int prim, int bucketIndex) {
//     std::hash<int> hash_int;
//     size_t combinedHash = hash_int(prim) ^ (hash_int(bucketIndex) << 1);
//     return combinedHash;
// }

int getPrimBucketKey(int prim, int bucketIndex)
{
    std::cout << "prim = " << prim << ", bucketIndex = " << bucketIndex << " key = " << prim * NUM_BUCKETS + bucketIndex << std::endl;
    return prim * NUM_BUCKETS + bucketIndex;
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

    if (map.count(key) > 1) {
        std::cerr << "Hash collision detected for key " << key << std::endl;
    }

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

void printProfilingData(const std::vector<comm_all>& array) {
    for (const auto& comm : array) {  // Iterate over each comm_all element
        std::cout << "Communicator Name: " << comm.name << "\n";
        std::cout << "  Size: " << comm.size << "\n";

        std::cout << "\n  Profiling Data:\n";
        for (const auto& entry : comm.data) {  // Iterate over each comm_data entry
            std::cout << "    Primitive: " << prim_names[entry.prim] << std::endl;
            std::cout << "    Bucket Index: " << entry.bucketIndex << std::endl;
            std::cout << "    Time: " << entry.time << std::endl;
            std::cout << "    Num Messages: " << entry.num_messages << std::endl;
            std::cout << "    Volume: " << entry.volume << std::endl;
        }

        std::cout << "--------------------\n";
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
    int i,j;
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

        strcpy(met->name,"WORLD");
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
    insertOrUpdatePrimBucketInfo(comm_prof->map, getPrimBucketKey(1, 2), 1.0,  500);
    insertOrUpdatePrimBucketInfo(comm_prof->map, getPrimBucketKey(1, 2), 2.0,  1000);
    insertOrUpdatePrimBucketInfo(comm_prof->map, getPrimBucketKey(0, 0), 5.0, 100);

    primBucketInfo data;
    data = comm_prof->map[getPrimBucketKey(1, 2)];
    printf("Rank %d: Updated data in map: time = %f, num_messages = %d, volume = %llu\n", rank, data.time, data.num_messages, data.volume);
    data = comm_prof->map[getPrimBucketKey(0, 0)];
    printf("Rank %d: Updated data in map: time = %f, num_messages = %d, volume = %llu\n", rank, data.time, data.num_messages, data.volume);

    MPI_Comm newcomm;
    MPI_Comm_split(MPI_COMM_WORLD, rank, rank, &newcomm);
    init_comm_prof(newcomm, 'S');
    PMPI_Comm_get_attr(newcomm,keyval[0],&met,&flag);
    if (flag){
        printf("Rank %d: nsize = %d, comms = %d, id = %c\n",
               rank,met->size,met->comms,met->id);
        strcpy(met->name,"SPLIT");
    }
    else {
        printf("Attribute not found\n");
    }
    comm_profiler *split_comm;
    PMPI_Comm_get_attr(newcomm, keyval[1], &split_comm, &flag);
    if (flag) {
        printf("Rank %d: Found the map\n", rank);
    } else {
        printf("Rank %d: Map not found\n", rank);
    }
    insertOrUpdatePrimBucketInfo(split_comm->map, getPrimBucketKey(1, 2), 0.5,  10);
    insertOrUpdatePrimBucketInfo(split_comm->map, getPrimBucketKey(1, 2), 5.0, 100);
    data = split_comm->map[getPrimBucketKey(1, 2)];
    printf("Rank %d: Updated data in map: time = %f, num_messages = %d, volume = %lu\n", rank, data.time, data.num_messages, data.volume);

    std::vector<comm_all> array;
    comm_all comm_meta;
    comm_profiler *world_prof;
    PMPI_Comm_get_attr(MPI_COMM_WORLD, keyval[1], &world_prof, &flag);
    if ( !flag ){
        printf("Map not found\n");
    }
    printf("WORLD\n");
    for (i=0; i<NUM_OF_PRIMS; ++i) {
        for ( j =0; j<NUM_BUCKETS; ++j ){
            // Check if the key exists in the map
            auto it = world_prof->map.find(getPrimBucketKey(i, j));

            if (it != world_prof->map.end()) {
                // Allocate comm_data struct and copy the data from map
                comm_data data;
                data.prim = i;
                data.bucketIndex = j;
                data.num_messages = it->second.num_messages;
                data.time = it->second.time;
                data.volume = it->second.volume;
                comm_meta.data.push_back(data);
                printf("Rank %d: Primitive = %d, Bucket = %d, Time = %f, Num Messages = %d, Volume = %lu\n", rank, i, j, it->second.time, it->second.num_messages, it->second.volume);
            }
        }
    }
    PMPI_Comm_get_attr(MPI_COMM_WORLD,keyval[0],&met,&flag);
    if (flag) {
        strcpy(comm_meta.name, met->name);
        comm_meta.size = met->size;
        array.push_back(comm_meta);
    }
    else{
        printf("Rank %d: Metadata not found\n", rank);
    }


    comm_all comm_meta2;
    PMPI_Comm_get_attr(newcomm, keyval[1], &comm_prof, &flag);
    if ( !flag ){
        printf("Map not found\n");
    }
    printf("SPLIT\n");
    for (i=0; i<NUM_OF_PRIMS; ++i) {
        for ( j =0; j<NUM_BUCKETS; ++j ){
            // Check if the key exists in the map
            auto it = comm_prof->map.find(getPrimBucketKey(i, j));

            if (it != comm_prof->map.end()) {
                // Allocate comm_data struct and copy the data from map
                comm_data data;
                data.prim = i;
                data.bucketIndex = j;
                data.num_messages = it->second.num_messages;
                data.time = it->second.time;
                data.volume = it->second.volume;
                comm_meta2.data.push_back(data);
                printf("Rank %d: Primitive = %d, Bucket = %d, Time = %f, Num Messages = %d, Volume = %lu\n", rank, i, j, it->second.time, it->second.num_messages, it->second.volume);
            }
        }
    }
    PMPI_Comm_get_attr(newcomm,keyval[0],&met,&flag);
    if (flag) {
        strcpy(comm_meta2.name, met->name);
        comm_meta2.size = met->size;
        array.push_back(comm_meta2);
    }
    else{
        printf("Rank %d: Metadata not found\n", rank);
    }

    printProfilingData(array);

    MPI_Finalize();
    return 0;
}
