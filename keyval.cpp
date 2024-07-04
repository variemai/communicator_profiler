#include <mpi.h>
#include <cstdlib>
#include <unistd.h>
#include <cstring>
#include <unordered_map>
#include <cstdint>
#include <cstdio>
#include <utility>
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
    int comm_id;
    int prim;
    int bucketIndex;
    double time;
    int num_messages;
    uint64_t volume;
} comm_data;

typedef struct comm_all{
    char name[NAMELEN];
    int size;
} comm_all;


std::vector<profiler_metadata> metadata_freelist;
std::vector<comm_profiler> profiler_freelist;
std::vector<MPI_Comm> comms;

typedef struct comm_profiler_meta_pair{
    comm_profiler prof;
    prof_metadata meta;
} prof_meta_pair;

std::vector<std::pair<prof_meta_pair*, MPI_Group>> free_array;


// Helper function to create unique key
// int getPrimBucketKey(int prim, int bucketIndex) {
//     std::hash<int> hash_int;
//     size_t combinedHash = hash_int(prim) ^ (hash_int(bucketIndex) << 1);
//     return combinedHash;
// }

int getPrimBucketKey(int prim, int bucketIndex)
{
    //std::cout << "prim = " << prim << ", bucketIndex = " << bucketIndex << " key = " << prim * NUM_BUCKETS + bucketIndex << std::endl;
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

// void printProfilingData(const std::vector<comm_all>& array) {
//     for (const auto& comm : array) {  // Iterate over each comm_all element
//         std::cout << "Communicator Name: " << comm.name << "\n";
//         std::cout << "  Size: " << comm.size << "\n";

//         std::cout << "\n  Profiling Data:\n";
//         for (const auto& entry : comm.data) {  // Iterate over each comm_data entry
//             std::cout << "    Primitive: " << prim_names[entry.prim] << std::endl;
//             std::cout << "    Bucket Index: " << entry.bucketIndex << std::endl;
//             std::cout << "    Time: " << entry.time << std::endl;
//             std::cout << "    Num Messages: " << entry.num_messages << std::endl;
//             std::cout << "    Volume: " << entry.volume << std::endl;
//         }

//         std::cout << "--------------------\n";
//     }
// }


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
    comms.push_back(MPI_COMM_WORLD);

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
    comms.push_back(newcomm);
    prof_metadata *met2;
    PMPI_Comm_get_attr(newcomm,keyval[0],&met2,&flag);
    if (flag){
        printf("Rank %d: nsize = %d, comms = %d, id = %c\n",
               rank,met2->size,met2->comms,met2->id);
        strcpy(met2->name,"SPLIT");
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
    if ( rank !=0 ){
        insertOrUpdatePrimBucketInfo(split_comm->map, getPrimBucketKey(3, 1), 2.0,  50);
        data = split_comm->map[getPrimBucketKey(3, 1)];
        printf("Rank %d: Updated data in map: time = %f, num_messages = %d, volume = %lu\n", rank, data.time, data.num_messages, data.volume);
    }
    //metadata_freelist.push_back(*met2);
    //profiler_freelist.push_back(*split_comm);
    MPI_Group group;
    MPI_Comm_group(newcomm, &group);
    prof_meta_pair *free_pair = new prof_meta_pair();
    free_pair->meta = *met2;
    free_pair->prof = *split_comm;
    free_array.push_back(std::make_pair(free_pair, group));

    // Find newcomm in comms and remove it
    for (i=0; i<comms.size(); i++){
        if (comms[i] == newcomm){
            comms.erase(comms.begin() + i);
            break;
        }
    }
    comms.shrink_to_fit();
    MPI_Comm_free(&newcomm);
    std::cout << "comm size = " << comms.size() << std::endl;
    // Place the profiling data into a vector and gather it to rank 0
    std::vector<comm_all> array;
    std::vector<comm_data> data_array;
    comm_all comm_meta;
    comm_profiler *world_prof;
    int commid;
    for (i=0; i<free_array.size(); i++){
        PMPI_Comm_create_group(MPI_COMM_WORLD, free_array[i].second, 0, &newcomm);
        PMPI_Comm_set_attr(newcomm, keyval[0], &free_array[i].first->meta);
        PMPI_Comm_set_attr(newcomm, keyval[1], &free_array[i].first->prof);
        comms.push_back(newcomm);
    }

    for (commid =0; commid<comms.size(); commid++){
        PMPI_Comm_get_attr(comms[commid], keyval[1], &world_prof, &flag);
        if ( !flag ){
            printf("Map not found\n");
        }
        for (i=0; i<NUM_OF_PRIMS; ++i) {
            for ( j =0; j<NUM_BUCKETS; ++j ){
                // Check if the key exists in the map
                auto it = world_prof->map.find(getPrimBucketKey(i, j));

                if (it != world_prof->map.end()) {
                    // Allocate comm_data struct and copy the data from map
                    comm_data data;
                    data.comm_id = commid;
                    data.prim = i;
                    data.bucketIndex = j;
                    data.num_messages = it->second.num_messages;
                    data.time = it->second.time;
                    data.volume = it->second.volume;
                    data_array.push_back(data);
                    // printf("Rank %d: Primitive = %d, Bucket = %d, Time = %f, Num Messages = %d, Volume = %lu\n", rank, i, j, it->second.time, it->second.num_messages, it->second.volume);
                }
            }
        }
        PMPI_Comm_get_attr(comms[commid],keyval[0],&met,&flag);
        if (flag) {
            strcpy(comm_meta.name, met->name);
            comm_meta.size = met->size;
            array.push_back(comm_meta);
        }
        else{
            printf("Rank %d: Metadata not found\n", rank);
        }
    }

    // Gather the profiling metadata from all ranks to rank 0
    int local_size = array.size();
    int *c_recvcounts = (int*)malloc(size * sizeof(int));
    int *c_displs = (int*)malloc(size * sizeof(int));
    int total_num_of_comms = 0;
    comm_all dummy;
    MPI_Gather(&local_size, 1, MPI_INT, c_recvcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // Root process: Allocate receiving buffer and calculate displacements
    std::vector<comm_all> recv_buffer;
    if (rank == 0) {
        c_displs[0] = 0;
        total_num_of_comms = c_recvcounts[0];
        for (i = 1; i < size; ++i) {
          c_displs[i] = c_displs[i - 1] + c_recvcounts[i - 1];
          total_num_of_comms += c_recvcounts[i];
        }
        std::cout << "mpisee: total number of communicators = " << total_num_of_comms << std::endl;
        recv_buffer.resize(total_num_of_comms);

    }
    // Gather the data
    MPI_Datatype MPI_COMM_ALL;
    MPI_Datatype types[2] = {MPI_CHAR, MPI_INT};
    int blocklengths[2] = {NAMELEN, 1};
    MPI_Aint displacements[2];
    MPI_Aint base;
    MPI_Get_address(&dummy, &base);
    MPI_Get_address(&dummy.name, &displacements[0]);
    MPI_Get_address(&dummy.size, &displacements[1]);
   // Convert addresses to displacements
    for (i = 0; i < 2; i++) {
        displacements[i] = MPI_Aint_diff(displacements[i], base);
    }

    MPI_Type_create_struct(2, blocklengths, displacements, types, &MPI_COMM_ALL);
    MPI_Type_commit(&MPI_COMM_ALL);

    MPI_Gatherv(array.data(), local_size, MPI_COMM_ALL,
                recv_buffer.data(), c_recvcounts, c_displs, MPI_COMM_ALL, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        for (i = 0; i < total_num_of_comms; ++i) {
            std::cout << "Communicator Name: " << recv_buffer[i].name << "\n";
            std::cout << "  Size: " << recv_buffer[i].size << "\n";
        }
    }

    // Free the MPI datatype
    MPI_Type_free(&MPI_COMM_ALL);
    // Free the "array" vector
    array.clear();
    array.shrink_to_fit(); // recv_buffer has all data now


// 1. Create MPI datatype for comm_data
    comm_data dummy_data;
    MPI_Datatype MPI_COMM_DATA;
    MPI_Datatype datatypes[6] = {MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_INT, MPI_UINT64_T};
    int blocklengths2[6] = {1, 1, 1, 1, 1, 1};
    MPI_Aint displacements2[6]; // Use more descriptive names for clarity
    MPI_Aint comm_data_base;    // More descriptive name

    MPI_Get_address(&dummy_data, &comm_data_base);

    MPI_Get_address(&dummy_data.comm_id, &displacements2[0]);
    MPI_Get_address(&dummy_data.prim, &displacements2[1]);
    MPI_Get_address(&dummy_data.bucketIndex, &displacements2[2]);
    MPI_Get_address(&dummy_data.time, &displacements2[3]);
    MPI_Get_address(&dummy_data.num_messages, &displacements2[4]);
    MPI_Get_address(&dummy_data.volume, &displacements2[5]);

    for (i = 0; i < 6; ++i) {
        displacements2[i] = MPI_Aint_diff(displacements2[i], comm_data_base);
    }

    MPI_Type_create_struct(6, blocklengths2, displacements2, datatypes, &MPI_COMM_DATA);

    MPI_Type_commit(&MPI_COMM_DATA);

    // 2. Gather the profiling data from all ranks to rank 0
    int local_data_size = data_array.size();
    int total_num_of_data = 0;
    int *recvcounts = (int*)malloc(size * sizeof(int));
    int *displs = (int*)malloc(size * sizeof(int));
    MPI_Gather(&local_data_size, 1, MPI_INT, recvcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // Root process: Allocate receiving buffer and calculate displacements
    std::vector<comm_data> recv_data_buffer;
    if (rank == 0) {
        displs[0] = 0;
        total_num_of_data = recvcounts[0];
        for (i = 1; i < size; ++i) {
          displs[i] = displs[i - 1] + recvcounts[i - 1];
          total_num_of_data += recvcounts[i];
        }
        std::cout << "mpisee: total number of data = " << total_num_of_data << std::endl;
        recv_data_buffer.resize(total_num_of_data);
    }
    // Gather the data
    MPI_Gatherv(data_array.data(), local_data_size, MPI_COMM_DATA,
                recv_data_buffer.data(), recvcounts, displs, MPI_COMM_DATA, 0, MPI_COMM_WORLD);

    // Clear the data_array vector
    data_array.clear();
    data_array.shrink_to_fit(); // recv_data_buffer has all data now

    int procs;
    int comms_per_proc;
    int data_per_procs;
    if (rank == 0) {
        for ( procs = 0; procs < size; ++procs ){
            comms_per_proc = c_recvcounts[procs];
            data_per_procs = recvcounts[procs];
            for (i = 0; i < comms_per_proc; ++i) {
                std::cout << "Communicator Name: " << recv_buffer[c_displs[procs] + i].name << "\n";
                std::cout << "  Size: " << recv_buffer[c_displs[procs] + i].size << "\n";
                for (j = 0; j < data_per_procs; ++j) {
                    if ( recv_data_buffer[displs[procs] + j].comm_id != i ) {
                        continue;
                    }
                    std::cout << "    Primitive: " << prim_names[recv_data_buffer[displs[procs] + j].prim] << "\n";
                    std::cout << "    Bucket Index: " << recv_data_buffer[displs[procs] + j].bucketIndex << "\n";
                    std::cout << "    Time: " << recv_data_buffer[displs[procs] + j].time << "\n";
                    std::cout << "    Num Messages: " << recv_data_buffer[displs[procs] + j].num_messages << "\n";
                    std::cout << "    Volume: " << recv_data_buffer[displs[procs] + j].volume << "\n";
                }
            }
        }
    }

    MPI_Type_free(&MPI_COMM_DATA);
    // Free the arrays you allocated
    free(recvcounts);
    free(displs);
    free(c_recvcounts);
    free(c_displs);
    MPI_Finalize();
    return 0;
}
