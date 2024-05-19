#include <mpi.h>
#include <iostream>
#include <functional>
#include <cstdlib>
#include <unistd.h>

#define NAMELEN 32
#define NUM_OF_PRIMS 6 // Simulate 6 primitives
#define NUM_BUCKETS 8
#define KEY_START_VALUE 100

const char prim_names[][NUM_OF_PRIMS] = {"Send", "Recv", "Allreduce", "Bcast", "Alltoall", "Reduce"};

// Metadata associated with an MPI communicator
typedef struct profiler_metadata {
    char name[NAMELEN];
    int size;
    int comms;
    char id;
} prof_metadata;

// Structure to store profiling data for a primitive and bucket
struct PrimBucketInfo {
    double time;
    int num_messages;
    uint64_t volume;
};

std::unordered_map<int, PrimBucketInfo*> bucket_table;

extern "C" {
    int namekey(void) {
        static int baseKeyval = MPI_KEYVAL_INVALID;

        if (baseKeyval == MPI_KEYVAL_INVALID) {
            // Reserve keys for metadata and PrimBucketInfo
            for (int i = 0; i < NUM_OF_PRIMS * NUM_BUCKETS + 1; ++i) { // +1 for metadata
                int tempKey;
                PMPI_Comm_create_keyval(MPI_COMM_NULL_COPY_FN, MPI_COMM_NULL_DELETE_FN, &tempKey, NULL);
            }
            baseKeyval = KEY_START_VALUE; // Set the base key value
        }
        return baseKeyval;
    }
}
// Helper function to create unique key
int getPrimBucketKey(int prim, int bucketIndex) {
    std::hash<int> hash_int;

    size_t combinedHash = hash_int(prim) ^ (hash_int(bucketIndex) << 1);
    int keyRangeSize = NUM_OF_PRIMS * NUM_BUCKETS;
    int offset = combinedHash % keyRangeSize;

    return namekey() + offset;
}

extern "C" {
void profile_this(MPI_Comm comm, int64_t count, MPI_Datatype datatype, int prim, double t_elapsed, int v) {
        int bucketIndex = choose_bucket(sum);
        int key = getPrimBucketKey(prim, bucketIndex);

        PrimBucketInfo* info = bucket_table[key];
        if (!info) {
            info = new PrimBucketInfo();
            bucket_table[key] = info;
        }

        // Update info values based on whether it's a collective operation
        if (v == 0) {
            info->time += t_elapsed;
            info->num_messages += 1;
            info->volume += sum;
        } else {
            info->time += t_elapsed;
            info->num_messages += 1;
            info->volume += sum;
        }

        // Attach the PrimBucketInfo to the communicator only once
        if (!bucket_table_attached[comm]) {
            PMPI_Comm_set_attr(comm, namekey(), &bucket_table); // Attaching the map
            bucket_table_attached[comm] = true;
        }
}

// Helper function to create unique key
int getPrimBucketKey(int prim, int bucketIndex) {
    std::hash<int> hash_int;

    size_t combinedHash = hash_int(prim) ^ (hash_int(bucketIndex) << 1);
    int keyRangeSize = NUM_OF_PRIMS * NUM_BUCKETS;
    int offset = combinedHash % keyRangeSize;

    return namekey() + offset;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    prof_metadata metadata;
    strncpy(metadata.name, "TestComm", NAMELEN);
    metadata.size = size;
    metadata.comms = 1;
    metadata.id = 'T';

    // Attach metadata to the communicator
    PMPI_Comm_set_attr(MPI_COMM_WORLD, namekey(), &metadata);

    srand(rank); // Seed the random number generator

    // Simulate profiling calls without actual communication
    for (int prim = 0; prim < NUM_OF_PRIMS; ++prim) {
        for (int i = 0; i < 10; ++i) { // Make 10 calls for each primitive
            int64_t count = rand() % 1000000 + 1; // Random message size
            double t_elapsed = (double)(rand() % 100) / 1000.0; // Random elapsed time (up to 100ms)
            profile_this(MPI_COMM_WORLD, count, MPI_INT, prim, t_elapsed, 0);
        }
        sleep(1); // For demonstration, not necessary in real usage
    }

    // Retrieve and print metadata and profiling data (for rank 0 only)
    if (rank == 0) {
        prof_metadata* retrieved_metadata;
        int flag;
        PMPI_Comm_get_attr(MPI_COMM_WORLD, namekey(), &retrieved_metadata, &flag);

        if (flag) {
            std::cout << "Retrieved Metadata:\n";
            std::cout << "  Name: " << retrieved_metadata->name << std::endl;
            std::cout << "  Size: " << retrieved_metadata->size << std::endl;

            std::cout << "\nProfiling Data:\n";
            for (int p = 0; p < NUM_OF_PRIMS; ++p) {
                for (int b = 0; b < NUM_BUCKETS; ++b) {
                    int key = getPrimBucketKey(p, b);
                    PrimBucketInfo* info;
                    PMPI_Comm_get_attr(MPI_COMM_WORLD, key, &info, &flag);
                    if (flag) {
                        std::cout << "  " << prim_names[p] << " Bucket " << b << ":\n";
                        std::cout << "    Time: " << info->time << std::endl;
                        std::cout << "    Messages: " << info->num_messages << std::endl;
                        // ... print other PrimBucketInfo fields ...
                    }
                }
            }
        }
    }

    MPI_Finalize();
    return 0;
}
