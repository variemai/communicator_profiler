#include <iostream>
#include <sqlite3.h>

static int callback(void *NotUsed, int argc, char **argv, char **azColName) {
   int i;
   for(i = 0; i<argc; i++) {
      printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
   }
   printf("\n");
   return 0;
}

// Function to execute SQL command
void executeSQL(sqlite3* db, const char* sql, const char* name) {
    char* errMsg = nullptr;
    int rc = sqlite3_exec(db, sql, nullptr, nullptr, &errMsg);

    if (rc != SQLITE_OK) {
        std::cerr << "SQL error: " << errMsg << std::endl;
        sqlite3_free(errMsg);
    } else {
        std::cout << name << " successfully" << std::endl;
    }
}

// Function to execute SQL command
void executeSQL(sqlite3* db, const std::string& sql, const char* name) {
    char* errMsg = nullptr;
    int rc = sqlite3_exec(db, sql.c_str(), nullptr, nullptr, &errMsg);

    if (rc != SQLITE_OK) {
        std::cerr << "SQL error: " << errMsg << std::endl;
        sqlite3_free(errMsg);
    } else {
        std::cout << name << " successfully" << std::endl;
    }
}



void createTables(sqlite3* db) {

    // Create Mapping Table
    const char* MappingTable =
        "CREATE TABLE IF NOT EXISTS mappings ("
        "id INTEGER PRIMARY KEY AUTOINCREMENT, "
        "machine TEXT);";
    executeSQL(db, MappingTable, "Mappings Table created");

    // Create MPI Operations Table
    const char* MPIOpsTable =
        "CREATE TABLE IF NOT EXISTS operations ("
        "id INTEGER PRIMARY KEY AUTOINCREMENT, "
        "operation TEXT);";
    executeSQL(db, MPIOpsTable, "MPI Operations Table");

    // Create Comms Table
    const char* CommsTable =
        "CREATE TABLE IF NOT EXISTS comms ("
        "id INTEGER PRIMARY KEY AUTOINCREMENT, "
        "name TEXT, "
        "size TEXT);";
    executeSQL(db, CommsTable, "Communicator Table created");

    // Create Data Table
    const char* DataTable =
        "CREATE TABLE IF NOT EXISTS data ("
        "id INTEGER PRIMARY KEY AUTOINCREMENT, "
        "rank INTEGER, "
        "mapping_id INTEGER, "
        "comm_id INTEGER, "
        "operation_id INTEGER, "
        "buffer_size_max INTEGER, "
        "buffer_size_min INTEGER, "
        "calls INTEGER, "
        "time REAL, "
        "FOREIGN KEY (operation_id) REFERENCES operations (id), "
        "FOREIGN KEY (comm_id) REFERENCES comms (id), "
        "FOREIGN KEY (mapping_id) REFERENCES mappings (id));";
    executeSQL(db, DataTable, "Data");
}

// Function to insert into mappings
void insertIntoMappings(sqlite3* db, const std::string& machine) {
    std::string sql = "INSERT INTO mappings (machine) VALUES ('" + machine + "')";
    executeSQL(db, sql, "INSERT INTO mappings");
}

// Function to insert into comms
void insertIntoComms(sqlite3* db, const std::string& name, const std::string& size) {
    std::string sql = "INSERT INTO comms (name, size) VALUES ('" + name + "', '" + size + "')";
    executeSQL(db, sql,"INSERT INTO comms" );
}

// Function to insert into operations
void insertIntoOperations(sqlite3* db, const std::string& operation) {
    std::string sql = "INSERT INTO operations (operation) VALUES ('" + operation + "')";
    executeSQL(db, sql, "INSERT INTO operations");
}

// Function to insert into data
void insertIntoData(sqlite3* db, int rank, int mappingId, int commId, int operationId, int bufferSizeMax, int bufferSizeMin, int calls, double time) {
    std::string sql = "INSERT INTO data (rank, mapping_id, comm_id, operation_id, buffer_size_max, buffer_size_min, calls, time) VALUES ("
                      + std::to_string(rank) + ", " + std::to_string(mappingId) + ", " + std::to_string(commId) + ", "
                      + std::to_string(operationId) + ", " + std::to_string(bufferSizeMax) + ", "
                      + std::to_string(bufferSizeMin) + ", " + std::to_string(calls) + ", " + std::to_string(time) + ")";
    executeSQL(db, sql, "INSERT INTO data");
}



int main(int argc, char* argv[]) {
  sqlite3 *db;
  char *zErrMsg = 0;
  int rc;

  // Open database
  rc = sqlite3_open("test_mpi_data.db", &db);
  if (rc) {
      std::cerr << "Can't open database: " << sqlite3_errmsg(db) << std::endl;
      return 1;
  } else {
      std::cout << "Opened database successfully" << std::endl;
  }

  createTables(db);

      insertIntoMappings(db, "machine1");
    insertIntoMappings(db, "machine2");

    // Insert test data into comms
    insertIntoComms(db, "comm1", "1024");
    insertIntoComms(db, "comm2", "2048");

    // Insert test data into operations
    insertIntoOperations(db, "MPI_Send");
    insertIntoOperations(db, "MPI_Recv");
    insertIntoOperations(db, "MPI_Allreduce");
    insertIntoOperations(db, "MPI_Bcast");

    // Insert test data into data for rank 0
    insertIntoData(db, 0, 1, 1, 1, 128, 64, 10, 0.001); // MPI_Send
    insertIntoData(db, 0, 2, 2, 2, 256, 128, 5, 0.002); // MPI_Recv
    insertIntoData(db, 0, 1, 1, 3, 512, 256, 8, 0.003); // MPI_Allreduce
    insertIntoData(db, 0, 2, 2, 4, 1024, 512, 3, 0.004); // MPI_Bcast

    // Insert test data into data for rank 1
    insertIntoData(db, 1, 1, 1, 1, 128, 64, 15, 0.005); // MPI_Send
    insertIntoData(db, 1, 2, 2, 2, 256, 128, 7, 0.006); // MPI_Recv
    insertIntoData(db, 1, 1, 1, 3, 512, 256, 10, 0.007); // MPI_Allreduce
    insertIntoData(db, 1, 2, 2, 4, 1024, 512, 5, 0.008); // MPI_Bcast

    sqlite3_close(db);



  sqlite3_close(db);
  return 0;
}

// int main(int argc, char* argv[]) {
//    sqlite3 *db;
//    char *zErrMsg = 0;
//    int rc;
//    char *sql;


//    if( rc ) {
//       fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
//       return(0);
//    } else {
//       fprintf(stderr, "Opened database successfully\n");
//    }


//    /* Execute SQL statement */
//    rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);

//    if( rc != SQLITE_OK ){
//       fprintf(stderr, "SQL error: %s\n", zErrMsg);
//       sqlite3_free(zErrMsg);
//    } else {
//       fprintf(stdout, "Records created successfully\n");
//    }
//    sqlite3_close(db);
//    return 0;
// }
