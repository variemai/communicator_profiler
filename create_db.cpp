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
        std::cout << name << " table created successfully" << std::endl;
    }
}


// Function to create database and tables
void createTables(sqlite3* db) {

    // Create Mapping Table
    const char *MappingTable = "CREATE TABLE IF NOT EXISTS mapping ("
                               "id INTEGER PRIMARY KEY AUTOINCREMENT, "
                               "machine TEXT);";

    executeSQL(db, MappingTable, "Mapping");

    // Create Comms Table
    const char *CommsTable = "CREATE TABLE IF NOT EXISTS comms ("
                             "id INTEGER PRIMARY KEY AUTOINCREMENT, "
                             "name TEXT, "
                             "size TEXT);";

    executeSQL(db, CommsTable, "Communicator");

    // Create Data Table
    const char *DataTable = "CREATE TABLE IF NOT EXISTS data ("
                            "id INTEGER PRIMARY KEY AUTOINCREMENT, "
                            "rank INTEGER, "
                            "mapping_id INTEGER, "
                            "comm_id INTEGER, "
        "operation_id, INTEGER, "
        "buffer_size_max INTEGER, "
        "buffer_size_min INTEGER, "
        "calls INTEGER, "
        "time REAL, "
        "FOREIGN KEY (operation_id) REFERENCES operations (id), "
        "FOREIGN KEY (comm_id) REFERENCES comms (id), "
        "FOREIGN KEY (mapping_id) REFERENCES mapping (id));";

    executeSQL(db, DataTable, "Data");

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
