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

    // Create Metrics Table
    const char* createMetricsTable =
        "CREATE TABLE IF NOT EXISTS metrics ("
        "id INTEGER PRIMARY KEY AUTOINCREMENT, "
        "metric_name TEXT, "
        "metric_unit TEXT);";
    executeSQL(db, createMetricsTable, "Metrics");

    // Create Mapping Table
    const char* createMappingTable =
        "CREATE TABLE IF NOT EXISTS mapping ("
        "id INTEGER PRIMARY KEY AUTOINCREMENT, "
        "mapping_description TEXT);";
    executeSQL(db, createMappingTable, "Mapping");

    // Create Data Table
    const char* createDataTable =
        "CREATE TABLE IF NOT EXISTS data ("
        "id INTEGER PRIMARY KEY AUTOINCREMENT, "
        "rank INTEGER, "
        "mapping_id INTEGER, "
        "metric_id INTEGER, "
        "value REAL, "
        "FOREIGN KEY (mapping_id) REFERENCES mapping (id), "
        "FOREIGN KEY (metric_id) REFERENCES metrics (id));";
    executeSQL(db, createDataTable, "Data");

    // Close database
    sqlite3_close(db);
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
