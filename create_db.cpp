#include <iostream>
#include <sqlite3.h>
#include <string>
#include <iomanip>
#include <sstream>
#include "utils.h"
#include "create_db.h"


int countTable(sqlite3 *db, const std::string& tablename) {
    sqlite3_stmt* stmt;
    int count = 0;
    std::string checkSql = "SELECT COUNT(*) FROM " + tablename;
    sqlite3_prepare_v2(db, checkSql.c_str(), -1, &stmt, nullptr);

    if (sqlite3_step(stmt) == SQLITE_ROW) {
        count = sqlite3_column_int(stmt, 0);
    }
    sqlite3_finalize(stmt);
    return count;
}


void printDataDetails(sqlite3* db) {
    sqlite3_stmt* stmt;

    std::string sql = "SELECT d.rank, m.machine, c.name, c.size, o.operation, "
                      "d.buffer_size_min, d.buffer_size_max, d.calls, d.time "
                      "FROM data d "
                      "JOIN mappings m ON d.rank = m.id "
                      "JOIN comms c ON d.comm_id = c.id "
                      "JOIN operations o ON d.operation_id = o.id";

    if (sqlite3_prepare_v2(db, sql.c_str(), -1, &stmt, nullptr) == SQLITE_OK) {
        while (sqlite3_step(stmt) == SQLITE_ROW) {
            int rank = sqlite3_column_int(stmt, 0);
            const char* machine = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 1));
            const char* commName = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 2));
            const char* commSize = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 3));
            const char* operation = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 4));
            int bufferSizeMin = sqlite3_column_int(stmt, 5);
            int bufferSizeMax = sqlite3_column_int(stmt, 6);
            int calls = sqlite3_column_int(stmt, 7);
            double time = sqlite3_column_double(stmt, 8);

            std::cout << "Rank: " << rank
                      << ", Machine: " << (machine ? machine : "NULL")
                      << ", Comm Name: " << (commName ? commName : "NULL")
                      << ", Comm Size: " << (commSize ? commSize : "NULL")
                      << ", MPI Operation: " << (operation ? operation : "NULL")
                      << ", Buffer Size: " << bufferSizeMin << " - " << bufferSizeMax
                      << ", Calls: " << calls
                      << ", Time: " << time << std::endl;
        }
        sqlite3_finalize(stmt);
    } else {
        std::cerr << "Failed to prepare statement: " << sqlite3_errmsg(db) << std::endl;
    }
}


void printData(sqlite3* db) {
    sqlite3_stmt* stmt;

    std::string sql = "SELECT c.name, c.size, d.rank, o.operation, "
                      "d.buffer_size_min, d.buffer_size_max, d.calls, d.time "
                      "FROM data d "
                      "JOIN comms c ON d.comm_id = c.id "
                      "JOIN operations o ON d.operation_id = o.id "
                      "ORDER BY c.name";  // Order by Comm Name

    if (sqlite3_prepare_v2(db, sql.c_str(), -1, &stmt, nullptr) == SQLITE_OK) {
        // Print header
        std::cout << std::left << std::setw(15) << "Comm Name" << std::setw(15) << "Comm Size"
                  << std::setw(10) << "Rank" << std::setw(20) << "Operation"
                  << std::setw(25) << "Buffer Size Range" << std::setw(15) << "Calls"
                  << std::setw(20) << "Time" << std::endl;

        // Print rows
        while (sqlite3_step(stmt) == SQLITE_ROW) {
          std::stringstream bufferSizeStream;
            bufferSizeStream << sqlite3_column_int(stmt, 4) << " - " << sqlite3_column_int(stmt, 5);
            std::string bufferSize = bufferSizeStream.str();

            std::cout << std::left << std::setw(15) << reinterpret_cast<const char*>(sqlite3_column_text(stmt, 0))
                      << std::setw(15) << reinterpret_cast<const char*>(sqlite3_column_text(stmt, 1))
                      << std::setw(10) << sqlite3_column_int(stmt, 2)
                      << std::setw(20) << reinterpret_cast<const char*>(sqlite3_column_text(stmt, 3))
                      << std::setw(25) << bufferSize //sqlite3_column_int(stmt, 4) << " - " << sqlite3_column_int(stmt, 5)
                      << std::setw(15) << sqlite3_column_int(stmt, 6)
                      << std::setw(20) << sqlite3_column_double(stmt, 7) << std::endl;
        }
        sqlite3_finalize(stmt);
    } else {
        std::cerr << "Failed to prepare statement: " << sqlite3_errmsg(db) << std::endl;
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


int getCommId(sqlite3* db, const std::string& commName) {
    sqlite3_stmt* stmt;
    int commId = -1;  // Default to an invalid ID
    // Prepare SQL query to select the ID from comms where name and size match
    std::string sql = "SELECT id FROM comms WHERE name = '" + commName + "'";
    if (sqlite3_prepare_v2(db, sql.c_str(), -1, &stmt, nullptr) == SQLITE_OK) {
        sqlite3_bind_text(stmt, 1, commName.c_str(), -1, SQLITE_STATIC);

        if (sqlite3_step(stmt) == SQLITE_ROW) {
            commId = sqlite3_column_int(stmt, 0);  // Get the id from the query result
        }
        sqlite3_finalize(stmt);
    } else {
        std::cerr << "Failed to prepare statement: " << sqlite3_errmsg(db) << std::endl;
    }
    return commId; // return the retrieved ID
}

int getMappingId(sqlite3* db, const std::string& machineName) {
    sqlite3_stmt* stmt;
    int mappingId = -1;  // Default to an invalid ID

    std::string sql = "SELECT id FROM mappings WHERE machine = '" + machineName + "'";
    if (sqlite3_prepare_v2(db, sql.c_str(), -1, &stmt, nullptr) == SQLITE_OK) {
        sqlite3_bind_text(stmt, 1, machineName.c_str(), -1, SQLITE_STATIC);

        if (sqlite3_step(stmt) == SQLITE_ROW) {
            mappingId = sqlite3_column_int(stmt, 0);  // Get the id from the query result
        }
        sqlite3_finalize(stmt);
    } else {
        std::cerr << "Failed to prepare statement: " << sqlite3_errmsg(db) << std::endl;
    }
    return mappingId;
}

void setMetadata(sqlite3 *db, const std::string &key,
                 const std::string &value) {


  std::string sql = "INSERT INTO metadata (key, value) VALUES ('" + key +
                    "', '" + value + "')";

    executeSQL(db, sql, "Metadata");
}

void insertMetadata(sqlite3 *db, char *mpi_lib, int size,
                    char *cmd[MAX_ARGS], int ac, int mpisee_major_v, int mpisee_minor_v,
                    char *build_date, char *build_time, const char *env) {
  std::string lib = mpi_lib;
  std::string sz = std::to_string(size);
  std::string mpisee_version =
      std::to_string(mpisee_major_v) + "." + std::to_string(mpisee_minor_v);

  std::string mpisee_date =
      std::string(build_date) + ", " + std::string(build_time);

  std::string env_var = env;
  std::string combinedString;
  for (int i = 0; i < ac && i< MAX_ARGS && cmd[i] != NULL; ++i) {
    if (i > 0) {
      combinedString += " ";
    }
    combinedString += cmd[i];
  }

  auto now = std::chrono::system_clock::now();

  // Convert time_point to time_t for easier manipulation
  std::time_t now_time_t = std::chrono::system_clock::to_time_t(now);

  // Convert time_t to string representation
  std::string profile_date = std::ctime(&now_time_t);

  setMetadata(db, "MPI Library", lib);
  setMetadata(db, "Processes", sz);
  setMetadata(db, "Run command", combinedString);
  setMetadata(db, "mpisee Version", mpisee_date);
  setMetadata(db, "Profile date", profile_date);

}


void printMetadata(sqlite3* db) {
    sqlite3_stmt* stmt;
    std::string sql = "SELECT key, value FROM metadata";
    std::cout << "Printing Metadata\n";
    if (sqlite3_prepare_v2(db, sql.c_str(), -1, &stmt, nullptr) == SQLITE_OK) {
      while (sqlite3_step(stmt) == SQLITE_ROW) {
        const char *key =
            reinterpret_cast<const char *>(sqlite3_column_text(stmt, 0));

        const char *value =
            reinterpret_cast<const char *>(sqlite3_column_text(stmt, 1));
        std::cout << key << ": " << value << '\n';
      }
    }
    sqlite3_finalize(stmt);
}

void createTables(sqlite3* db) {

    // Create Mapping Table
    const char* MappingTable =
        "CREATE TABLE IF NOT EXISTS mappings ("
        "id INTEGER PRIMARY KEY AUTOINCREMENT, "
        "machine TEXT);";
    executeSQL(db, MappingTable, "Mappings Table created");

    // Create Metadata Table
    const char* Metadata =
        "CREATE TABLE IF NOT EXISTS metadata ("
        "key TEXT PRIMARY KEY, "
        "value TEXT);";
    executeSQL(db, Metadata, "Metadata Table created");

    // Create MPI Operations Table
    const char* MPIOpsTable =
        "CREATE TABLE IF NOT EXISTS operations ("
        "id INTEGER PRIMARY KEY AUTOINCREMENT, "
        "operation TEXT);";
    executeSQL(db, MPIOpsTable, "MPI Operations Table created ");

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
        "comm_id INTEGER, "
        "operation_id INTEGER, "
        "buffer_size_max INTEGER, "
        "buffer_size_min INTEGER, "
        "calls INTEGER, "
        "time REAL, "
        "FOREIGN KEY (operation_id) REFERENCES operations (id), "
        "FOREIGN KEY (comm_id) REFERENCES comms (id), "
        "FOREIGN KEY (rank) REFERENCES mappings (id));";
    executeSQL(db, DataTable, "Data Table created");
}

// Function to insert into mappings
void insertIntoMappings(sqlite3 *db, const std::string &machine) {
  int count;
  count = countTable(db, "mappings");
  std::string insertSql;
  if (count == 0) {
    // The table is empty, insert with id = 0
    insertSql = "INSERT INTO mappings (id, machine) VALUES (0, '" + machine + "')";
  } else {
    // The table is not empty, let SQLite auto-increment the id
    insertSql = "INSERT INTO mappings (machine) VALUES ('" + machine + "')";
  }

  executeSQL(db, insertSql, "INSERT INTO mappings");
}



// Functions to insert into comms
// void insertIntoComms(sqlite3 *db, const std::string &name,
//                      const std::string &size) {
//   int count;

//   count = countTable(db, "comms");
//   std::string insertSql;
//   if (count == 0) {
//     // The table is empty, insert with id = 0
//     insertSql = "INSERT INTO comms (id, name, size) VALUES (0, '" + name + "', '" + size + "')";
//   } else {
//     // The table is not empty, let SQLite auto-increment the id
//     sqlite3_stmt *stmt;
//     std::string checkSql = "SELECT COUNT(*) FROM comms WHERE name = '" + name + "'";
//     sqlite3_prepare_v2(db, checkSql.c_str(), -1, &stmt, nullptr);
//     sqlite3_bind_text(stmt, 1, name.c_str(), -1, SQLITE_STATIC);
//     count = -1;
//     count = sqlite3_column_int(stmt, 0);
//     // If the entry does not exist, insert it
//     if ( count == 0 ){
//       insertSql = "INSERT INTO comms (name, size) VALUES ('" + name + "', '" +
//                   size + "')";
//     } else {
//       std::cout << "Entry with name '" << name
//                 << "' already exists in comms table." << std::endl;
//       return;
//     }
//   }
//     executeSQL(db, insertSql,"INSERT INTO comms" );
// }

void insertIntoComms(sqlite3 *db, const std::string &name, const std::string &size) {
    sqlite3_stmt *stmt;
    std::string checkSql = "SELECT COUNT(*) FROM comms WHERE name = ?";
    sqlite3_prepare_v2(db, checkSql.c_str(), -1, &stmt, nullptr);
    sqlite3_bind_text(stmt, 1, name.c_str(), -1, SQLITE_STATIC);

    int count = -1;
    if (sqlite3_step(stmt) == SQLITE_ROW) {
        count = sqlite3_column_int(stmt, 0);
    }
    sqlite3_finalize(stmt);

    // If the entry does not exist, insert it
    if (count == 0) {
        std::string insertSql;
        // Check if it's the first entry to set id to 0
        if (countTable(db, "comms") == 0) {
            insertSql = "INSERT INTO comms (id, name, size) VALUES (0, '" + name + "', '" + size + "')";
        } else {
            insertSql = "INSERT INTO comms (name, size) VALUES ('" + name + "', '" + size + "')";
        }
        executeSQL(db, insertSql, "INSERT INTO comms");
    } else {
        std::cout << "Entry with name '" << name << "' already exists in comms table." << std::endl;
    }
}


// Functions to insert into operations
void insertIntoOperations(sqlite3* db, const std::string& operation) {
  int count;
  count = countTable(db, "operations");
  std::string insertSql;
  if (count == 0) {
    insertSql = "INSERT INTO operations (id, operation) VALUES (0, '" + operation + "')";
  } else {
    insertSql = "INSERT INTO operations (operation) VALUES ('" + operation + "')";

  }
  executeSQL(db, insertSql, "INSERT INTO operations");
}

// Functions to insert into data
void insertIntoData(sqlite3* db, int rank, int commId, int operationId, int bufferSizeMax, int bufferSizeMin, int calls, double time) {
  int count;
  count = countTable(db, "operations");
  std::string insertSql;
  if (count == 0) {
    insertSql = "INSERT INTO data (id, rank, comm_id, operation_id, buffer_size_max, buffer_size_min, calls, time) VALUES (0, "
                      + std::to_string(rank) + ", " + std::to_string(commId) + ", "
                      + std::to_string(operationId) + ", " + std::to_string(bufferSizeMax) + ", "
                      + std::to_string(bufferSizeMin) + ", " + std::to_string(calls) + ", " + std::to_string(time) + ")";

  } else {
    insertSql = "INSERT INTO data (rank, comm_id, operation_id, buffer_size_max, buffer_size_min, calls, time) VALUES ("
                      + std::to_string(rank) + ", " + std::to_string(commId) + ", "
                      + std::to_string(operationId) + ", " + std::to_string(bufferSizeMax) + ", "
                      + std::to_string(bufferSizeMin) + ", " + std::to_string(calls) + ", " + std::to_string(time) + ")";
  }
  executeSQL(db, insertSql, "INSERT INTO data");
}

/*
int main(int argc, char* argv[]) {
  sqlite3 *db;
  int rc,i;

  // Open database
  rc = sqlite3_open("test_mpi_data.db", &db);
  if (rc) {
      std::cerr << "Can't open database: " << sqlite3_errmsg(db) << std::endl;
      return 1;
  } else {
      std::cout << "Opened database successfully" << std::endl;
  }

  createTables(db);

  setMetadata(db, "MPI Library", "OpenMPI 4.1.4");
  setMetadata(db, "Run command", "./miniAMR 4 4 2 10 2 5");
  // Get current time as a time_point
  auto now = std::chrono::system_clock::now();

  // Convert time_point to time_t for easier manipulation
  std::time_t now_time_t = std::chrono::system_clock::to_time_t(now);

  // Convert time_t to string representation
  std::string now_str = std::ctime(&now_time_t);

  // Print the date and time
  std::cout << now_str; // This will print the date and time in a standard format
  setMetadata(db, "Profile Date", now_str);

  insertIntoMappings(db, "machine1");
  insertIntoMappings(db, "machine2");

  // Insert test data into comms
  std::string world = "W";
  std::string comm1 = "W1";
  std::string comm2 = "W2";
  insertIntoComms(db, world, "2");
  insertIntoComms(db, comm1, "1");
  insertIntoComms(db, comm2, "1");
  insertIntoComms(db, world, "2");

  // Insert test data into operations
  for (i = 0; i<NUM_OF_PRIMS ; i++ ) {
    insertIntoOperations(db, prim_names[i]);
  }
  // insertIntoOperations(db, "MPI_Send");
  // insertIntoOperations(db, "MPI_Recv");
  // insertIntoOperations(db, "MPI_Allreduce");
  // insertIntoOperations(db, "MPI_Bcast");

  // Insert test data into data for rank 0
  insertIntoData(db, 0, 0, Send, 128, 64, 10, 0.001); // MPI_Send
  insertIntoData(db, 0, 0, Recv, 256, 128, 5, 0.002); // MPI_Recv
  insertIntoData(db, 0, 1, Allreduce, 512, 256, 8, 0.003); // MPI_Allreduce
  insertIntoData(db, 0, 1, Bcast, 1024, 512, 3, 0.004); // MPI_Bcast

  // Insert test data into data for rank 1
  insertIntoData(db, 1, 0, Send, 128, 64, 15, 0.005); // MPI_Send
  insertIntoData(db, 1, 0, Recv, 256, 128, 7, 0.006); // MPI_Recv
  insertIntoData(db, 1, 2, Allreduce, 512, 256, 10, 0.007); // MPI_Allreduce
  insertIntoData(db, 1, 2, Bcast, 1024, 512, 5, 0.008); // MPI_Bcast

  // printDataDetails(db);
  printMetadata(db);
  printData(db);

  // std::cout << getCommId(db, world) << '\n';
  // char *machname = "machine2";
  // std::cout << getMappingId(db, machname) << '\n';

  sqlite3_close(db);
  return 0;
}
*/
