#ifndef SQLITEDB_H_
#define SQLITEDB_H_
#include <sqlite3.h>
#include <string>
#include "utils.h"

void printDataDetails(sqlite3 *db);

void setMetadata(sqlite3 *db, const std::string &key, const std::string &value);

void printMetadata(sqlite3 *db);

void createTables(sqlite3 *db);

void insertIntoMappings(sqlite3 *db, const std::string &machine);

void insertIntoComms(sqlite3 *db, const std::string &name, int size);

void insertIntoOperations(sqlite3 *db, const std::string &operation);

void insertIntoData(sqlite3 *db, int rank, int commId, int operationId,
                    int bufferSizeMax, int bufferSizeMin, int calls,
                    double time);

void insertMetadata(sqlite3 *db, char *mpi_lib, int size, char *cmd[MAX_ARGS],
                    int ac, int mpisee_major_v, int mpisee_minor_v,
                    char *build_date, char *build_time, const char *env);

int getCommId(sqlite3 *db, const std::string &commName);

void printData(sqlite3* db);

void printCommsTable(sqlite3 *db);

#endif
