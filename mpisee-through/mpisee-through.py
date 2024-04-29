#!/usr/bin/env python3
import argparse
import sqlite3
import re
import os
import sys
import csv


from signal import signal, SIGPIPE, SIG_DFL

signal(SIGPIPE, SIG_DFL)

# Global variables
RED = "\033[1;31m"
BLUE = "\033[1;34m"
CYAN = "\033[1;36m"
GREEN = "\033[0;32m"
RESET = "\033[0;0m"
BOLD = "\033[;1m"
REVERSE = "\033[;7m"



# Function declarations
def print_decoration(decoration):
    if sys.stdout.isatty():
        # the output is not redirected, we can use a fancy style:
        sys.stdout.write(decoration)

def is_almost_equal(a, b, epsilon=1e-6):
    return abs(a - b) <= epsilon

def parse_enum_from_header(header_path):
    with open(header_path, 'r') as file:
        content = file.read()

    # Regular expression to find enum definition
    enum_pattern = re.compile(r'enum primitives\{([^}]+)\};', re.MULTILINE | re.DOTALL)
    match = enum_pattern.search(content)

    enum_dict = {}
    if match:
        enum_content = match.group(1)
        # Extract individual enum items
        enum_items = enum_content.split(',')
        value = 0  # Assuming enum starts at 0
        for item in enum_items:
            item = item.strip()
            if item:  # Non-empty string
                if '=' in item:
                    name, val = item.split('=')
                    value = int(val.strip())
                    enum_dict[name.strip()] = value
                else:
                    enum_dict[item] = value
                value += 1
    else:
        print("Enum 'primitives' not found in the header file.")
    return enum_dict

def get_exec_time_by_rank(db_path, rank):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    total_time = 0.0

    try:
        sql = "SELECT time FROM exectimes WHERE id = ?"

        cursor.execute(sql, (rank,))
        result = cursor.fetchone()

        if result:
            total_time = result[0]
            return total_time
        else:
            print(f"No data found for rank {rank}.")
            return None

    except sqlite3.Error as e:
        print("An error occurred:", e)
        return None
    finally:
        conn.close()


def get_mpi_time_by_rank(db_path, rank):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    try:
        sql = "SELECT total_time FROM mpi_time_sum WHERE rank = ?"

        cursor.execute(sql, (rank,))
        result = cursor.fetchone()

        if result:
            total_time = result[0]
            return total_time
        else:
            print(f"No data found for rank {rank}.")
            return None

    except sqlite3.Error as e:
        print("An error occurred:", e)
        return None
    finally:
        conn.close()


def query_data_by_rank(db_path, rank):
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT * FROM data WHERE rank = ?", (rank,))
        return cursor.fetchall()

def exec_query_and_print(db_path,sql,order,num_of_rows,ranks,comms,*args):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    params = args
    if len(ranks) > 0:
        placeholders = ','.join('?' * len(ranks))
        sql += f" AND d.rank IN ({placeholders})"
        params += tuple(ranks)
    if len(comms) > 0:
        placeholders = ','.join('?' *len(comms))
        sql += f" AND c.name IN ({placeholders})"
        params += tuple(comms)

    sql = select_order(sql,order)

    try:
        # Execute the query
        cursor.execute(sql,params)

        # Print header
        print_decoration(BOLD)
        print(f"{'Comm Name':<15}{'Comm Size':<15}{'Rank':<10}{'MPI Operation':<20}"
              f"{'Buffer Size (Bytes)':<30}{'Calls':<15}{'Time (s)':<15}{'% of MPI Time':<20}{'% of Total Time':<10}")
        print_decoration(RESET)
        # Print rows
        r = 0
        prev_rank = -1
        percentage_mpi_time = 0.0
        percentage_exec_time = 0.0
        exec_time = -1.0
        mpi_time = -1.0
        for row in cursor.fetchall():
            name, size, rank, operation, buf_min, buf_max, calls, time = row
            buffer_size = f"{buf_min} - {buf_max}"
            if rank != prev_rank:
                exec_time = get_exec_time_by_rank(db_path,rank)
                mpi_time = get_mpi_time_by_rank(db_path,rank)
                prev_rank = rank

            percentage_exec_time = (time/exec_time)*100
            percentage_mpi_time = (time/mpi_time)*100
            print(f"{name:<15}{size:<15}{rank:<10}{operation:<20}"
                  f"{buffer_size:<30}{calls:<15}{time:<15.3f}{percentage_mpi_time:<20.3f}{percentage_exec_time:<10.3f}")
            r+=1
            if num_of_rows > 0 and r >= num_of_rows:
                break

    except sqlite3.Error as e:
        print("Failed to read data from SQLite table", e)
    finally:
        # Close the database connection
        if conn:
            conn.close()

def select_order(sql,order):
    if order == 0:
        sql += """
        ORDER BY c.name"""
    elif order == 1:
        sql += """
        ORDER BY d.time DESC"""
    elif order == 2:
        sql += """
        ORDER BY d.time ASC"""
    elif order == 3:
        sql += """
        ORDER BY d.operation_id DESC"""
    elif order == 4:
        sql += """
        ORDER BY d.buffer_size_min DESC"""
    elif order == 5:
        sql += """
        ORDER BY d.buffer_size_min ASC"""
    elif order == 6:
        sql += """
        ORDER BY d.calls DESC"""
    elif order == 7:
        sql += """
        ORDER BY d.calls ASC"""
    return sql


def print_all_data(db_path):
    # Connect to the SQLite database
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # SQL query
    sql = """
    SELECT c.name, c.size, d.rank, o.operation, d.buffer_size_min, d.buffer_size_max,
           d.calls, d.time
    FROM data d
    JOIN comms c ON d.comm_id = c.id
    JOIN operations o ON d.operation_id = o.id
    ORDER BY c.name
    """

    try:
        # Execute the query
        cursor.execute(sql)

        # Print header
        print(f"{'Comm Name':<15}{'Comm Size':<15}{'Rank':<10}{'Operation':<20}"
              f"{'Buffer Size Range (Bytes)':<25}{'Calls':<15}{'Time':<20}")

        # Print rows
        for row in cursor.fetchall():
            name, size, rank, operation, buf_min, buf_max, calls, time = row
            buffer_size = f"{buf_min} - {buf_max}"
            print(f"{name:<15}{size:<15}{rank:<10}{operation:<20}"
                  f"{buffer_size:<25}{calls:<15}{time:<20}")

    except sqlite3.Error as e:
        print("Failed to read data from SQLite table", e)
    finally:
        # Close the database connection
        if conn:
            conn.close()

def print_data_by_rank(db_path, rank):
    # Connect to the SQLite database
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # SQL query
    sql = """
    SELECT c.name, c.size, d.rank, o.operation, d.buffer_size_min, d.buffer_size_max,
           d.calls, d.time
    FROM data d
    JOIN comms c ON d.comm_id = c.id
    JOIN operations o ON d.operation_id = o.id
    WHERE d.rank = ?
    ORDER BY c.name
    """

    try:
        # Execute the query
        cursor.execute(sql,(rank,))

        # Print header
        print(f"{'Comm Name':<15}{'Comm Size':<15}{'Rank':<10}{'Operation':<20}"
              f"{'Buffer Size Range':<25}{'Calls':<15}{'Time':<20}")

        # Print rows
        for row in cursor.fetchall():
            name, size, rank, operation, buf_min, buf_max, calls, time = row
            buffer_size = f"{buf_min} - {buf_max}"
            print(f"{name:<15}{size:<15}{rank:<10}{operation:<20}"
                  f"{buffer_size:<25}{calls:<15}{time:<20}")

    except sqlite3.Error as e:
        print("Failed to read data from SQLite table", e)
    finally:
        # Close the database connection
        if conn:
            conn.close()

def print_data_by_comm(db_path, comm):
    # Connect to the SQLite database
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # SQL query
    sql = """
    SELECT c.name, c.size, d.rank, o.operation, d.buffer_size_min, d.buffer_size_max,
           d.calls, d.time
    FROM data d
    JOIN comms c ON d.comm_id = c.id
    JOIN operations o ON d.operation_id = o.id
    WHERE c.name = ?
    ORDER BY d.time DESC
    """
    print(f"Data for communicator: {comm}")
    try:
        # Execute the query
        cursor.execute(sql,(comm))

        # Print header
        print(f"{'Comm Name':<15}{'Comm Size':<15}{'Rank':<10}{'Operation':<20}"
              f"{'Buffer Size Range':<25}{'Calls':<15}{'Time':<20}")

        # Print rows
        for row in cursor.fetchall():
            name, size, rank, operation, buf_min, buf_max, calls, time = row
            buffer_size = f"{buf_min} - {buf_max}"
            print(f"{name:<15}{size:<15}{rank:<10}{operation:<20}"
                  f"{buffer_size:<25}{calls:<15}{time:<20}")

    except sqlite3.Error as e:
        print("Failed to read data from SQLite table", e)
    finally:
        # Close the database connection
        if conn:
            conn.close()

def print_data_by_prim(db_path, MPI_prim,nresults=10):
    # Connect to the SQLite database
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # SQL query
    sql = """
    SELECT c.name, c.size, d.rank, o.operation, d.buffer_size_min, d.buffer_size_max,
           d.calls, d.time
    FROM data d
    JOIN comms c ON d.comm_id = c.id
    JOIN operations o ON d.operation_id = o.id
    WHERE o.operation = ?
    ORDER BY d.time DESC
    """
    j = 0
    try:
        # Execute the query
        cursor.execute(sql,(MPI_prim,))

        # Print header
        print(f"{'Comm Name':<15}{'Comm Size':<15}{'Rank':<10}{'Operation':<20}"
              f"{'Buffer Size Range':<25}{'Calls':<15}{'Time':<20}")

        # Print rows
        for row in cursor.fetchall():
            if (j >= nresults):
                break
            name, size, rank, operation, buf_min, buf_max, calls, time = row
            buffer_size = f"{buf_min} - {buf_max}"
            print(f"{name:<15}{size:<15}{rank:<10}{operation:<20}"
                  f"{buffer_size:<25}{calls:<15}{time:<20}")
            j+=1

    except sqlite3.Error as e:
        print("Failed to read data from SQLite table", e)
    finally:
        # Close the database connection
        if conn:
            conn.close()

def print_execution_time(dpath,ranks=[]):
    conn = sqlite3.connect(dpath)
    cursor = conn.cursor()
    sql_exectime = """
    SELECT t.id, t.time
    FROM exectimes t
    """
    params_exectime = ""
    if len(ranks) > 0:
        placeholders = ','.join('?' * len(ranks))
        sql_exectime += f"WHERE t.id IN ({placeholders})"
        params_exectime = tuple(ranks)

    # elif order == 1:
    #     sql_exectime += """
    #     ORDER BY t.time DESC"""
    # elif order == 2:
    #     sql_exectime+= """
    #     ORDER BY t.time ASC"""

    sql_exectime += """
    ORDER BY t.id ASC"""


    sql_mpitime = """
    SELECT rank, total_time as mpi_time
    FROM mpi_time_sum
    """
    params_mpitime = ""
    if len(ranks) > 0:
        placeholders = ','.join('?' * len(ranks))
        sql_mpitime += f" WHERE rank IN ({placeholders})"
        params_mpitime = tuple(ranks)

    sql_mpitime += f"GROUP BY rank "
    sql_mpitime += """
    ORDER BY rank ASC"""

    # if order == 1:
    #     sql_mpitime += f" ORDER BY mpi_time DESC"
    # else:
    #     sql_mpitime += f" ORDER BY mpi_time ASC"



    try:
        # Execute the query
        cursor.execute(sql_exectime,(params_exectime))
        exectimes = cursor.fetchall()
        cursor.execute(sql_mpitime,(params_mpitime))
        mpitimes = cursor.fetchall()

        # Print header
        print_decoration(BOLD)
        string_ranks = ",".join(str(num) for num in ranks)
        print(f"Time statistics for MPI ranks: {string_ranks}")
        print(f"{'MPI Rank':<10}{'MPI Time (s)':<15}{'Execution Time (s)':<20}{'MPI time to Execution Time Ratio (%)':<10}")
        print_decoration(RESET)

        # Print rows
        for i in range(len(exectimes)):
            id,time = exectimes[i]
            id_mpi,total_time = mpitimes[i]
            print(f"{id:<10}{total_time:<15.3f}{time:<20.3f}{(total_time/time)*100:.2f}")

    except sqlite3.Error as e:
        print("Failed to read data from SQLite exectime table:", e)
    finally:
        # Close the database connection
        if conn:
            conn.close()

def mpi_time(dbpath,order=1,ranks=[]):
    conn = sqlite3.connect(dbpath)
    cursor = conn.cursor()

    sql ="""
    SELECT rank, total_time as mpi_time
    FROM mpi_time_sum
    """
    params = ""
    if len(ranks) > 0:
        placeholders = ','.join('?' * len(ranks))
        sql += f" WHERE rank IN ({placeholders})"
        params = tuple(ranks)

    sql += f"GROUP BY rank "

    if order == 1:
        sql += f" ORDER BY mpi_time DESC"
    else:
        sql += f" ORDER BY mpi_time ASC"

    try:
        cursor.execute(sql,params)
        rows = cursor.fetchall()
        print_decoration(BOLD)
        print(f"{'Rank':<8}{'MPI Time':<15}")
        print_decoration(RESET)
        for row in rows:
            rank, total_time = row
            print(f"{rank:<10}{total_time:.3f}")
    except sqlite3.Error as e:
        print("Failed to read data from SQLite exectime table:", e)
    finally:
        conn.close()

def print_data_by_time(dbpath,order=1,num_of_rows=0,rank_list=[],comms=[],*args):
    # SQL query
    sql = """
    SELECT c.name, c.size, d.rank, o.operation, d.buffer_size_min, d.buffer_size_max,
           d.calls, d.time
    FROM data d
    JOIN comms c ON d.comm_id = c.id
    JOIN operations o ON d.operation_id = o.id
    WHERE d.time >= ? AND d.time <= ?
    """
    exec_query_and_print(dbpath,sql,order,num_of_rows,rank_list,comms,*args)

def print_data_by_bufsize(dbpath,order=1,num_of_rows=0,rank_list=[],comms=[],*args):
    # SQL query
    sql = """
    SELECT c.name, c.size, d.rank, o.operation, d.buffer_size_min, d.buffer_size_max,
           d.calls, d.time
    FROM data d
    JOIN comms c ON d.comm_id = c.id
    JOIN operations o ON d.operation_id = o.id
    WHERE d.buffer_size_min >= ? AND d.buffer_size_max <= ?
    """
    exec_query_and_print(dbpath,sql,order,num_of_rows,rank_list,comms,*args)

def  print_data_collectives(dbpath,order=1,num_of_rows=0,rank_list=[],comms=[],*args):
    sql = """
        SELECT c.name, c.size, d.rank, o.operation, d.buffer_size_min, d.buffer_size_max,
           d.calls, d.time
        FROM data d
        JOIN comms c ON d.comm_id = c.id
        JOIN operations o ON d.operation_id = o.id
        WHERE d.buffer_size_min >= ? AND d.buffer_size_max <= ? AND d.operation_id >= ? """
    exec_query_and_print(dbpath,sql,order,num_of_rows,rank_list,comms,*args)


def print_data_pt2pt(dbpath,order=1,num_of_rows=0,rank_list=[],comms=[],*args):
    sql = """
        SELECT c.name, c.size, d.rank, o.operation, d.buffer_size_min, d.buffer_size_max,
           d.calls, d.time
        FROM data d
        JOIN comms c ON d.comm_id = c.id
        JOIN operations o ON d.operation_id = o.id
        WHERE d.buffer_size_min >= ? AND d.buffer_size_max <= ? AND d.operation_id <= ?  """
    exec_query_and_print(dbpath,sql,order,num_of_rows,rank_list,comms,*args)

def query_all_data(dbpath,order=1,num_of_rows=0,rank_list=[],comms=[],*args):
    sql = """SELECT c.name, c.size, d.rank, o.operation,
                      d.buffer_size_min, d.buffer_size_max, d.calls, d.time
                      FROM data d
                      JOIN comms c ON d.comm_id = c.id
                      JOIN operations o ON d.operation_id = o.id """
    exec_query_and_print(dbpath,sql,order,num_of_rows,rank_list,comms,*args)


def sort_by(results,order):
    if ( order == 0):
        results = dict(sorted(results.items(), key=lambda item: item[1]['comm_name']))
    elif ( order == 1):
        results = dict(sorted(results.items(), key=lambda item: item[1]['time_s'], reverse=True))
    elif ( order == 2):
        results = dict(sorted(results.items(), key=lambda item: item[1]['time_s']))
    elif ( order == 3):
        results = dict(sorted(results.items(), key=lambda item: item[1]['operation']))
    elif ( order == 4):
        results = dict(sorted(results.items(), key=lambda item: item[1]['buffer_size_min'], reverse=True))
    elif ( order == 5):
        results = dict(sorted(results.items(), key=lambda item: item[1]['buffer_size_min']))
    elif ( order == 6):
        results = dict(sorted(results.items(), key=lambda item: item[1]['calls'], reverse=True))
    elif ( order == 7):
        results = dict(sorted(results.items(), key=lambda item: item[1]['calls']))
    return results


def default_query(dbpath,enum_primitives,order=1):
    sql = """
    SELECT
    d.comm_id AS cid,
    c.name AS comm_name,
    c.size AS comm_size,
    d.rank,
    d.operation_id AS opid,
    o.operation,
    d.buffer_size_min,
    d.buffer_size_max,
    SUM(d.calls) AS calls,
    MAX(d.time) AS time_s,
    AVG(d.time) AS avg_time,
    SUM(volume) As total_volume
FROM data d
JOIN comms c ON d.comm_id = c.id
JOIN operations o ON d.operation_id = o.id
GROUP BY c.name, c.size, o.operation, d.buffer_size_min, d.buffer_size_max
        """

    conn = sqlite3.connect(dbpath)
    cursor = conn.cursor()
    try:
        # Execute the query
        cursor.execute(sql)

        # Print header
        print_decoration(BOLD)
        print(f"{'Comm Name':<15}{'Processes':<25}{'Comm Size':<12}{'MPI Operation':<15}"
              f"{'Min Buffer':<12}{'Max Buffer':<12}{'Calls':<12}{'Max Time(s)':<13}{'Avg Time(s)':<13}{'Total Volume(Bytes)':<15}")
        print_decoration(RESET)

        data = cursor.fetchall()  # Retrieve all data

        results = {}

        for row in data:
            cid,comm_name, comm_size, rank, opid, operation, buf_min, buf_max, calls, time, avg_time, volume = row
            key = (cid, comm_name, comm_size, opid, operation, buf_min, buf_max)
            if key not in results:
                results[key] = {
                'cid': cid,
                'comm_name': comm_name,
                'comm_size': comm_size,
                'rank': rank,
                'opid': opid,
                'operation': operation,
                'buffer_size_min': buf_min,
                'buffer_size_max': buf_max,
                'calls': calls,
                'time_s': time,
                'avg_time': avg_time,
                'total_volume': volume

            }

        results = sort_by(results,order)

        for result in results.values():
            calls = result['calls']
            if result['opid'] == enum_primitives['Sendrecv']:
                calls = calls // 2
            elif result['opid'] >= enum_primitives['Bcast']:
                calls = calls // result['comm_size']


            procs = list_truncate_tostr(get_ranks_by_comm(dbpath,result['cid']))
            print(f"{result['comm_name']:<15}{procs:<25}{result['comm_size']:<12}{result['operation']:<15}"
                  f"{result['buffer_size_min']:<12}{result['buffer_size_max']:<12}{calls:<12}{result['time_s']:<13.3f}{result['avg_time']:<13.3f}{result['total_volume']}")  # Adjust formatting as needed

    except sqlite3.Error as e:
        print("Failed to read data from SQLite table", e)
    finally:
        # Close the database connection
        if conn:
            conn.close()



def query_colls_pt2pt(dbpath,enum_primitives,mpi_op,bufmin,bufmax,tmin,tmax,order=1,outfile=None):

    main_query = """
SELECT
  d.comm_id AS cid,
  c.name AS comm_name,
  c.size AS comm_size,
  d.rank,
  d.operation_id AS opid,
  o.operation,
  d.buffer_size_min,
  d.buffer_size_max,
  SUM(d.calls) AS calls,
  MAX(d.time) AS time_s,
  AVG(d.time) AS avg_time,
  SUM(volume) AS total_volume
FROM data d
JOIN comms c ON d.comm_id = c.id
JOIN operations o ON d.operation_id = o.id
"""

    if mpi_op == 'Bcast':
        where_clause = "WHERE d.operation_id >= ?"
    elif mpi_op == 'Ibsend':
        where_clause = "WHERE d.operation_id <= ?"
    else:
        where_clause = ""  # No WHERE clause needed

    group_query = " GROUP BY c.name, c.size, o.operation, d.buffer_size_min, d.buffer_size_max"


    sql = main_query + " " + where_clause + group_query


    conn = sqlite3.connect(dbpath)
    cursor = conn.cursor()
    try:
        # Execute the query
        if ( mpi_op != None):
            cursor.execute(sql,(enum_primitives[mpi_op],))
        else:
            cursor.execute(sql)
        # Print header

        if outfile:
              csvfile=open(outfile, 'w', newline='')
              csv_writer = csv.writer(csvfile)
              # Write header row
              csv_writer.writerow(['Comm Name', 'Processes', 'Comm Size', 'MPI Operation','Min Buffer', 'Max Buffer', 'Calls', 'Max Time(s)', 'Avg Time(s)','Volume(Bytes)'])
        else:
            print_decoration(BOLD)
            print(f"{'Comm Name':<15}{'Processes':<25}{'Comm Size':<12}{'MPI Operation':<15}"
                  f"{'Min Buffer':<12}{'Max Buffer':<12}{'Calls':<12}{'Max Time(s)':<13}{'Avg Time(s)':<13}{'Total Volume(Bytes)':<15}")
            print_decoration(RESET)

        data = cursor.fetchall()  # Retrieve all data

        results = {}

        for row in data:
            cid,comm_name, comm_size, rank, opid, operation, buf_min, buf_max, calls, time, avg_time, volume = row
            key = (cid, comm_name, comm_size, opid, operation, buf_min, buf_max)
            if key not in results:
                results[key] = {
                'cid': cid,
                'comm_name': comm_name,
                'comm_size': comm_size,
                'rank': rank,
                'opid': opid,
                'operation': operation,
                'buffer_size_min': buf_min,
                'buffer_size_max': buf_max,
                'calls': calls,
                'time_s': time,
                'avg_time': avg_time,
                'total_volume': volume

            }

        results = sort_by(results,order)

        for result in results.values():
            calls = result['calls']
            if result['opid'] == enum_primitives['Sendrecv']:
                calls = calls // 2
            elif result['opid'] >= enum_primitives['Bcast']:
                calls = calls // result['comm_size']

            if (result['buffer_size_min'] < bufmin or result['buffer_size_max'] > bufmax):
                continue
            if (result['time_s'] < tmin or result['time_s'] > tmax): # Fix this comparison
                continue
            procs = list_truncate_tostr(get_ranks_by_comm(dbpath,result['cid']))
            if outfile:
                # Write output to csv file instead
                csv_writer.writerow([result['comm_name'], procs,
                                     result['comm_size'], result['operation'],
                                     result['buffer_size_min'],
                                     result['buffer_size_max'], calls,
                                     "{:.3f}".format(result['time_s']),
                                     "{:.3f}".format(result['avg_time']),
                                     result['total_volume']])
            else:
                print(f"{result['comm_name']:<15}{procs:<25}{result['comm_size']:<12}{result['operation']:<15}"
                      f"{result['buffer_size_min']:<12}{result['buffer_size_max']:<12}{calls:<12}{result['time_s']:<13.3f}{result['avg_time']:<13.3f}{result['total_volume']}")  # Adjust formatting as needed

    except sqlite3.Error as e:
        print("Failed to read data from SQLite table", e)
    finally:
        # Close the database connection
        if conn:
            conn.close()
        if outfile:
            csvfile.close()


def query_ranks(dbpath,enum_primitives,mpi_op,bufmin,bufmax,tmin,tmax,ranks,order=1,outfile=None):

    main_query = """
SELECT
  d.comm_id AS cid,
  c.name AS comm_name,
  c.size AS comm_size,
  d.rank,
  d.operation_id AS opid,
  o.operation,
  d.buffer_size_min,
  d.buffer_size_max,
  d.calls AS calls,
  d.time AS time_s,
  volume AS total_volume
FROM data d
JOIN comms c ON d.comm_id = c.id
JOIN operations o ON d.operation_id = o.id
"""

    if mpi_op == 'Bcast':
        where_clause = "WHERE d.operation_id >= ?"
    elif mpi_op == 'Ibsend':
        where_clause = "WHERE d.operation_id <= ?"
    else:
        where_clause = ""  # No WHERE clause needed

    sql = main_query + " " + where_clause


    conn = sqlite3.connect(dbpath)
    cursor = conn.cursor()
    try:
        # Execute the query
        if ( mpi_op != None):
            cursor.execute(sql,(enum_primitives[mpi_op],))
        else:
            cursor.execute(sql)
        # Print header
        if outfile:
              csvfile=open(outfile, 'w', newline='')
              csv_writer = csv.writer(csvfile)
              # Write header row
              csv_writer.writerow(['Comm Name', 'Processes', 'Comm Size', 'Rank', 'MPI Operation',
                                   'Min Buffer', 'Max Buffer', 'Calls', 'Time(s)', 'Volume(Bytes)'])
        else:
            print_decoration(BOLD)
            print(f"{'Comm Name':<15}{'Processes':<25}{'Comm Size':<12}{'Rank':<10}{'MPI Operation':<15}"
              f"{'Min Buffer':<12}{'Max Buffer':<12}{'Calls':<12}{'Time(s)':<13}{'Volume(Bytes)':<15}")
            print_decoration(RESET)

        data = cursor.fetchall()  # Retrieve all data

        results = {}

        for row in data:
            cid,comm_name, comm_size, rank, opid, operation, buf_min, buf_max, calls, time, volume = row
            if rank not in ranks:
                continue
            key = (cid, comm_name, comm_size, opid, operation, rank, buf_min, buf_max)
            if key not in results:
                results[key] = {
                'cid': cid,
                'comm_name': comm_name,
                'comm_size': comm_size,
                'rank': rank,
                'opid': opid,
                'operation': operation,
                'buffer_size_min': buf_min,
                'buffer_size_max': buf_max,
                'calls': calls,
                'time_s': time,
                'total_volume': volume

            }

        results = sort_by(results,order)

        for result in results.values():
            calls = result['calls']
            if result['opid'] == enum_primitives['Sendrecv']:
                calls = calls // 2

            if (result['buffer_size_min'] < bufmin or result['buffer_size_max'] > bufmax):
                continue
            if (result['time_s'] < tmin or result['time_s'] > tmax): # Fix this comparison
                 continue
            procs = list_truncate_tostr(get_ranks_by_comm(dbpath,result['cid']))
            if outfile:
                 csv_writer.writerow([result['comm_name'], procs, result['comm_size'], result['rank'], result['operation'],
                                      result['buffer_size_min'], result['buffer_size_max'], result['calls'], "{:.3f}".format(result['time_s']), result['total_volume']])
            else:
                print(f"{result['comm_name']:<15}{procs:<25}{result['comm_size']:<12}{result['rank']:<10}{result['operation']:<15}"
                  f"{result['buffer_size_min']:<12}{result['buffer_size_max']:<12}{calls:<12}{result['time_s']:<13.3f}{result['total_volume']}")  # Adjust formatting as needed

    except sqlite3.Error as e:
        print("Failed to read data from SQLite table", e)
    finally:
        # Close the database connection
        if conn:
            conn.close()
        if outfile:
            csvfile.close()

def query_summarize_time(dbpath,order=1,outfile=None):
    sql = """
    SELECT
    c.name AS comm_name,
    c.size AS comm_size,
    SUM(d.time) AS total_time_s
FROM data d
JOIN comms c ON d.comm_id = c.id
JOIN operations o ON d.operation_id = o.id
GROUP BY c.name, c.size
    """

    if order == 0:
        sql += """
        ORDER BY c.name"""
    elif order == 1:
        sql += """
        ORDER BY total_time_s DESC"""
    elif order == 2:
        sql += """
        ORDER BY total_time_s ASC"""
    elif order == 3:
        sql += """
        ORDER BY d.operation_id DESC"""
    elif order == 4:
        sql += """
        ORDER BY d.buffer_size_min DESC"""
    elif order == 5:
        sql += """
        ORDER BY d.buffer_size_min ASC"""
    elif order == 6:
        sql += """
        ORDER BY d.calls DESC"""
    elif order == 7:
        sql += """
        ORDER BY d.calls ASC"""

    conn = sqlite3.connect(dbpath)
    cursor = conn.cursor()
    try:
        cursor.execute(sql)
        # Print header
        if outfile:
            # Write data to csv file
            csvfile = open(outfile, 'w', newline='')
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(['Comm Name', 'Comm Size', 'Total Time(s)'])
        else:
            print_decoration(BOLD)
            print(f"{'Comm Name':<15}{'Comm Size':<15}{'Total Time (s)':<20}")
            print_decoration(RESET)

        data = cursor.fetchall()  # Retrieve all data

        for row in data:
            comm_name, comm_size, total_time = row
            if outfile:
                csv_writer.writerow([comm_name, comm_size, "{:.4f}".format(total_time)])
            else:
                print(f"{comm_name:<15}{comm_size:<15}{total_time:<20.4f}")

    except sqlite3.Error as e:
        print("Failed to read data from SQLite table", e)
    finally:
        # Close the database connection
        if conn:
            conn.close()
        if outfile:
            csvfile.close()



def default_query_summary(dbpath,MPI_prim):
    sql = """
        SELECT
        c.name AS comm_name,
        c.size AS comm_size,
        o.operation,
        d.buffer_size_min,
        d.buffer_size_max,
        SUM(d.calls) AS calls,
        MAX(d.time) AS time_s
        FROM data d
        JOIN comms c ON d.comm_id = c.id
        JOIN operations o ON d.operation_id = o.id
        WHERE o.operation = ?
        GROUP BY c.name, c.size, d.rank, o.operation, d.buffer_size_min, d.buffer_size_max;
        """

    conn = sqlite3.connect(dbpath)
    cursor = conn.cursor()
    time = 0.0
    try:
        cursor.execute(sql,(MPI_prim,))

        # Print header
        print_decoration(BOLD)
        print(f"{'Comm Name':<15}{'Comm Size':<15}{'MPI Operation':<20}"
              f"{'Buffer Size (Bytes)':<25}{'Calls':<15}{'Time (s)':<15}")
        print_decoration(RESET)

        data = cursor.fetchall()  # Retrieve all data

        results = {}

        for row in data:
            comm_name, comm_size, operation, buf_min, buf_max, calls, time = row
            key = (comm_name, comm_size, operation, buf_min, buf_max)
            if key not in results or results[key]['time_s'] < time:
                results[key] = {
                'comm_name': comm_name,
                'comm_size': comm_size,
                'operation': operation,
                'buffer_size_min': buf_min,
                'buffer_size_max': buf_max,
                'calls': calls,
                'time_s': time
            }
        results = dict(sorted(results.items(), key=lambda item: item[1]['time_s'], reverse=True))
        for result in results.values():
            time += result['time_s']
            buffer_size = f"{result['buffer_size_min']} - {result['buffer_size_max']}"
            print(f"{result['comm_name']:<15}{result['comm_size']:<15}{result['operation']:<20}"
                  f"{buffer_size:<25}{result['calls']:<15}{result['time_s']:<15.3f}")


    except sqlite3.Error as e:
        print("Failed to read data from SQLite table", e)
    finally:
        # Close the database connection
        print_decoration(BOLD)
        print(f"\nTotal time for {MPI_prim}: {time:.3f} s")
        print_decoration(RESET)
        if conn:
            conn.close()

def clear_table_if_exists(db_path, table_name):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    try:
        # Check if the table exists
        cursor.execute("SELECT count(name) FROM sqlite_master WHERE type='table' AND name=?", (table_name,))
        if cursor.fetchone()[0] == 1:
            # Table exists, so clear it
            cursor.execute(f"DELETE FROM {table_name}")
            conn.commit()
            #print(f"Table '{table_name}' cleared.")
        #else:
        #    print(f"Table '{table_name}' does not exist.")

    except sqlite3.Error as e:
        print("An error occurred:", e)
    finally:
        conn.close()

def create_and_populate_summary_table(db_path):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    clear_table_if_exists(db_path,"mpi_time_sum")

    # Step 1: Create a new table
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS mpi_time_sum (
            rank INTEGER PRIMARY KEY,
            total_time REAL
        )
    """)

    # Step 2: Aggregate and insert data
    cursor.execute("""
        INSERT INTO mpi_time_sum (rank, total_time)
        SELECT d.rank, SUM(d.time) as total_time
        FROM data d
        GROUP BY d.rank
    """)

    conn.commit()
    conn.close()

def print_comms_table(db_path):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute("SELECT * FROM comms")

    rows = cursor.fetchall()
    print_decoration(BLUE)
    print(f"{'ID':<5}{'Name':<15}{'Size':<15}")
    print_decoration(RESET)
    for row in rows:
        id, name, size = row
        print(f"{id:<5}{name:<15}{size:<15}")

    conn.close()

def print_metadata_table(db_path):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute("SELECT key, value FROM metadata")

    rows = cursor.fetchall()
    print_decoration(BLUE)
    for row in rows:
        key, value = row
        print(f"{key}: {value}")

    print_decoration(RESET)
    conn.close()

def print_operations_table(db_path):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute("SELECT * FROM operations")

    rows = cursor.fetchall()
    print_decoration(BLUE)
    print(f"{'ID':<5}{'Operation':<15}")
    print_decoration(RESET)
    for row in rows:
        id, operation = row
        print(f"{id:<5}{operation:<15}")

    conn.close()

def print_data_table(db_path):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute("SELECT * FROM data")

    rows = cursor.fetchall()
    print_decoration(BLUE)
    print(f"{'ID':<5}{'Rank':<5}{'Comm ID':<10}{'Operation ID':<15}{'Buffer Size Min':<20}"
          f"{'Buffer Size Max':<20}{'Calls':<10}{'Time':<10}")
    print_decoration(RESET)
    for row in rows:
        id, rank, comm_id, operation_id, buf_min, buf_max, calls, time = row
        print(f"{id:<5}{rank:<5}{comm_id:<10}{operation_id:<15}{buf_min:<20}{buf_max:<20}{calls:<10}{time:<10}")

    conn.close()


def get_max_time_rank(db_path,sql):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    max_time = -1.0
    rank = -1
    try:
        cursor.execute(sql)
        result = cursor.fetchone()

        if result:
            rank, max_time = result
        else:
            print("No data found in table.")

    except sqlite3.Error as e:
        print("An error occurred:", e)
    finally:
        conn.close()
    return rank,max_time

def get_avg_time(db_path,sql):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    time = -1.0
    try:
        cursor.execute(sql)
        result = cursor.fetchone()[0]

        if result:
            time = result
        else:
            print("No data found in table.")

    except sqlite3.Error as e:
        print("An error occurred:", e)
    finally:
        conn.close()
    return time

def get_all_times(db_path,sql):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    all_times_dict = {}
    try:
        cursor.execute(sql)
        for row in cursor.fetchall():
            all_times_dict[row[0]] = row[1]

    except sqlite3.Error as e:
        print("An error occurred:", e)
    finally:
        conn.close()

    return all_times_dict

def max_value_in_dict(d):
    if not d:
        return None,None

    max_key = max(d, key=lambda k: d[k])

    return max_key,d[max_key]

def avg_value_in_dict(d):
    average = -1.0
    if not d:
        return None

    total = sum(d.values())
    average = total/len(d)

    return average

def dict_ratios(dict_mpi,dict_exec):
    #dict_mpi is dictionary with the MPI times
    #dict_exec is the dictionary of Execution times
    #ratios will contain the MPI to Execution time Ratio for every MPI rank

    ratios = {}
    # Keys of both dictionaries must be the same
    for k in dict_mpi.keys():
        ratios[k] = float((dict_mpi[k] / dict_exec[k]))*100

    return ratios

def print_general_stats(db_path):
    size = -1
    sql = """
    SELECT value
    FROM metadata
    WHERE key = 'Processes'
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    try:
        cursor.execute(sql)
        result = cursor.fetchone()

        if result:
            size = int(result[0])
        else:
            print("No data found in tabe.")

    except sqlite3.Error as e:
        print("An error occurred:", e)
    finally:
        conn.close()

    print_decoration(GREEN)
    print("Overall Statistics")


    sql = """
    SELECT id, time
    FROM exectimes
    """
    exec_times_dict = get_all_times(db_path,sql)
    rank,max_exec_time = max_value_in_dict(exec_times_dict)
    if max_exec_time == None or rank == None:
         print("Error occured in max exec time")
    else:
         print(f"Maximum Execution time: {max_exec_time:.3f} s, Rank: {rank}")

    avg_exec = avg_value_in_dict(exec_times_dict)
    if avg_exec != None:
        print(f"Average Execution time across {size} MPI Ranks: {avg_exec:.3f} s")

    sql = """
    SELECT rank, total_time
    FROM mpi_time_sum
    """
    mpi_times_dict = get_all_times(db_path,sql)
    rank,max_mpi_time = max_value_in_dict(mpi_times_dict)
    if max_mpi_time == None or rank == None:
         print("Error occured in max MPI time")
    else:
         print(f"Maximum MPI time: {max_mpi_time:.3f} s, Rank: {rank}")

    avg_mpi = avg_value_in_dict(mpi_times_dict)
    if avg_exec != None and avg_mpi != None:
        print(f"Average MPI time across {size} MPI Ranks: {avg_mpi:.3f} s")
        print(f"Average Ratio of MPI time to Execution time across {size} MPI Ranks: {(avg_mpi/avg_exec)*100:.2f}%")
    ratios = dict_ratios(mpi_times_dict,exec_times_dict)
    rank,max_ratio = max_value_in_dict(ratios)
    print(f"Maximum Ratio of MPI time to Execution time: {max_ratio:.2f}%, Rank: {rank}\n")
    print_decoration(RESET)


def output_to_csv(plot_data,csv_file):
    # Prepare the data for writing to a CSV file
    rows = []
    for operation, comm_data in plot_data.items():
        for comm, avg_time in comm_data.items():
            rows.append([operation, comm, avg_time])

    # Write the data to a CSV file
    with open(csv_file,'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Call', 'Communicator', 'Time'])
        writer.writerows(rows)



def get_average_time_per_operation_top(db_path, n):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # SQL query to group by MPI operation, buffer size range, and calculate average time
    sql = """
    SELECT o.operation, d.buffer_size_min, d.buffer_size_max, c.name, c.size, AVG(d.time) as avg_time
    FROM data d
    JOIN comms c ON d.comm_id = c.id
    JOIN operations o ON d.operation_id = o.id
    GROUP BY o.operation, d.buffer_size_min, d.buffer_size_max, c.name
    HAVING AVG(d.time) > 0.001
    ORDER BY avg_time DESC
    """
    plot_data = {}
    try:
        cursor.execute(sql)
        rows = cursor.fetchall()
        # Aggregate data into a structure suitable for plotting
        for row in rows:
            operation, buf_min, buf_max, comm_name, comm_size, avg_time = row
            op_with_buf = f"{operation} ({buf_min}-{buf_max})"
            name_with_size = f"{comm_name} ({comm_size})"
            if op_with_buf not in plot_data:
                plot_data[op_with_buf] = {}
            plot_data[op_with_buf][name_with_size] = avg_time

        # Get the top N MPI operations by the total average time
        sorted_ops = sorted(plot_data.items(), key=lambda item: sum(item[1].values()), reverse=True)
        if ( n > 0 ):
            top_ops_data = dict(sorted_ops[:n])
        else:
           top_ops_data = dict(sorted_ops)
        return top_ops_data

    except sqlite3.Error as e:
        print("An error occurred:", e)
    finally:
        conn.close()

def get_average_time_per_communicator_top(db_path,n):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # SQL query to group by communicator, MPI operation, and buffer size range
    # and calculate average time
    sql = """
    SELECT c.name, AVG(d.time) as avg_time
    FROM data d
    JOIN comms c ON d.comm_id = c.id
    JOIN operations o ON d.operation_id = o.id
    GROUP BY c.name, o.operation, d.buffer_size_min, d.buffer_size_max
    ORDER BY avg_time DESC
    """

    communicator_list = []
    try:
        cursor.execute(sql)
        rows = cursor.fetchall()
        communicator_totals = {}
        for row in rows:
            comm_name, avg_time = row

            if comm_name in communicator_totals:
                communicator_totals[comm_name] += avg_time
            else:
                communicator_totals[comm_name] = avg_time

        communicator_list = sorted(communicator_totals.items(), key=lambda x: x[1], reverse=True)

    except sqlite3.Error as e:
        print("An error occurred:", e)
        return []
    finally:
        conn.close()

    return communicator_list[:n]
        




def get_all_comms(db):
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    try:
        cursor.execute("""
        SELECT c.name
        FROM comms c
        """)

        result = cursor.fetchall()
        if not result:
            print("No data found.")
            return
        for item in result:
            print(item)

    except sqlite3.Error as e:
        print("An error occurred:", e)
    finally:
        conn.close()


def print_comms_ranks(db):
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    try:
        cursor.execute("""
        SELECT c.name, d.rank
        FROM data d
        JOIN comms c ON d.comm_id = c.id
        GROUP BY c.name, d.rank
        """)

        result = cursor.fetchall()
        if not result:
            print("No data found.")
            return
        for item in result:
            print(item)

    except sqlite3.Error as e:
        print("An error occurred:", e)
    finally:
        conn.close()

def get_ranks_by_comm(db,comm):
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    ranks = []
    try:
        cursor.execute("""
        SELECT DISTINCT d.rank
        FROM data d
        JOIN comms c ON d.comm_id = c.id
        WHERE c.id = ?;
        """, (comm,))
        ranks = [item[0] for item in cursor.fetchall()]


    except sqlite3.Error as e:
        print("An error occurred:", e)
    finally:
        conn.close()
        return ranks

# Return a new truncated list of 4 elements: [1,2,3,4,5] -> [1,2,...,5]
def list_truncate_tostr(l):
    if len(l) <= 4:
        return str(l).replace(" ", "")
    return f"[{l[0]},{l[1]},...,{l[-1]}]"

# Print the volume table
def print_volume_table(db_path):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute("""
    SELECT
    v.comm_id AS cid,
    c.name AS comm_name,
    v.rank,
    v.operation_id,
    o.operation,
    v.volume
FROM ops_volume v
JOIN comms c ON v.comm_id = c.id
JOIN operations o ON v.operation_id = o.id
    """)


    rows = cursor.fetchall()
    print_decoration(BOLD)
    print(f"{'Comm Name':<15}{'Operation':<25}{'Rank':<10}{'Volume':<20}")
    print_decoration(RESET)
    for row in rows:
        comm_id, comm_name, rank, operation_id, operation, volume = row
        print(f"{comm_name:<15}{operation:<25}{rank:<10}{volume:<20}")

    conn.close()


# Dump the volume table
def dump_volume_table(db_path):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute("SELECT * FROM ops_volume")


    rows = cursor.fetchall()
    print_decoration(BOLD)
    print(f"{'ID':<5}{'Operation ID':<15}{'Rank':<5}{'Comm ID':<10}{'Volume':<20}")
    print_decoration(RESET)
    for row in rows:
        id, operation_id, rank, comm_id, volume = row
        print(f"{id:<5}{operation_id:<15}{rank:<5}{comm_id:<10}{volume:<20}")

    conn.close()

def summarize_volume_by_comm_operation(db_path):
    sql = """
    SELECT
    v.comm_id AS cid,
    v.operation_id AS oid,
    c.name AS comm_name,
    o.operation,
    SUM(v.volume) AS total_volume
    FROM ops_volume v
    JOIN comms c ON v.comm_id = c.id
    JOIN operations o ON v.operation_id = o.id
    GROUP BY c.name, o.operation
          """

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute(sql)
    rows = cursor.fetchall()

    volume_summary = {}
    for row in rows:
        cid,oid,comm_name, operation, total_volume = row
        volume_summary[(oid,cid)] = total_volume

    conn.close()
    return volume_summary




def main():
    parser = argparse.ArgumentParser(description="mpisee SQLite database basic query tool. **Requires input database file.**  Summarizes MPI data across ranks by default. Use filters to refine the output.")
    parser.add_argument("-i", "--inputdb", required=True, help="REQUIRED: path to the input SQLite database file of mpisee.")
    display_group = parser.add_argument_group("Main Output Options: Summarizes MPI data across ranks by default if none of the below options are selected")
    display_group.add_argument("-a", "--all",  action='store_true', required=False, help="Displays data among all ranks.")
    display_group.add_argument("-p", "--pt2pt",action='store_true', required=False, help="Displays only point to point MPI operations,")
    display_group.add_argument("-c", "--collectives", action='store_true', required=False, help="Displays only collective MPI operations.")
    filter_group = parser.add_argument_group("Filtering Options: Combined with the main output options to refine the output")
    filter_group.add_argument("-r", "--ranks", type=str, required=False, help="Select the data of specific MPI ranks.  Use a comma-separated list of ranks. Cannot be combined with the default query.")
    filter_group.add_argument("-b", "--buffsize", type=str, required=False, help="Select the data for a specific buffer size range defined as min:max.")
    filter_group.add_argument("-t", "--time", type=str, required=False, help="Select the data for a specific time range in seconds defined as min:max.")
    filter_group.add_argument("-s", "--sort", required=False,  type=int, default=1, help="Sort the results using one of the following integers: 0 by communicator, 1 descending by time(default), 2 ascending by time, 3 by MPI operation, 4 ascending by buffer size, 5 descending by buffer size, 6 ascending by number of calls, 7 descending by number of calls.")
    output_group = parser.add_argument_group("Output Options")
    output_group.add_argument("-o","--outputcsv", required=False, type=str, help="Output to a csv file")
    second_display_group = parser.add_argument_group("Secondary Output Options")
    second_display_group.add_argument("-e", "--exectime", required=False, action='store_true', help="Display the net MPI time and the total Execution time for each process.")
    second_display_group.add_argument("-u", "--ctime", action='store_true', required=False, help="Display the time summary for each communicator.")
    parser.add_argument("--debug", required=False, action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()

    #header_path = '../utils.h'

    # Path to the script file (this script)
    script_path = os.path.abspath(__file__)

    # Directory where the script is located
    script_dir = os.path.dirname(script_path)

    # Path to the header file, relative to the script location
    header_path = os.path.join(script_dir, '../utils.h')

    # Normalize the path to resolve any ".." components
    header_path = os.path.normpath(header_path)
    enum_primitives = parse_enum_from_header(header_path)

    db_path = args.inputdb

    if args.ranks:
        rank_list = [int(rank) for rank in args.ranks.split(',')]
    else:
        rank_list = []

    if args.buffsize:
        tmp = args.buffsize.split(':')[0]
        if ( len(tmp) > 0 ):
            buffsizemin = int(args.buffsize.split(':')[0])
        else:
            buffsizemin = 0
        tmp = args.buffsize.split(':')[1]
        if ( len(tmp) > 0 ):
            buffsizemax = int(args.buffsize.split(':')[1])
        else:
            buffsizemax = 2147483647
    else:
        buffsizemax = 2147483647
        buffsizemin = 0

    if args.time:
        tmp = args.time.split(':')[0]
        if ( len(tmp) > 0 ):
            timemin = float(args.time.split(':')[0])
        else:
            timemin = -1.0
        tmp = args.time.split(':')[1]
        if ( len(tmp) > 0 ):
            timemax = float(args.time.split(':')[1])
        else:
            timemax = sys.float_info.max
    else:
        timemax = sys.float_info.max
        timemin = -1.0



    create_and_populate_summary_table(db_path)

    print_metadata_table(db_path)

    print_general_stats(db_path)

    print_decoration(RESET)

    if args.pt2pt:
        if ( args.ranks ):
            query_ranks(db_path,enum_primitives,'Ibsend',buffsizemin,buffsizemax,timemin,timemax,rank_list,args.sort,args.outputcsv)
        else:
            query_colls_pt2pt(db_path,enum_primitives,'Ibsend',buffsizemin,buffsizemax,timemin,timemax,args.sort,args.outputcsv)
    elif args.collectives:
        if ( args.ranks ):
            query_ranks(db_path,enum_primitives,'Bcast',buffsizemin,buffsizemax,timemin,timemax,rank_list,args.sort,args.outputcsv)
        else:
            query_colls_pt2pt(db_path,enum_primitives,'Bcast',buffsizemin,buffsizemax,timemin,timemax,args.sort,args.outputcsv)
    elif args.exectime:
        print_execution_time(db_path,rank_list)
    elif args.all:
        query_ranks(db_path,enum_primitives,None,buffsizemin,buffsizemax,timemin,timemax,rank_list,args.sort,args.outputcsv)
    elif args.debug:
        print_comms_table(db_path)
        print_operations_table(db_path)
        print_data_table(db_path)
        print_volume_table(db_path)
        summarize_volume_by_comm_operation(db_path)
    elif args.ctime:
        query_summarize_time(db_path,args.sort,args.outputcsv)
    else:
        query_colls_pt2pt(db_path,enum_primitives,None,buffsizemin,buffsizemax,timemin,timemax,args.sort,args.outputcsv)


if __name__ == "__main__":
    main()

