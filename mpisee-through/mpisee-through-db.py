#!/usr/bin/env python3
import argparse
import sqlite3
import re
import os
import sys

def parse_enum_from_header(header_path):
    with open(header_path, 'r') as file:
        content = file.read()

    # Regular expression to find enum definition
    enum_pattern = re.compile(r'enum primitives\{([^}]+)\};', re.MULTILINE | re.DOTALL)
    match = enum_pattern.search(content)

    if match:
        enum_content = match.group(1)
        # Extract individual enum items
        enum_items = enum_content.split(',')
        enum_dict = {}
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
        return enum_dict
    else:
        print("Enum 'primitives' not found in the header file.")
        return None

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

def exec_query_and_print(db_path,sql,order,num_of_rows,ranks,*args):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    if len(ranks) > 0:
            placeholders = ','.join('?' * len(ranks))
            sql += f" AND d.rank IN ({placeholders})"
            params = args + tuple(ranks)
    else:
        params = args
    sql = select_order(sql,order)
    try:
        # Execute the query
        cursor.execute(sql,params)

        # Print header
        print(f"{'Comm Name':<15}{'Comm Size':<15}{'Rank':<10}{'MPI Operation':<20}"
              f"{'Buffer Size Range':<20}{'Calls':<15}{'Time (s)':<15}{'% of MPI Time':<20}{'% of Total Time':<10}")

        # Print rows
        r = 0
        prev_rank = -1
        percentage_mpi_time = 0.0
        percentage_exec_time = 0.0
        for row in cursor.fetchall():
            name, size, rank, operation, buf_min, buf_max, calls, time = row
            buffer_size = f"{buf_min} - {buf_max}"
            if rank != prev_rank:
                percentage_exec_time = (time/get_exec_time_by_rank(db_path,rank))*100
                percentage_mpi_time = (time/get_mpi_time_by_rank(db_path,rank))*100
                prev_rank = rank

            print(f"{name:<15}{size:<15}{rank:<10}{operation:<20}"
                  f"{buffer_size:<20}{calls:<15}{time:<15.3f}{percentage_mpi_time:<20.3f}{percentage_exec_time:<10.3f}")
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
    ORDER BY c.name
    """

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

def print_execution_time(dpath,order=1,ranks=[]):
    conn = sqlite3.connect(dpath)
    cursor = conn.cursor()

    sql = """
    SELECT t.id, t.time
    FROM exectimes t
    """
    params = ""
    if len(ranks) > 0:
            placeholders = ','.join('?' * len(ranks))
            sql += f"WHERE t.id IN ({placeholders})"
            params = tuple(ranks)

    elif order == 1:
        sql += """
        ORDER BY t.time DESC"""
    elif order == 2:
        sql += """
        ORDER BY t.time ASC"""


    try:
        # Execute the query
        cursor.execute(sql,(params))

        # Print header
        print(f"{'MPI Rank':<10}{'Execution Time (s)':<10}")

        # Print rows
        for row in cursor.fetchall():
            id,time = row
            print(f"{id:<10}{time:<10.4f}")

    except sqlite3.Error as e:
        print("Failed to read data from SQLite exectime table", e)
    finally:
        # Close the database connection
        if conn:
            conn.close()


def print_data_by_time(dbpath,order=1,num_of_rows=0,rank_list=[],*args):
    # SQL query
    sql = """
    SELECT c.name, c.size, d.rank, o.operation, d.buffer_size_min, d.buffer_size_max,
           d.calls, d.time
    FROM data d
    JOIN comms c ON d.comm_id = c.id
    JOIN operations o ON d.operation_id = o.id
    WHERE d.time >= ? AND d.time <= ?
    """
    exec_query_and_print(dbpath,sql,order,num_of_rows,rank_list,*args)

def print_data_by_bufsize(dbpath,order=1,num_of_rows=0,rank_list=[],*args):
    # SQL query
    sql = """
    SELECT c.name, c.size, d.rank, o.operation, d.buffer_size_min, d.buffer_size_max,
           d.calls, d.time
    FROM data d
    JOIN comms c ON d.comm_id = c.id
    JOIN operations o ON d.operation_id = o.id
    WHERE d.buffer_size_min >= ? AND d.buffer_size_max <= ?
    """
    exec_query_and_print(dbpath,sql,order,num_of_rows,rank_list,*args)

def  print_data_collectives(dbpath,order=1,num_of_rows=0,rank_list=[],*args):
    sql = """
        SELECT c.name, c.size, d.rank, o.operation, d.buffer_size_min, d.buffer_size_max,
           d.calls, d.time
        FROM data d
        JOIN comms c ON d.comm_id = c.id
        JOIN operations o ON d.operation_id = o.id
        WHERE d.buffer_size_min >= ? AND d.buffer_size_max <= ? AND d.operation_id >= ?  """
    exec_query_and_print(dbpath,sql,order,num_of_rows,rank_list,*args)


def print_data_pt2pt(dbpath,order=1,num_of_rows=0,rank_list=[],*args):
    sql = """
        SELECT c.name, c.size, d.rank, o.operation, d.buffer_size_min, d.buffer_size_max,
           d.calls, d.time
        FROM data d
        JOIN comms c ON d.comm_id = c.id
        JOIN operations o ON d.operation_id = o.id
        WHERE d.buffer_size_min >= ? AND d.buffer_size_max <= ? AND d.operation_id <= ?  """
    exec_query_and_print(dbpath,sql,order,num_of_rows,rank_list,*args)

def query_all_data(dbpath,order=1,num_of_rows=0,rank_list=[],*args):
    sql = """SELECT c.name, c.size, d.rank, o.operation,
                      d.buffer_size_min, d.buffer_size_max, d.calls, d.time
                      FROM data d
                      JOIN comms c ON d.comm_id = c.id
                      JOIN operations o ON d.operation_id = o.id """
    exec_query_and_print(dbpath,sql,order,num_of_rows,rank_list,*args)

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

def mpi_time(dbpath,order=1,ranks=[]):
    conn = sqlite3.connect(dbpath)
    cursor = conn.cursor()
    sql ="""
    SELECT rank, total_time as mpi_time
    FROM mpi_time_sum
    """

    if len(ranks) > 0:
        placeholders = ','.join('?' * len(ranks))
        sql += f" WHERE rank IN ({placeholders})"

    sql += f"GROUP BY rank "

    if order == 1:
        sql += f" ORDER BY mpi_time DESC"
    else:
        sql += f" ORDER BY mpi_time ASC"

    try:
        cursor.execute(sql)
        rows = cursor.fetchall()
        print(f"{'Rank':<8}{'MPI Time':<15}")
        for row in rows:
            rank, total_time = row
            print(f"{rank:<10}{total_time:.3f}")
    except sqlite3.Error as e:
        print("An error occurred:", e)
    finally:
        conn.close()

def print_metadata_table(db_path):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute("SELECT key, value FROM metadata")

    rows = cursor.fetchall()
    for row in rows:
        key, value = row
        print(f"{key}: {value}")

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

def print_general_stats(db_path):
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

    print("Overall statistics")


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

    avg_exec = avg_value_in_dict(exec_times_dict)
    avg_mpi = avg_value_in_dict(mpi_times_dict)
    print(f"Average Execution time across {size} processes: {avg_exec:.3f} s")
    print(f"Average MPI time across {size} processes: {avg_mpi:.3f} s")
    print(f"MPI Time to Execution time Ratio: {(avg_mpi/avg_exec)*100:.2f}%")

    # sql = """
    # SELECT rank, total_time
    # FROM mpi_time_sum
    # WHERE total_time = (SELECT MAX(total_time) FROM mpi_time_sum)
    # """
    # rank, time_max = get_max_time_rank(db_path,sql)
    # if rank < 0 or time_max < 0.0:
    #      print("Error occured in max mpi time")
    # else:
    #      print(f"Maximum MPI time: {time_max:.3f} s, Rank: {rank}")

    # sql = """
    # SELECT AVG(total_time)
    # FROM mpi_time_sum
    # """
    # avg_mpi_time = get_avg_time(db_path,sql)
    # if avg_mpi_time < 0.0:
    #      print("Error occured in average mpi time calculation")
    # else:
    #      print(f"Average MPI time across {size} processes: {avg_mpi_time:.3f} s")
    # sql = """
    # SELECT id, time
    # FROM exectimes
    # WHERE time = (SELECT MAX(time) FROM exectimes)
    # """
    # rank, time_max = get_max_time_rank(db_path,sql)
    # if rank < 0 or time_max < 0.0:
    #      print("Error occured in max exec time time")
    # else:
    #     print(f"Maximum Execution time: {time_max:.3f} s, Rank: {rank}")

    # sql = """
    # SELECT AVG(time)
    # FROM exectimes
    # """
    # avg_exec_time = get_avg_time(db_path,sql)
    # if avg_exec_time < 0.0:
    #      print("Error occured in average exec time calculation")
    # else:
    #      print(f"Average Execution time across {size} processes: {avg_exec_time:.3f} s")


    print()


def main():
    parser = argparse.ArgumentParser(description="Query the mpisee SQLite database.")
    parser.add_argument("-d", "--db_path", required=True, help="Path to the mpisee SQLite database file.")
    parser.add_argument("-e", "--exectime", required=False, action='store_true', help="Print the execution time for each process.")
    parser.add_argument("-p", "--pt2pt",action='store_true', required=False, help="Show only point to point MPI operations,")
    parser.add_argument("-c", "--collectives", action='store_true', required=False, help="Show only collective MPI operations.")
    parser.add_argument("-r", "--ranks", type=str, required=False, help="Show the data of specific MPI ranks.")
    parser.add_argument("-b", "--buffsize", type=str, required=False, help="Show the data for a specific buffer size range defined as min:max.")
    parser.add_argument("-t", "--time", type=str, required=False, help="Show the data for a specific time range in seconds defined as min:max.")
    parser.add_argument("-m", "--mpitime", action='store_true', required=False, help="Show the data for a specific time range in seconds defined as min:max.")
    parser.add_argument("-n", "--nresults", required=False, type=int, default=0, help="Show the first N results. By default all are printed.")
    parser.add_argument("-s", "--sort", required=False,  type=int, default=1, help="Sort the results: 0 by communicator, 1 descending by time(default), 2 ascending by time, 3 by MPI operation, 4 ascending by buffer size, 5 descending by buffer size, 6 ascending by number of calls, 7 descending by number of calls.")
    args = parser.parse_args()

    header_path = '../utils.h'

    # Path to the script file (this script)
    script_path = os.path.abspath(__file__)

    # Directory where the script is located
    script_dir = os.path.dirname(script_path)

    # Path to the header file, relative to the script location
    header_path = os.path.join(script_dir, '../utils.h')

    # Normalize the path to resolve any ".." components
    header_path = os.path.normpath(header_path)
    enum_primitives = parse_enum_from_header(header_path)

    db_path = args.db_path

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
            timemin = 0
        tmp = args.time.split(':')[1]
        if ( len(tmp) > 0 ):
            timemax = float(args.time.split(':')[1])
        else:
            timemax = sys.float_info.max
    else:
        timemax = -1
        timemin = sys.float_info.max

    create_and_populate_summary_table(db_path)

    print_metadata_table(db_path)

    print_general_stats(db_path)

    if args.pt2pt:
        print_data_pt2pt(db_path,args.sort,args.nresults,rank_list,buffsizemin,buffsizemax,enum_primitives['Issend'])
    elif args.collectives:
        print_data_collectives(db_path,args.sort,args.nresults,rank_list,buffsizemin,buffsizemax,enum_primitives['Bcast'])
    elif args.buffsize:
        print_data_by_bufsize(db_path,args.sort,args.nresults,rank_list,buffsizemin,buffsizemax)
    elif args.time:
        print_data_by_time(db_path,args.sort,args.nresults,rank_list,timemin,timemax)
    elif args.exectime:
        print_execution_time(db_path,args.sort,rank_list)
    elif args.mpitime:
        mpi_time(db_path,args.sort,rank_list)
    else:
        query_all_data(db_path,args.sort,args.nresults,rank_list)

if __name__ == "__main__":
    main()

