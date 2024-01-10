#!/usr/bin/env python3
import argparse
import sqlite3
import re

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


def query_all_data(db_path):
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT c.name, c.size, d.rank, o.operation, "
                      "d.buffer_size_min, d.buffer_size_max, d.calls, d.time "
                      "FROM data d "
                      "JOIN comms c ON d.comm_id = c.id "
                      "JOIN operations o ON d.operation_id = o.id "
                       "ORDER BY c.name;")
        return cursor.fetchall()

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
    #print(sql)
    try:
        # Execute the query
        cursor.execute(sql,params)

        # Print header
        print(f"{'Comm Name':<15}{'Comm Size':<15}{'Rank':<10}{'Operation':<20}"
              f"{'Buffer Size Range':<25}{'Calls':<15}{'Time':<20}")

        # Print rows
        r = 0
        for row in cursor.fetchall():
            name, size, rank, operation, buf_min, buf_max, calls, time = row
            buffer_size = f"{buf_min} - {buf_max}"
            print(f"{name:<15}{size:<15}{rank:<10}{operation:<20}"
                  f"{buffer_size:<25}{calls:<15}{time:<20}")
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


def print_data_by_bufsize(dbpath, min=0, max=2147483647):
    # Connect to the SQLite database
    conn = sqlite3.connect(dbpath)
    cursor = conn.cursor()

    # SQL query
    sql = """
    SELECT c.name, c.size, d.rank, o.operation, d.buffer_size_min, d.buffer_size_max,
           d.calls, d.time
    FROM data d
    JOIN comms c ON d.comm_id = c.id
    JOIN operations o ON d.operation_id = o.id
    WHERE d.buffer_size_min >= ? AND d.buffer_size_max <= ?
    ORDER BY c.name
    """

    try:
        # Execute the query
        cursor.execute(sql,(min,max))

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

def print_data_by_time(dbpath, min=0.0, max=0.0):
    # Connect to the SQLite database
    conn = sqlite3.connect(dbpath)
    cursor = conn.cursor()

    if max == 0.0:
        cursor.execute("SELECT MAX(time) FROM data")
        max = cursor.fetchone()[0]

    # SQL query
    sql = """
    SELECT c.name, c.size, d.rank, o.operation, d.buffer_size_min, d.buffer_size_max,
           d.calls, d.time
    FROM data d
    JOIN comms c ON d.comm_id = c.id
    JOIN operations o ON d.operation_id = o.id
    WHERE d.time >= ? AND d.time <= ?
    ORDER BY c.name
    """

    try:
        # Execute the query
        cursor.execute(sql,(min,max))

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

def print_data_by_time_sort_desc(dbpath, min=0.0, max=0.0):
    # Connect to the SQLite database
    conn = sqlite3.connect(dbpath)
    cursor = conn.cursor()

    if max == 0.0:
        cursor.execute("SELECT MAX(time) FROM data")
        max = cursor.fetchone()[0]

    # SQL query
    sql = """
    SELECT c.name, c.size, d.rank, o.operation, d.buffer_size_min, d.buffer_size_max,
           d.calls, d.time
    FROM data d
    JOIN comms c ON d.comm_id = c.id
    JOIN operations o ON d.operation_id = o.id
    WHERE d.time >= ? AND d.time <= ?
    ORDER BY c.name
    """

    try:
        # Execute the query
        cursor.execute(sql,(min,max))

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

def main():
    parser = argparse.ArgumentParser(description="Query the mpisee SQLite database.")
    parser.add_argument("-d", "--db_path", required=True, help="Path to the mpisee SQLite database file.")
    parser.add_argument("-p","--pt2pt",action='store_true', required=False, help="Show only point to point MPI operations,")
    parser.add_argument("-c","--collectives", action='store_true', required=False, help="Show only collective MPI operations.")
    parser.add_argument("-r","--ranks", type=str, required=False, help="Show the data of specific MPI ranks.")
    parser.add_argument("-b","--buffsize", type=str, required=False, help="Show the data for a specific buffer size range defined as min:max.")
    parser.add_argument("-n","--nresults", required=False, type=int, default=0, help="Show the first N results. By default all are printed.")
    parser.add_argument("-s","--sort", required=False,  type=int, default=1, help="Sort the results: 0 by communicator, 1 descending by time(default), 2 ascending by time, 3 by MPI operation, 4 ascending by buffer size, 5 descending by buffer size, 6 ascending by number of calls, 7 descending by number of calls.")
    args = parser.parse_args()

    header_path = '../utils.h'
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


    #print_data_by_rank(db_path,0)

    #print_data_by_comm(db_path,"W")

    #print_data_by_bufsize(db_path, min=128)

    #print_data_by_time(db_path, min=0.01)

    #print_data_by_time_sort_desc(db_path, min=0.01)

    if args.pt2pt:
        print_data_pt2pt(db_path,args.sort,args.nresults,rank_list,buffsizemin,buffsizemax,enum_primitives['Issend'])
    elif args.collectives:
        print_data_collectives(db_path,args.sort,args.nresults,rank_list,buffsizemin,buffsizemax,enum_primitives['Bcast'])
    elif args.buffsize:
        print_data_by_bufsize(db_path, max=buffsizemax, min=buffsizemin)
    else:
        print_all_data(db_path)

if __name__ == "__main__":
    main()

