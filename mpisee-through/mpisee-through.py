#!/usr/bin/env python3
"""
Module documentation.
"""

# Imports
import sys
import csv
from operator import add
import argparse
#from art import *

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

# Class declarations

# Function declarations

def print_decoration(decoration):
    if sys.stdout.isatty():
        # the output is not redirected, we can use a fancy style:
        sys.stdout.write(decoration)


def ratio_to_percentage(ratio):
    return f"{ratio * 100:.2f}%"

def get_values(inlist):
    ranks = []
    values = []

    # Flatten the list if its first element is a list
    if isinstance(inlist[0], list):
        inlist = inlist[0] + inlist[1:]

    for entry in inlist[1:]:  # skipping the first informational string
        rank, value = entry.split()
        ranks.append(int(rank))
        values.append(float(value))

    max_value = max(values)
    min_value = min(values)

    max_rank = ranks[values.index(max_value)]
    min_rank = ranks[values.index(min_value)]

    return max_rank,max_value,min_rank,min_value


def print_times(elapsed_time, mpi_times):
    max_rank, max_value, min_rank, min_value = get_values(elapsed_time)
    max_mpi_rank, max_mpi_time, min_mpi_rank, min_mpi_time = get_values(mpi_times)

    ranks = []
    ratios = []

    mpi_times = mpi_times[0]+mpi_times[1:]
    for i in range(1,len(elapsed_time)):  # skipping the first informational string
        time = elapsed_time[i].split()[1]
        mpi_time = mpi_times[i].split()[1]
        ranks.append(int(i))
        if time == 0:
            return "Cannot compute ratio: Division by zero"
        ratios.append(float(mpi_time)/float(time))
    
    max_ratio = max(ratios)
    min_ratio = min(ratios)
    max_ratio_rank = ranks[ratios.index(max_ratio)]
    min_ratio_rank = ranks[ratios.index(min_ratio)]

    print_decoration(GREEN)
    print(f"Overall Timing Statistics for {len(ratios)} MPI Processes (Ranks in MPI_COMM_WORLD)")
    print_decoration(RESET)
    print(f"Maximum Total Time: {max_value}s (MPI Rank: {max_rank})")
    print(f"Minimum Total Time: {min_value}s (MPI Rank: {min_rank})")
    print(f"Maximum MPI Time: {max_mpi_time}s (MPI Rank: {max_mpi_rank})")
    print(f"Minimum MPI Time: {min_mpi_time}s (MPI Rank: {min_mpi_rank})")
    print(f"Maximum Percentage of MPI Time to Total Time: {ratio_to_percentage(max_ratio)} (MPI Rank: {max_ratio_rank})")
    print(f"Minimum Percentage of MPI Time to Total Time: {ratio_to_percentage(min_ratio)} (MPI Rank: {min_ratio_rank})")
    print("\n")



def compact_proc_list(proc_list):
    proc_list.sort()
    proc_str = ""
    start_proc = proc_list[0]
    for i in range(1, len(proc_list)):
        if proc_list[i - 1] + 1 != proc_list[i]:
            proc_str += (
                f"{start_proc}-{proc_list[i-1]}"
                if start_proc != proc_list[i - 1]
                else f"{start_proc}"
            )
            proc_str += ", "
            start_proc = proc_list[i]

    proc_str += (
        f"{start_proc}-{proc_list[-1]}"
        if start_proc != proc_list[-1]
        else f"{start_proc}"
    )

    return proc_str


def print_mapping(mapping):
    mapping[0] = mapping[0].split(":")[1].strip()
    node_to_procs = {}
    for proc in mapping:
        proc_rank, node = proc.split()
        proc_rank = int(proc_rank)
        if node in node_to_procs:
            node_to_procs[node].append(proc_rank)
        else:
            node_to_procs[node] = [proc_rank]

    print_decoration(GREEN)
    print("Mapping of MPI ranks to Compute Nodes")
    print_decoration(RESET)
    for node, proc_list in node_to_procs.items():

        print(f"{node}: {compact_proc_list(proc_list)}")

    print(end="\n")


def prepare_data(file_path):
    with open(file_path, "r") as file:
        csv_file = csv.reader(file)
        for i in range(7): #skip those lines for now
            next(csv_file)
        mapping = next(csv_file)
        time_elapsed = next(csv_file)                    #skip also this line
        index_to_colname = next(csv_file)
        colname_to_index = {
            index_to_colname[i]: i for i in range(0, len(index_to_colname))
        }
        raw_data = []
        mpi_data = []
        for row in csv_file:
            row_parsed = []
            if '#' in row[0]:
                mpi_data.append(row)
            else:
                for j in range(0, len(row)):
                    row_parsed.append(
                        float(row[j]) if colname_to_index["Comm"] != j else row[j]
                    )
                raw_data.append(row_parsed)

    data_groupBy_comm, comm_to_procs = groupBy_comm(raw_data, colname_to_index)

    table = create_table(data_groupBy_comm, colname_to_index, index_to_colname)

    return table, mapping, comm_to_procs, time_elapsed, mpi_data


def groupBy_comm(data, colname_to_index):
    data_groupBy_comm = {}
    comm_to_procs = {}
    for row in data:
        comm = row[colname_to_index["Comm"]]
        row[colname_to_index["Comm"]] = ""
        if comm in data_groupBy_comm:
            data_groupBy_comm[comm].append(row)
            comm_to_procs[comm].append(int(row[colname_to_index["Rank"]]))
        else:
            data_groupBy_comm[comm] = [row]
            comm_to_procs[comm] = [int(row[colname_to_index["Rank"]])]

    return data_groupBy_comm, comm_to_procs


def create_table(data_groupBy_comm, colname_to_index, index_to_colname):
    table = []

    time_indexes = [kv[1] for kv in colname_to_index.items() if kv[0].endswith("Time")]

    for comm in data_groupBy_comm:
        maximum = list(data_groupBy_comm[comm][0])
        minimum = list(data_groupBy_comm[comm][0])
        cum_sum = list(data_groupBy_comm[comm][0])

        for i in range(1, len(data_groupBy_comm[comm])):
            row = data_groupBy_comm[comm][i]
            maximum = list(map(max, maximum, row))
            minimum = list(map(min, minimum, row))
            cum_sum = list(map(add, cum_sum, row))

        average = list(
            map(
                lambda x: x / len(data_groupBy_comm[comm])
                if not isinstance(x, str)
                else x,
                cum_sum,
            )
        )

        comm_size = int(data_groupBy_comm[comm][0][colname_to_index["Size"]])
        for i in time_indexes:
            call = "_".join(index_to_colname[i].split("_")[:-1])
            call_mean, call_min, call_max = average[i], minimum[i], maximum[i]
            call_tot_vol = int(cum_sum[colname_to_index[call + "_Volume"]])
            nb_tot_call = int(cum_sum[colname_to_index[call + "_Calls"]])

            if nb_tot_call > 0:
                table.append(
                    {
                        "call": call,
                        "call_mean": call_mean,
                        "call_min": call_min,
                        "call_max": call_max,
                        "call_tot_vol": call_tot_vol,
                        "nb_tot_call": nb_tot_call,
                        "comm": comm,
                    }
                )

    return table


def print_cct(table, comm_to_procs, comm_limit):
    comm_limit = (
        comm_limit
        if comm_limit != None and comm_limit < len(comm_to_procs)
        else len(comm_to_procs)
    )
    nb_comm_printed = 0
    print_decoration(GREEN)
    print("Statistics Per Communicator")
    print_decoration(RESET)
    for comm in comm_to_procs.keys():
        if nb_comm_printed >= comm_limit:
            break

        comm_calls = list(filter(lambda x: x["comm"] == comm, table))
        comm_calls.sort(key=lambda x: x["call_mean"], reverse=True)

        print_decoration(BOLD)
        print("COMM".ljust(20) + "SIZE".ljust(10) + "PROCS")
        print_decoration(RESET)
        print(
            f"{comm}".ljust(20)
            + f"{len(comm_to_procs[comm])}".ljust(10)
            + f"{compact_proc_list(comm_to_procs[comm])}",
            end="\n\n",
        )
        print_decoration(RESET)
        print(
            "\t"
            + "Call".ljust(20)
            + "Mean[s]".rjust(10)
            + "Min[s]".rjust(10)
            + "Max[s]".rjust(10)
            + "Volume".rjust(15)
            + "#Calls".rjust(10)
        )

        print("\t" + "".join(["-"] * 75))

        for i in range(0, len(comm_calls)):
            call = comm_calls[i]["call"]
            call_mean = comm_calls[i]["call_mean"]
            call_min = comm_calls[i]["call_min"]
            call_max = comm_calls[i]["call_max"]
            call_tot_vol = comm_calls[i]["call_tot_vol"]
            nb_tot_call = comm_calls[i]["nb_tot_call"]
            print(
                "\t"
                + f"{call}".ljust(20)
                + f"{call_mean:.4f}".rjust(10)
                + f"{call_min:.4f}".rjust(10)
                + f"{call_max:.4f}".rjust(10)
                + f"{call_tot_vol}".rjust(15)
                + f"{nb_tot_call}".rjust(10)
            )

        print("\n")
        nb_comm_printed += 1

    return


def print_tc(table, comm_to_procs, limit):

    table.sort(key=lambda x: x["call_mean"], reverse=True)

    print("".join(["-"] * 100))
    print(
        "Call".ljust(20)
        + "Mean[s]".rjust(10)
        + "Min[s]".rjust(10)
        + "Max[s]".rjust(10)
        + "Volume".rjust(15)
        + "#Calls".rjust(10)
        + "Comm".rjust(15)
        + "Size".rjust(10)
    )

    print("".join(["-"] * 100))

    end_range = limit if limit != None and limit < len(table) else len(table)
    for i in range(0, end_range):
        call = table[i]["call"]
        call_mean = table[i]["call_mean"]
        call_min = table[i]["call_min"]
        call_max = table[i]["call_max"]
        call_tot_vol = table[i]["call_tot_vol"]
        nb_tot_call = table[i]["nb_tot_call"]
        comm = table[i]["comm"]
        comm_size = len(comm_to_procs[comm])
        print(
            f"{call}".ljust(20)
            + f"{call_mean:.4f}".rjust(10)
            + f"{call_min:.4f}".rjust(10)
            + f"{call_max:.4f}".rjust(10)
            + f"{call_tot_vol}".rjust(15)
            + f"{nb_tot_call}".rjust(10)
            + f"{comm}".rjust(15)
            + f"{comm_size}".rjust(10)
        )

    return


def print_header():
     print("\n" + "".join(["*"] * 70))
     for i in range(0, 5):
        if i == 2:
            print("*" + "mpisee".center(68) + "*")
        else:
            print("*" + f"*".rjust(69))

     print("".join(["*"] * 70), end="\n\n")
    #tprint("mpisee")


def main():
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument("-f", "--full", action="store_true", help="print full info")
    my_parser.add_argument("-cct", action="store_true", help="print cct")
    my_parser.add_argument("-tc", action="store_true", help="print tc")
    my_parser.add_argument("-tc_limit", action="store", type=int, help="print tc")
    my_parser.add_argument(
        "-cct_limit",
        action="store",
        type=int,
        help="Limit number of communicators to ptint in cct view",
    )
    my_parser.add_argument(
        "file_path", metavar="path", type=str, help="Path to the csv file"
    )
    args = my_parser.parse_args()

    table, mapping, comm_to_procs,elapsed_time,mpi_times = prepare_data(args.file_path)

    print_header()
    print_mapping(mapping)
    print_times(elapsed_time,mpi_times)
    #print_elapsed_time(elapsed_time)
    #print_mpi_times(mpi_times)

    if args.cct:
        print_cct(table, comm_to_procs, args.cct_limit)

    if args.tc:
        print_tc(table, comm_to_procs, args.tc_limit)


# Main body
if __name__ == "__main__":
    main()
