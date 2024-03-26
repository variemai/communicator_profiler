#include "commprof.h"
#include "symbols.h"
#include "utils.h"
#include <mpi.h>


int
MPI_Win_create(void *base, MPI_Aint size, int disp_unit, MPI_Info info,
               MPI_Comm comm, MPI_Win *win)
{
    int ret;
    ret = PMPI_Win_create(base, size, disp_unit, info, comm, win);
    if (ret == MPI_SUCCESS)
    {
        MPI_Win_set_attr(*win, win_namekey(), comm);
    }
    return ret;
}

extern "C"{
void mpi_win_create_(void *base, MPI_Aint *size, int *disp_unit, MPI_Fint *info,
                     MPI_Fint *comm, MPI_Fint *win, MPI_Fint *ierr)
{
    int ret;
    MPI_Info c_info = MPI_Info_f2c(*info);
    MPI_Comm c_comm = MPI_Comm_f2c(*comm);
    MPI_Win c_win;

    ret = MPI_Win_create(base, *size, *disp_unit, c_info, c_comm, &c_win);
    *win = MPI_Win_c2f(c_win);
    *ierr = ret;
}
}


int
MPI_Put(const void *origin_addr, int origin_count,
        MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
        int target_count, MPI_Datatype target_datatype, MPI_Win win)
{
    int ret, flag;
    double t_elapsed;
    MPI_Comm comm;
    if (prof_enabled == 1){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Put(origin_addr, origin_count, origin_datatype, target_rank,
                   target_disp, target_count, target_datatype, win);
        t_elapsed = MPI_Wtime() - t_elapsed;
        MPI_Win_get_attr(win, win_namekey(), &comm, &flag);
        if (flag == 0){
            mcpt_abort("Window communicator not found\n");
        }
        else{
            profile_this(comm, origin_count, origin_datatype, Put ,t_elapsed, 0);
        }
    }
    else{
        ret = PMPI_Put(origin_addr, origin_count, origin_datatype, target_rank,
                   target_disp, target_count, target_datatype, win);
    }
    return ret;

}

extern "C" {
 void mpi_put_(const void *origin_addr, int *origin_count,
        MPI_Fint *origin_datatype, int *target_rank, int *target_disp,
        int *target_count, MPI_Fint *target_datatype, MPI_Fint *win, MPI_Fint *ierr)
{
    int ret;
    MPI_Win c_win = MPI_Win_f2c(*win);
    MPI_Datatype c_origin_datatype = MPI_Type_f2c(*origin_datatype);
    MPI_Datatype c_target_datatype = MPI_Type_f2c(*target_datatype);

    ret = MPI_Put(origin_addr, *origin_count, c_origin_datatype, *target_rank,
                  *target_disp, *target_count, c_target_datatype, c_win);
    *ierr = ret;
}

}

int
MPI_Accumulate(const void *origin_addr, int origin_count,
               MPI_Datatype origin_datatype, int target_rank,
               MPI_Aint target_disp, int target_count,
               MPI_Datatype target_datatype, MPI_Op op, MPI_Win win)
{
    int ret, flag;
    double t_elapsed;
    MPI_Comm comm;
    if (prof_enabled == 1){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Accumulate(origin_addr, origin_count, origin_datatype, target_rank,
                   target_disp, target_count, target_datatype, op, win);
        t_elapsed = MPI_Wtime() - t_elapsed;
        MPI_Win_get_attr(win, win_namekey(), &comm, &flag);
        if (flag == 0){
            mcpt_abort("Window communicator not found\n");
        }
        else{
            profile_this(comm, origin_count, origin_datatype, Accumulate ,t_elapsed, 0);
        }
    }
    else{
        ret = PMPI_Accumulate(origin_addr, origin_count, origin_datatype, target_rank,
                   target_disp, target_count, target_datatype, op, win);
    }
    return ret;
}

extern "C" {
void mpi_accumulate_(const void *origin_addr, int *origin_count,
               MPI_Fint *origin_datatype, int *target_rank,
               int *target_disp, int *target_count,
               MPI_Fint *target_datatype, MPI_Fint *op, MPI_Fint *win, MPI_Fint *ierr)
{
    int ret;
    MPI_Win c_win = MPI_Win_f2c(*win);
    MPI_Datatype c_origin_datatype = MPI_Type_f2c(*origin_datatype);
    MPI_Datatype c_target_datatype = MPI_Type_f2c(*target_datatype);
    MPI_Op c_op = MPI_Op_f2c(*op);

    ret = MPI_Accumulate(origin_addr, *origin_count, c_origin_datatype, *target_rank,
                  *target_disp, *target_count, c_target_datatype, c_op, c_win);
    *ierr = ret;
}
}


int
MPI_Get(void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
        int target_rank, MPI_Aint target_disp, int target_count,
        MPI_Datatype target_datatype, MPI_Win win)
{
    int ret, flag;
    double t_elapsed;
    MPI_Comm comm;
    if (prof_enabled == 1){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Get(origin_addr, origin_count, origin_datatype, target_rank,
                   target_disp, target_count, target_datatype, win);
        t_elapsed = MPI_Wtime() - t_elapsed;
        MPI_Win_get_attr(win, win_namekey(), &comm, &flag);
        if (flag == 0){
            mcpt_abort("Window communicator not found\n");
        }
        else{
            profile_this(comm, origin_count, origin_datatype, Get ,t_elapsed, 0);
        }
    }
    else{
        ret = PMPI_Get(origin_addr, origin_count, origin_datatype, target_rank,
                   target_disp, target_count, target_datatype, win);
    }
    return ret;
}

extern "C" {
void mpi_get_(void *origin_addr, int *origin_count, MPI_Fint *origin_datatype,
        int *target_rank, int *target_disp, int *target_count,
        MPI_Fint *target_datatype, MPI_Fint *win, MPI_Fint *ierr)
{
    int ret;
    MPI_Win c_win = MPI_Win_f2c(*win);
    MPI_Datatype c_origin_datatype = MPI_Type_f2c(*origin_datatype);
    MPI_Datatype c_target_datatype = MPI_Type_f2c(*target_datatype);

    ret = MPI_Get(origin_addr, *origin_count, c_origin_datatype, *target_rank,
                  *target_disp, *target_count, c_target_datatype, c_win);
    *ierr = ret;
}
}
