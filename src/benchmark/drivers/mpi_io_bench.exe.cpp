#include "smesh_alloc.hpp"

#include "matrixio_array.h"

#include <mpi.h>

#include <algorithm>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>


// MPI_Info info;
// MPI_Info_create(&info);
// MPI_Info_set(info, "romio_cb_read", "enable");
// MPI_Info_set(info, "cb_nodes", "64");          // start with one aggregator per node
// MPI_Info_set(info, "cb_buffer_size", "16777216"); // 16 MiB, benchmark 16–64 MiB too

namespace {

struct ReadBuffer {
  void *data{nullptr};
  MPI_Win win{MPI_WIN_NULL};
  bool shared{false};
};

struct Range {
  ptrdiff_t begin{0};
  ptrdiff_t count{0};
};

Range split_range(const ptrdiff_t n, const int size, const int rank) {
  const ptrdiff_t base = n / size;
  const ptrdiff_t remainder = n - base * size;
  Range ret;
  ret.count = base + (rank < remainder ? 1 : 0);
  ret.begin = rank * base + std::min<ptrdiff_t>(rank, remainder);
  return ret;
}

MPI_Datatype parse_type(const char *name) {
  if (std::strcmp(name, "float32") == 0 || std::strcmp(name, "f32") == 0) {
    return MPI_FLOAT;
  }

  if (std::strcmp(name, "float64") == 0 || std::strcmp(name, "f64") == 0) {
    return MPI_DOUBLE;
  }

  if (std::strcmp(name, "int32") == 0 || std::strcmp(name, "i32") == 0) {
    return MPI_INT32_T;
  }

  if (std::strcmp(name, "int64") == 0 || std::strcmp(name, "i64") == 0) {
    return MPI_INT64_T;
  }

  if (std::strcmp(name, "uint8") == 0 || std::strcmp(name, "u8") == 0) {
    return MPI_UINT8_T;
  }

  if (std::strcmp(name, "uint16") == 0 || std::strcmp(name, "u16") == 0) {
    return MPI_UINT16_T;
  }

  if (std::strcmp(name, "uint32") == 0 || std::strcmp(name, "u32") == 0) {
    return MPI_UINT32_T;
  }

  if (std::strcmp(name, "uint64") == 0 || std::strcmp(name, "u64") == 0) {
    return MPI_UINT64_T;
  }

  return MPI_DATATYPE_NULL;
}

void *alloc_local_bytes(const ptrdiff_t count, const int type_size) {
  const size_t nbytes =
      static_cast<size_t>(std::max<ptrdiff_t>(count, 1)) * (size_t)type_size;
  return SMESH_ALLOC(nbytes);
}

long long sample_bytes(const void *data, const ptrdiff_t count,
                       const int type_size) {
  if (count <= 0) {
    return 0;
  }

  const unsigned char *bytes = static_cast<const unsigned char *>(data);
  const size_t nbytes = static_cast<size_t>(count) * (size_t)type_size;
  long long ret = 0;
  ret += bytes[0];
  ret += bytes[nbytes / 2];
  ret += bytes[nbytes - 1];
  return ret;
}

int read_mpi_file(MPI_Comm comm, const char *path, MPI_Datatype type,
                  const bool collective, const ptrdiff_t segment_elems,
                  ReadBuffer *buffer_out, ptrdiff_t *nlocal_out,
                  ptrdiff_t *nglobal_out) {
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  MPI_File file = MPI_FILE_NULL;
  MPI_Offset nbytes = 0;
  int type_size = 0;
  MPI_Status status;

  if (MPI_File_open(comm, path, MPI_MODE_RDONLY, MPI_INFO_NULL, &file) !=
      MPI_SUCCESS) {
    return 1;
  }

  if (MPI_File_get_size(file, &nbytes) != MPI_SUCCESS) {
    MPI_File_close(&file);
    return 1;
  }

  if (MPI_Type_size(type, &type_size) != MPI_SUCCESS) {
    MPI_File_close(&file);
    return 1;
  }

  const ptrdiff_t nglobal = static_cast<ptrdiff_t>(nbytes / type_size);
  if (nglobal * type_size != nbytes) {
    MPI_File_close(&file);
    return 1;
  }

  const Range local = split_range(nglobal, size, rank);
  void *data = alloc_local_bytes(local.count, type_size);

  const ptrdiff_t max_segment =
      segment_elems > 0 ? std::min<ptrdiff_t>(segment_elems, INT_MAX) : INT_MAX;

  if (collective) {
    const int local_rounds =
        static_cast<int>((local.count + max_segment - 1) / max_segment);
    int max_rounds = 0;
    MPI_Allreduce(&local_rounds, &max_rounds, 1, MPI_INT, MPI_MAX, comm);

    unsigned char dummy = 0;
    for (int round = 0; round < max_rounds; ++round) {
      const ptrdiff_t offset = static_cast<ptrdiff_t>(round) * max_segment;
      const ptrdiff_t remaining = local.count - offset;
      const int chunk = remaining > 0
                            ? static_cast<int>(
                                  std::min<ptrdiff_t>(max_segment, remaining))
                            : 0;
      MPI_Offset byte_offset = 0;
      void *buffer = &dummy;

      if (chunk > 0) {
        byte_offset = static_cast<MPI_Offset>(local.begin + offset) * type_size;
        buffer = &static_cast<unsigned char *>(data)[offset * type_size];
      }

      if (MPI_File_read_at_all(file, byte_offset, buffer, chunk, type,
                               &status) != MPI_SUCCESS) {
        SMESH_FREE(data);
        MPI_File_close(&file);
        return 1;
      }
    }
  } else {
    const int nrounds =
        static_cast<int>((local.count + max_segment - 1) / max_segment);
    for (int round = 0; round < nrounds; ++round) {
      const ptrdiff_t offset = static_cast<ptrdiff_t>(round) * max_segment;
      const ptrdiff_t remaining = local.count - offset;
      const int chunk =
          static_cast<int>(std::min<ptrdiff_t>(max_segment, remaining));
      const MPI_Offset byte_offset =
          static_cast<MPI_Offset>(local.begin + offset) * type_size;
      void *buffer = &static_cast<unsigned char *>(data)[offset * type_size];

      if (MPI_File_read_at(file, byte_offset, buffer, chunk, type, &status) !=
          MPI_SUCCESS) {
        SMESH_FREE(data);
        MPI_File_close(&file);
        return 1;
      }
    }
  }

  MPI_File_close(&file);
  buffer_out->data = data;
  buffer_out->win = MPI_WIN_NULL;
  buffer_out->shared = false;
  *nlocal_out = local.count;
  *nglobal_out = nglobal;
  return 0;
}

int read_node_aggregated(MPI_Comm comm, const char *path, MPI_Datatype type,
                         const ptrdiff_t segment_elems,
                         const int leaders_per_node, ReadBuffer *buffer_out,
                         ptrdiff_t *nlocal_out, ptrdiff_t *nglobal_out) {
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  MPI_Comm shared_comm = MPI_COMM_NULL;
  MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
                      &shared_comm);

  int local_rank = 0;
  int local_size = 1;
  MPI_Comm_rank(shared_comm, &local_rank);
  MPI_Comm_size(shared_comm, &local_size);

  if (leaders_per_node <= 0 || leaders_per_node > local_size) {
    MPI_Comm_free(&shared_comm);
    return 1;
  }

  int type_size = 0;
  MPI_Type_size(type, &type_size);

  MPI_File file = MPI_FILE_NULL;
  MPI_Offset nbytes = 0;
  if (MPI_File_open(comm, path, MPI_MODE_RDONLY, MPI_INFO_NULL, &file) !=
      MPI_SUCCESS) {
    MPI_Comm_free(&shared_comm);
    return 1;
  }

  if (MPI_File_get_size(file, &nbytes) != MPI_SUCCESS) {
    MPI_File_close(&file);
    MPI_Comm_free(&shared_comm);
    return 1;
  }
  MPI_File_close(&file);

  const ptrdiff_t nglobal = static_cast<ptrdiff_t>(nbytes / type_size);
  if (nglobal * type_size != nbytes) {
    MPI_Comm_free(&shared_comm);
    return 1;
  }

  const Range local = split_range(nglobal, size, rank);
  const MPI_Aint local_bytes =
      static_cast<MPI_Aint>(std::max<ptrdiff_t>(local.count, 1)) *
      static_cast<MPI_Aint>(type_size);

  MPI_Comm group_comm = MPI_COMM_NULL;
  MPI_Comm_split(shared_comm, local_rank % leaders_per_node, local_rank,
                 &group_comm);

  int group_rank = 0;
  int group_size = 1;
  MPI_Comm_rank(group_comm, &group_rank);
  MPI_Comm_size(group_comm, &group_size);

  MPI_Win win = MPI_WIN_NULL;
  void *data = nullptr;
  if (MPI_Win_allocate_shared(local_bytes, 1, MPI_INFO_NULL, group_comm, &data,
                              &win) != MPI_SUCCESS) {
    MPI_Comm_free(&group_comm);
    MPI_Comm_free(&shared_comm);
    return 1;
  }

  long long begin_ll = static_cast<long long>(local.begin);
  long long count_ll = static_cast<long long>(local.count);
  std::vector<long long> group_begin;
  std::vector<long long> group_count;
  if (group_rank == 0) {
    group_begin.resize(group_size);
    group_count.resize(group_size);
  }

  MPI_Gather(&begin_ll, 1, MPI_LONG_LONG, group_begin.data(), 1, MPI_LONG_LONG,
             0, group_comm);
  MPI_Gather(&count_ll, 1, MPI_LONG_LONG, group_count.data(), 1, MPI_LONG_LONG,
             0, group_comm);

  const ptrdiff_t max_segment =
      segment_elems > 0 ? std::min<ptrdiff_t>(segment_elems, INT_MAX) : INT_MAX;

  if (group_rank == 0) {
    MPI_File leader_file = MPI_FILE_NULL;
    if (MPI_File_open(MPI_COMM_SELF, path, MPI_MODE_RDONLY, MPI_INFO_NULL,
                      &leader_file) != MPI_SUCCESS) {
      MPI_Win_free(&win);
      MPI_Comm_free(&group_comm);
      MPI_Comm_free(&shared_comm);
      return 1;
    }

    std::vector<unsigned char *> group_ptr(group_size, nullptr);
    for (int peer = 0; peer < group_size; ++peer) {
      MPI_Aint peer_bytes = 0;
      int peer_disp_unit = 0;
      void *peer_ptr = nullptr;
      MPI_Win_shared_query(win, peer, &peer_bytes, &peer_disp_unit, &peer_ptr);
      group_ptr[peer] = static_cast<unsigned char *>(peer_ptr);
    }

    for (int peer = 0; peer < group_size; ++peer) {
      const ptrdiff_t peer_begin = static_cast<ptrdiff_t>(group_begin[peer]);
      const ptrdiff_t peer_count = static_cast<ptrdiff_t>(group_count[peer]);
      if (peer_count <= 0) {
        continue;
      }

      unsigned char *dest = group_ptr[peer];
      const int nrounds =
          static_cast<int>((peer_count + max_segment - 1) / max_segment);
      for (int round = 0; round < nrounds; ++round) {
        const ptrdiff_t off = static_cast<ptrdiff_t>(round) * max_segment;
        const int chunk = static_cast<int>(
            std::min<ptrdiff_t>(max_segment, peer_count - off));
        const MPI_Offset byte_offset =
            static_cast<MPI_Offset>(peer_begin + off) * type_size;
        if (MPI_File_read_at(leader_file, byte_offset, &dest[off * type_size],
                             chunk, type, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
          MPI_File_close(&leader_file);
          MPI_Win_free(&win);
          MPI_Comm_free(&group_comm);
          MPI_Comm_free(&shared_comm);
          return 1;
        }
      }
    }

    MPI_File_close(&leader_file);
  }

  MPI_Win_sync(win);
  MPI_Barrier(group_comm);
  MPI_Win_sync(win);

  MPI_Comm_free(&group_comm);
  MPI_Comm_free(&shared_comm);

  buffer_out->data = data;
  buffer_out->win = win;
  buffer_out->shared = true;
  *nlocal_out = local.count;
  *nglobal_out = nglobal;
  return 0;
}

int run_strategy(MPI_Comm comm, const char *path, MPI_Datatype type,
                 const std::string &strategy, const ptrdiff_t segment_elems,
                 const int leaders_per_node, ReadBuffer *buffer_out,
                 ptrdiff_t *nlocal_out, ptrdiff_t *nglobal_out) {
  if (strategy == "matrixio") {
    void *data = nullptr;
    const int err =
        array_create_from_file(comm, path, type, &data, nlocal_out, nglobal_out);
    buffer_out->data = data;
    buffer_out->win = MPI_WIN_NULL;
    buffer_out->shared = false;
    return err;
  }

  if (strategy == "matrixio_segmented") {
    const int segment = static_cast<int>(
        std::min<ptrdiff_t>(segment_elems > 0 ? segment_elems : INT_MAX, INT_MAX));
    void *data = nullptr;
    const int err = array_create_from_file_segmented(comm, path, type, &data,
                                                     segment, nlocal_out,
                                                     nglobal_out);
    buffer_out->data = data;
    buffer_out->win = MPI_WIN_NULL;
    buffer_out->shared = false;
    return err;
  }

  if (strategy == "mpi_collective") {
    return read_mpi_file(comm, path, type, true, segment_elems, buffer_out,
                         nlocal_out, nglobal_out);
  }

  if (strategy == "mpi_independent") {
    return read_mpi_file(comm, path, type, false, segment_elems, buffer_out,
                         nlocal_out, nglobal_out);
  }

  if (strategy == "node_aggregated") {
    return read_node_aggregated(comm, path, type, segment_elems,
                                leaders_per_node, buffer_out, nlocal_out,
                                nglobal_out);
  }

  return 1;
}

} // namespace

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (argc < 4 || argc > 6) {
    if (rank == 0) {
      std::fprintf(stderr,
                   "Usage: %s <input_path> <dtype> <strategy> [segment_elems] "
                   "[leaders_per_node]\n",
                   argv[0]);
    }
    MPI_Finalize();
    return 1;
  }

  MPI_Datatype type = parse_type(argv[2]);
  if (type == MPI_DATATYPE_NULL) {
    if (rank == 0) {
      std::fprintf(stderr, "Unsupported dtype: %s\n", argv[2]);
    }
    MPI_Finalize();
    return 1;
  }

  const std::string strategy = argv[3];
  const ptrdiff_t segment_elems =
      argc >= 5 ? static_cast<ptrdiff_t>(std::atoll(argv[4])) : INT_MAX;
  const int leaders_per_node = argc >= 6 ? std::atoi(argv[5]) : 1;

  MPI_Barrier(MPI_COMM_WORLD);
  const double t0 = MPI_Wtime();

  ReadBuffer buffer;
  ptrdiff_t nlocal = 0;
  ptrdiff_t nglobal = 0;
  const int err =
      run_strategy(MPI_COMM_WORLD, argv[1], type, strategy, segment_elems,
                   leaders_per_node, &buffer, &nlocal, &nglobal);

  MPI_Barrier(MPI_COMM_WORLD);
  const double local_dt = MPI_Wtime() - t0;

  int type_size = 0;
  MPI_Type_size(type, &type_size);

  const double local_max = local_dt;
  double elapsed = 0;
  MPI_Reduce(&local_max, &elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  long long local_bytes =
      static_cast<long long>(nlocal) * static_cast<long long>(type_size);
  long long total_bytes = 0;
  MPI_Reduce(&local_bytes, &total_bytes, 1, MPI_LONG_LONG, MPI_SUM, 0,
             MPI_COMM_WORLD);

  long long local_sample =
      err == 0 ? sample_bytes(buffer.data, nlocal, type_size) : 0;
  long long global_sample = 0;
  MPI_Reduce(&local_sample, &global_sample, 1, MPI_LONG_LONG, MPI_SUM, 0,
             MPI_COMM_WORLD);

  int global_err = 0;
  MPI_Allreduce(&err, &global_err, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  if (rank == 0) {
    const double gib_per_sec =
        elapsed > 0 ? (double)total_bytes / elapsed / (1024.0 * 1024.0 * 1024.0)
                    : 0.0;
    std::printf(
        "strategy=%s ranks=%d leaders_per_node=%d segment_elems=%lld "
        "nglobal=%lld total_bytes=%lld time_s=%.6f gib_s=%.6f sample=%lld "
        "status=%s\n",
        strategy.c_str(), size, leaders_per_node, (long long)segment_elems,
        (long long)nglobal, total_bytes, elapsed, gib_per_sec, global_sample,
        global_err == 0 ? "ok" : "error");
    std::fflush(stdout);
  }

  if (buffer.shared) {
    MPI_Win_free(&buffer.win);
  } else {
    SMESH_FREE(buffer.data);
  }
  MPI_Finalize();
  return global_err == 0 ? 0 : 1;
}
