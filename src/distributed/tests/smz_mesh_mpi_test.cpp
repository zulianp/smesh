// Compressed mesh folders and fields (Mesh::write(path, MeshWriteOptions), Output), at several rank
// counts (ctest *_np1.._np4):
//  * a folder written with lossless compression reads back exactly like the raw folder, through the
//    distributed array readers and through Mesh::create_from_file (partitioning, ghosts, ...);
//  * a distributed Mesh::write with compression gives the same arrays as a distributed raw write;
//  * meta.yaml names the compressed files; nothing raw is left next to them;
//  * Output writes compressed nodal fields within their bound.

#include "smesh_byte_codec.hpp"
#include "smesh_communicator.hpp"
#include "smesh_conversion.hpp"
#include "smesh_compression_options.hpp"
#include "smesh_distributed_read.hpp"
#include "smesh_file_extensions.hpp"
#include "smesh_mesh.hpp"
#include "smesh_output.hpp"
#include "smesh_sideset.hpp"
#include "smesh_read.hpp"
#include "smesh_test.hpp"

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <vector>

using namespace smesh;
using namespace smesh::compression;
namespace fs = std::filesystem;

static std::string kDir;

template <class T>
static bool same(const T *a, const T *b, ptrdiff_t n) {
    return n == 0 || memcmp(a, b, size_t(n) * sizeof(T)) == 0;
}

static std::string read_text(const std::string &p) {
    std::ifstream     in(p);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

static MeshWriteOptions lossless_small_chunks() {
    CompressionOptions o;
    CompressionOptions::parse("lossless,chunk:257", &o);
    return MeshWriteOptions::all(o);
}

// The same arrays from a raw and a compressed folder, read with the distributed readers.
static int compare_folders(const std::shared_ptr<Communicator> &comm, const Path &a, const Path &b) {
    int        nxe_a = 0, nxe_b = 0, dim_a = 0, dim_b = 0;
    ptrdiff_t  nle_a = 0, nge_a = 0, nle_b = 0, nge_b = 0, nln_a = 0, ngn_a = 0, nln_b = 0, ngn_b = 0;
    idx_t    **ea = nullptr, **eb = nullptr;
    geom_t   **pa = nullptr, **pb = nullptr;
    SMESH_TEST_EQ(mesh_block_from_folder<idx_t>(comm->get(), a, &nxe_a, &ea, &nle_a, &nge_a), SMESH_SUCCESS);
    SMESH_TEST_EQ(mesh_block_from_folder<idx_t>(comm->get(), b, &nxe_b, &eb, &nle_b, &nge_b), SMESH_SUCCESS);
    SMESH_TEST_EQ(mesh_coordinates_from_folder<geom_t>(comm->get(), a, &dim_a, &pa, &nln_a, &ngn_a), SMESH_SUCCESS);
    SMESH_TEST_EQ(mesh_coordinates_from_folder<geom_t>(comm->get(), b, &dim_b, &pb, &nln_b, &ngn_b), SMESH_SUCCESS);
    SMESH_TEST_EQ(nxe_a, nxe_b);
    SMESH_TEST_EQ(nle_a, nle_b);
    SMESH_TEST_EQ(nge_a, nge_b);
    SMESH_TEST_EQ(dim_a, dim_b);
    SMESH_TEST_EQ(nln_a, nln_b);
    for (int v = 0; v < nxe_a; v++) SMESH_TEST_ASSERT(same(ea[v], eb[v], nle_a));
    for (int d = 0; d < dim_a; d++) SMESH_TEST_ASSERT(same(pa[d], pb[d], nln_a));
    for (int v = 0; v < nxe_a; v++) free(ea[v]), free(eb[v]);
    for (int d = 0; d < dim_a; d++) free(pa[d]), free(pb[d]);
    free(ea), free(eb), free(pa), free(pb);
    return SMESH_TEST_SUCCESS;
}

static int test_compressed_folder_reads_like_raw() {
    auto       comm = Communicator::world();
    const Path raw(kDir + "/raw"), cmp(kDir + "/cmp");
    if (comm->rank() == 0) {
        auto mesh = Mesh::create_hex8_cube(Communicator::self(), 9, 7, 5, 0, 0, 0, 1, 2, 3);
        SMESH_TEST_ASSERT(mesh != nullptr);
        auto left = Sideset::create_from_selector(mesh, [](geom_t x, geom_t, geom_t) { return x < geom_t(1e-6); });
        SMESH_TEST_ASSERT(!left.empty());
        mesh->add_sidesets("left", left);
        SMESH_TEST_EQ(mesh->write(raw, MeshWriteOptions()), SMESH_SUCCESS);
        SMESH_TEST_EQ(mesh->write(cmp, lossless_small_chunks()), SMESH_SUCCESS);
        // Only compressed arrays in the compressed folder, and meta.yaml names them.
        for (const auto &e : fs::directory_iterator(cmp.to_string())) {
            const std::string name = e.path().filename().string();
            if (name == "meta.yaml" || e.is_directory()) continue;
            SMESH_TEST_ASSERT(is_compressed_path(std::string_view(name)));
        }
        for (const auto &e : fs::recursive_directory_iterator(cmp.to_string() + "/sidesets")) {
            const std::string name = e.path().filename().string();
            if (name == "meta.yaml" || e.is_directory()) continue;
            SMESH_TEST_ASSERT(is_compressed_path(std::string_view(name)));
        }
        const std::string meta = read_text(cmp.to_string() + "/meta.yaml");
        SMESH_TEST_ASSERT(meta.find("i0.int32.smz") != std::string::npos || meta.find("i0.int64.smz") != std::string::npos);
        SMESH_TEST_ASSERT(meta.find("x.float") != std::string::npos && meta.find(".smz") != std::string::npos);
    }
    comm->barrier();
    if (compare_folders(comm, raw, cmp) != SMESH_TEST_SUCCESS) return SMESH_TEST_FAILURE;

    // The full pipeline: partitioning, ownership, ghosts.
    auto m_raw = Mesh::create_from_file(comm, raw);
    auto m_cmp = Mesh::create_from_file(comm, cmp);
    SMESH_TEST_ASSERT(m_raw && m_cmp);
    SMESH_TEST_EQ(m_raw->n_nodes(), m_cmp->n_nodes());
    SMESH_TEST_EQ(m_raw->n_elements(), m_cmp->n_elements());
    SMESH_TEST_EQ(m_raw->sidesets().size(), m_cmp->sidesets().size());
    SMESH_TEST_ASSERT(!m_cmp->sidesets().empty());
    for (size_t k = 0; k < m_raw->sidesets().size(); k++) {
        auto a = m_raw->sidesets()[k].second, b = m_cmp->sidesets()[k].second;
        SMESH_TEST_EQ(a->size(), b->size());
        SMESH_TEST_ASSERT(same(a->parent()->data(), b->parent()->data(), a->size()));
        SMESH_TEST_ASSERT(same(a->lfi()->data(), b->lfi()->data(), a->size()));
    }
    for (int d = 0; d < m_raw->spatial_dimension(); d++)
        SMESH_TEST_ASSERT(same(m_raw->points()->data()[d], m_cmp->points()->data()[d], m_raw->n_nodes()));
    return SMESH_TEST_SUCCESS;
}

// Distributed writes: compressed and raw give the same arrays.
static int test_distributed_write() {
    auto       comm = Communicator::world();
    const Path raw(kDir + "/raw");
    auto       mesh = Mesh::create_from_file(comm, raw);
    SMESH_TEST_ASSERT(mesh != nullptr);
    const Path out_raw(kDir + "/dist_raw"), out_cmp(kDir + "/dist_cmp");
    SMESH_TEST_EQ(mesh->write(out_raw, MeshWriteOptions()), SMESH_SUCCESS);
    CompressionOptions z;
    SMESH_TEST_EQ(CompressionOptions::parse("lossless,chunk:100", &z), SMESH_SUCCESS);
    if (byte_stage_available(ByteStage::Zstd)) z.byte_stage = ByteStage::Zstd;
    SMESH_TEST_EQ(mesh->write(out_cmp, MeshWriteOptions::all(z)), SMESH_SUCCESS);
    comm->barrier();
    return compare_folders(comm, out_raw, out_cmp);
}

// Compressed nodal fields, raw for reference.
static int test_output_fields() {
    auto       comm = Communicator::world();
    const Path raw(kDir + "/raw");
    auto       mesh = Mesh::create_from_file(comm, raw);
    SMESH_TEST_ASSERT(mesh != nullptr);
    const ptrdiff_t  n = mesh->n_nodes();
    std::vector<f64> u(static_cast<size_t>(n));
    auto             p = mesh->points()->data();
    for (ptrdiff_t i = 0; i < n; i++) u[size_t(i)] = std::sin(3.0 * p[0][i]) + p[1][i] * p[2][i];
    auto out = Output::create(mesh, Path(kDir + "/fields"));
    SMESH_TEST_EQ(out->write_nodal("u", SMESH_FLOAT64, u.data()), SMESH_SUCCESS);
    CompressionOptions o;
    SMESH_TEST_EQ(CompressionOptions::parse("range:1e-4,chunk:200", &o), SMESH_SUCCESS);
    out->set_compression(o);
    SMESH_TEST_EQ(out->write_nodal("uc", SMESH_FLOAT64, u.data()), SMESH_SUCCESS);
    comm->barrier();
    if (comm->rank() == 0) {
        f64      *a = nullptr, *b = nullptr;
        ptrdiff_t na = 0, nb = 0;
        SMESH_TEST_EQ(array_read_convert_from_extension<f64>(Path(kDir + "/fields/u.float64"), &a, &na), SMESH_SUCCESS);
        SMESH_TEST_EQ(array_read_convert_from_extension<f64>(Path(kDir + "/fields/uc.float64.smz"), &b, &nb),
                      SMESH_SUCCESS);
        SMESH_TEST_EQ(na, nb);
        double lo = a[0], hi = a[0];
        for (ptrdiff_t i = 0; i < na; i++) lo = std::min(lo, a[i]), hi = std::max(hi, a[i]);
        for (ptrdiff_t i = 0; i < na; i++) SMESH_TEST_ASSERT(std::fabs(a[i] - b[i]) <= 1e-4 * (hi - lo));
        free(a), free(b);
    }
    comm->barrier();
    return SMESH_TEST_SUCCESS;
}

static std::vector<u8> slurp(const std::string &p) {
    std::ifstream in(p, std::ios::binary);
    return std::vector<u8>((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
}

// Bound mode h: node heights come from the elements of every rank (minimum over ranks), so a
// distributed write gives the same bytes as the serial one.
static int test_height_bounded_coordinates_rank_independent() {
    auto       comm = Communicator::world();
    const Path graded(kDir + "/graded"), ser(kDir + "/h_serial"), dist(kDir + "/h_dist");
    CompressionOptions h;
    SMESH_TEST_EQ(CompressionOptions::parse("h:1e-3,chunk:300", &h), SMESH_SUCCESS);
    MeshWriteOptions o = MeshWriteOptions::all(CompressionOptions::lossless());
    o.coordinates      = h;
    if (comm->rank() == 0) {
        auto     mesh = Mesh::create_hex8_cube(Communicator::self(), 11, 9, 7);
        geom_t **p    = mesh->points()->data();
        for (ptrdiff_t i = 0; i < mesh->n_nodes(); i++) {
            const double x = p[0][i], y = p[1][i];
            p[0][i]        = geom_t(x * x * x);
            p[1][i]        = geom_t(y + 0.2 * x * x);
        }
        SMESH_TEST_EQ(mesh->write(graded, MeshWriteOptions()), SMESH_SUCCESS);
        auto back = Mesh::create_from_file(Communicator::self(), graded);
        SMESH_TEST_EQ(back->write(ser, o), SMESH_SUCCESS);
    }
    comm->barrier();
    auto mesh = Mesh::create_from_file(comm, graded);
    SMESH_TEST_ASSERT(mesh != nullptr);
    SMESH_TEST_EQ(mesh->write(dist, o), SMESH_SUCCESS);
    comm->barrier();
    if (comm->rank() == 0) {
        for (const char *f : {"x.float32.smz", "y.float32.smz", "z.float32.smz"}) {
            const auto a = slurp(ser.to_string() + "/" + f), b = slurp(dist.to_string() + "/" + f);
            SMESH_TEST_ASSERT(!a.empty());
            if (a != b) {
                fprintf(stderr, "%s: distributed file differs from the serial one\n", f);
                return SMESH_TEST_FAILURE;
            }
        }
    }
    comm->barrier();
    return SMESH_TEST_SUCCESS;
}

// A HEX8 + TET4 multiblock mesh (the distributed reader takes the SFC-partitioned path): lossless
// compressed folders read like raw ones, h-bounded coordinates stay within their bound.
static int test_multiblock() {
    auto       comm = Communicator::world();
    const Path raw(kDir + "/mb_raw"), cmp(kDir + "/mb_cmp"), lossy(kDir + "/mb_h");
    CompressionOptions h;
    SMESH_TEST_EQ(CompressionOptions::parse("h:1e-3,chunk:500", &h), SMESH_SUCCESS);
    MeshWriteOptions oh = MeshWriteOptions::all(CompressionOptions::lossless());
    oh.coordinates      = h;
    if (comm->rank() == 0) {
        auto            cube   = Mesh::create_hex8_cube(Communicator::self(), 8, 4, 3, 0, 0, 0, 2, 1, 1);
        const ptrdiff_t n_all  = cube->n_elements();
        const ptrdiff_t n_keep = n_all / 2, n_conv = n_all - n_keep;
        auto            src    = cube->elements(0)->data();
        auto            hex    = create_host_buffer<idx_t>(8, size_t(n_keep));
        for (int d = 0; d < 8; ++d) std::memcpy(hex->data()[d], src[d], size_t(n_keep) * sizeof(idx_t));
        idx_t *tail[8];
        for (int d = 0; d < 8; ++d) tail[d] = src[d] + n_keep;
        auto tet = create_host_buffer<idx_t>(4, size_t(n_conv * 6));
        mesh_hex8_to_6x_tet4<idx_t>(n_conv, tail, tet->data());
        std::vector<std::shared_ptr<Mesh::Block>> blocks;
        blocks.push_back(std::make_shared<Mesh::Block>("hex", HEX8, hex));
        blocks.push_back(std::make_shared<Mesh::Block>("tet", TET4, tet));
        Mesh serial(Communicator::self(), blocks, cube->points());
        SMESH_TEST_EQ(serial.write(raw, MeshWriteOptions()), SMESH_SUCCESS);
        SMESH_TEST_EQ(serial.write(cmp, lossless_small_chunks()), SMESH_SUCCESS);
        SMESH_TEST_EQ(serial.write(lossy, oh), SMESH_SUCCESS);
        const std::string meta = read_text(cmp.to_string() + "/meta.yaml");
        SMESH_TEST_ASSERT(meta.find("blocks/hex/i0.") != std::string::npos && meta.find(".smz") != std::string::npos);
    }
    comm->barrier();
    auto a = Mesh::create_from_file(comm, raw);
    auto b = Mesh::create_from_file(comm, cmp);
    auto c = Mesh::create_from_file(comm, lossy);
    SMESH_TEST_ASSERT(a && b && c);
    SMESH_TEST_EQ(a->n_blocks(), b->n_blocks());
    SMESH_TEST_EQ(a->n_nodes(), b->n_nodes());
    SMESH_TEST_EQ(a->n_nodes(), c->n_nodes());
    for (size_t k = 0; k < a->n_blocks(); k++) {
        auto ea = a->block(k)->elements(), eb = b->block(k)->elements();
        SMESH_TEST_EQ(ea->extent(1), eb->extent(1));
        for (size_t v = 0; v < size_t(ea->extent(0)); v++)
            SMESH_TEST_ASSERT(same(ea->data()[v], eb->data()[v], ea->extent(1)));
    }
    for (int d = 0; d < 3; d++) {
        SMESH_TEST_ASSERT(same(a->points()->data()[d], b->points()->data()[d], a->n_nodes()));
        // h-bounded: within 1e-3 of the smallest height; the smallest height here is >= 1/(8 sqrt 3)
        for (ptrdiff_t i = 0; i < a->n_nodes(); i++)
            SMESH_TEST_ASSERT(std::fabs(double(c->points()->data()[d][i]) - double(a->points()->data()[d][i])) <=
                              1e-3 * 0.25);
    }
    // Distributed compressed write of the multiblock mesh == distributed raw write.
    const Path dr(kDir + "/mb_dist_raw"), dc(kDir + "/mb_dist_cmp");
    SMESH_TEST_EQ(a->write(dr, MeshWriteOptions()), SMESH_SUCCESS);
    SMESH_TEST_EQ(a->write(dc, lossless_small_chunks()), SMESH_SUCCESS);
    comm->barrier();
    auto x = Mesh::create_from_file(comm, dr), y = Mesh::create_from_file(comm, dc);
    SMESH_TEST_ASSERT(x && y);
    for (int d = 0; d < 3; d++) SMESH_TEST_ASSERT(same(x->points()->data()[d], y->points()->data()[d], x->n_nodes()));
    return SMESH_TEST_SUCCESS;
}

int main(int argc, char *argv[]) {
    SMESH_UNIT_TEST_INIT(argc, argv);
    auto comm = Communicator::world();
    kDir      = "smz_mesh_test_tmp_np" + std::to_string(comm->size());
    if (comm->rank() == 0) {
        fs::remove_all(kDir);
        fs::create_directories(kDir);
    }
    comm->barrier();
    SMESH_RUN_TEST(test_compressed_folder_reads_like_raw);
    SMESH_RUN_TEST(test_distributed_write);
    SMESH_RUN_TEST(test_output_fields);
    SMESH_RUN_TEST(test_height_bounded_coordinates_rank_independent);
    SMESH_RUN_TEST(test_multiblock);
    int any = err, all = 0;
    MPI_Allreduce(&any, &all, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (comm->rank() == 0 && !all) fs::remove_all(kDir);
    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}
