#include "smesh_mesh.hpp"
#include "smesh_adjacency.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_file_extensions.hpp"
#include "smesh_path.hpp"
#include "smesh_edgeset.hpp"
#include "smesh_edgesets.hpp"
#include "smesh_nodeset.hpp"
#include "smesh_sideset.hpp"
#include "smesh_sidesets.hpp"
#include "smesh_tracer.hpp"
#include "smesh_types.hpp"
#include "smesh_write.hpp"

#include <cstring>

#ifdef SMESH_ENABLE_MPI
#include "smesh_distributed_write.hpp"
#endif

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace smesh {

static const char *xdmf_endianess() {
#if defined(__BYTE_ORDER__) && defined(__ORDER_BIG_ENDIAN__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    return "Big";
#else
    return "Little";
#endif
}

static const char *xdmf_topology_type(const enum ElemType et) {
    switch (et) {
        case EDGE2:
        case EDGESHELL2:
        case BEAM2:
            return "Polyline";
        case NODE1:
            return "Polyvertex";
        case EDGE3:
        case EDGESHELL3:
            return "Edge_3";
        case TRI3:
        case TRISHELL3:
        case MACRO_TRI3:
        case MACRO_TRISHELL3:
            return "Triangle";
        case TRI6:
        case TRISHELL6:
            return "Triangle_6";
        case QUAD4:
        case QUADSHELL4:
        case PROTEUS_QUAD4:
        case PROTEUS_QUADSHELL4:
            return "Quadrilateral";
        case QUAD9:
        case QUADSHELL9:
        case PROTEUS_QUAD9:
        case PROTEUS_QUADSHELL9:
            return "Quadrilateral_9";
        case TET4:
        case MACRO_TET4:
        case PROTEUS_TET4:
            return "Tetrahedron";
        case TET10:
        case PROTEUS_TET10:
            return "Tetrahedron_10";
        case PYRAMID5:
        case PROTEUS_PYRAMID5:
            return "Pyramid";
        case WEDGE6:
        case PROTEUS_WEDGE6:
            return "Wedge";
        case HEX8:
        case PROTEUS_HEX8:
            return "Hexahedron";
        case HEX27:
        case PROTEUS_HEX27:
            return "Hexahedron_27";
        default:
            break;
    }
    if (is_semistructured_type(et)) {
        const int           nxe = elem_num_nodes(et);
        const enum ElemType fam = ss_source_family(et);
        if (fam == HEX8 && nxe == 8) {
            return "Hexahedron";
        }
        if (fam == HEX8 && nxe == 27) {
            return "Hexahedron_27";
        }
        if (fam == TET4 && nxe == 4) {
            return "Tetrahedron";
        }
        if (fam == TET4 && nxe == 10) {
            return "Tetrahedron_10";
        }
        if (fam == QUAD4 && nxe == 4) {
            return "Quadrilateral";
        }
        if (fam == QUAD4 && nxe == 9) {
            return "Quadrilateral_9";
        }
        if (fam == WEDGE6 && nxe == 6) {
            return "Wedge";
        }
        if (fam == PYRAMID5 && nxe == 5) {
            return "Pyramid";
        }
    }
    return nullptr;
}

static std::string xdmf_typed_name(const char *stem) {
    return std::string(stem) + "." + std::string(TypeToString<idx_t>::value());
}

static std::string xdmf_coord_name(const char axis) {
    return std::string(1, axis) + "." + std::string(TypeToString<geom_t>::value());
}

static Path connectivity_aos_path(const Path &block_folder) {
    return block_folder / Path(xdmf_typed_name("connectivity"));
}

struct XdmfSurfGrid {
    std::string name;
    std::string conn;
    const char *topo;
    ptrdiff_t   n;
    int         nxe;
};

static void write_xdmf_uniform_grid(FILE           *fp,
                                    const char     *name,
                                    const char     *topo,
                                    const ptrdiff_t ne,
                                    const int       nxe,
                                    const char     *conn,
                                    const char     *endian) {
    fprintf(fp, "    <Grid Name=\"%s\" GridType=\"Uniform\">\n", name);
    fprintf(fp, "      <Geometry Reference=\"/Xdmf/Domain/Geometry[1]\"/>\n");
    if (std::strcmp(topo, "Polyline") == 0 || std::strcmp(topo, "Polyvertex") == 0) {
        fprintf(fp,
                "      <Topology TopologyType=\"%s\" NumberOfElements=\"%ld\" NodesPerElement=\"%d\">\n",
                topo,
                (long)ne,
                nxe);
    } else {
        fprintf(fp,
                "      <Topology TopologyType=\"%s\" NumberOfElements=\"%ld\">\n",
                topo,
                (long)ne);
    }
    fprintf(fp,
            "        <DataItem Format=\"Binary\" Dimensions=\"%ld %d\" NumberType=\"Int\" "
            "Precision=\"%d\" Endian=\"%s\">%s</DataItem>\n",
            (long)ne,
            nxe,
            (int)sizeof(idx_t),
            endian,
            conn);
    fprintf(fp, "      </Topology>\n");
    fprintf(fp, "    </Grid>\n");
}

static int write_block_connectivity_aos(const Mesh       &mesh,
                                        const Path       &block_folder,
                                        const block_idx_t bid,
                                        const ptrdiff_t   n_file_elements,
                                        const int         from_memory) {
    auto block = mesh.block(bid);
    if (!block) {
        return SMESH_FAILURE;
    }
    const int  nxe = elem_num_nodes(block->element_type());
    const Path aos = connectivity_aos_path(block_folder);
    if (from_memory) {
        return mesh_write_soa_to_aos(aos, nxe, n_file_elements, block->elements()->data());
    }
    return mesh_write_soa_files_to_aos(block_folder, aos, nxe, n_file_elements);
}

static int write_mesh_xdmf(const Mesh                      &mesh,
                           const Path                      &folder,
                           const ptrdiff_t                 *n_elem_file,
                           const ptrdiff_t                  n_nodes_file,
                           const std::vector<XdmfSurfGrid> &surfs) {
    const int  dim       = mesh.spatial_dimension();
    const int  n_blocks  = (int)mesh.n_blocks();
    const Path xdmf_path = folder / "mesh.xdmf";
    FILE      *fp        = fopen(xdmf_path.c_str(), "w");
    if (!fp) {
        SMESH_ERROR("write_with_xdmf: Unable to write %s\n", xdmf_path.c_str());
        return SMESH_FAILURE;
    }

    const char *endian = xdmf_endianess();
    fprintf(fp,
            "<?xml version=\"1.0\" ?>\n"
            "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n"
            "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">\n"
            "  <Domain>\n");

    if (dim >= 3) {
        fprintf(fp, "    <Geometry name=\"geo\" GeometryType=\"X_Y_Z\">\n");
        for (int d = 0; d < 3; ++d) {
            const char axis = (char)('x' + d);
            fprintf(fp,
                    "      <DataItem Format=\"Binary\" Dimensions=\"%ld\" NumberType=\"Float\" "
                    "Precision=\"%d\" Endian=\"%s\">%s</DataItem>\n",
                    (long)n_nodes_file,
                    (int)sizeof(geom_t),
                    endian,
                    xdmf_coord_name(axis).c_str());
        }
        fprintf(fp, "    </Geometry>\n");
    } else {
        fprintf(fp, "    <Geometry name=\"geo\" GeometryType=\"XY\">\n");
        for (int d = 0; d < std::max(dim, 1); ++d) {
            const char axis = (char)('x' + d);
            fprintf(fp,
                    "      <DataItem Format=\"Binary\" Dimensions=\"%ld\" NumberType=\"Float\" "
                    "Precision=\"%d\" Endian=\"%s\">%s</DataItem>\n",
                    (long)n_nodes_file,
                    (int)sizeof(geom_t),
                    endian,
                    xdmf_coord_name(axis).c_str());
        }
        if (dim < 2) {
            fprintf(fp,
                    "      <DataItem Format=\"Binary\" Dimensions=\"%ld\" NumberType=\"Float\" "
                    "Precision=\"%d\" Endian=\"%s\">%s</DataItem>\n",
                    (long)n_nodes_file,
                    (int)sizeof(geom_t),
                    endian,
                    xdmf_coord_name('x').c_str());
        }
        fprintf(fp, "    </Geometry>\n");
    }

    const int collection = n_blocks > 1 || !surfs.empty();
    if (collection) {
        fprintf(fp, "    <Grid Name=\"mesh\" GridType=\"Collection\" CollectionType=\"Spatial\">\n");
    }

    for (int b = 0; b < n_blocks; ++b) {
        auto block = mesh.block(b);
        if (!block) {
            continue;
        }
        const enum ElemType et   = block->element_type();
        const char         *topo = xdmf_topology_type(et);
        if (!topo) {
            fprintf(stderr,
                    "write_with_xdmf: no XDMF topology for %s (block %s); connectivity file still written\n",
                    type_to_string(et),
                    block->name().c_str());
            continue;
        }
        const int         nxe  = elem_num_nodes(et);
        const ptrdiff_t   ne   = n_elem_file[b];
        const std::string conn = n_blocks == 1
                                         ? xdmf_typed_name("connectivity")
                                         : (std::string("blocks/") + block->name() + "/" +
                                            xdmf_typed_name("connectivity"));
        write_xdmf_uniform_grid(fp, block->name().c_str(), topo, ne, nxe, conn.c_str(), endian);
    }

    for (size_t i = 0; i < surfs.size(); ++i) {
        write_xdmf_uniform_grid(fp,
                                surfs[i].name.c_str(),
                                surfs[i].topo,
                                surfs[i].n,
                                surfs[i].nxe,
                                surfs[i].conn.c_str(),
                                endian);
    }

    if (collection) {
        fprintf(fp, "    </Grid>\n");
    }
    fprintf(fp, "  </Domain>\n</Xdmf>\n");
    fclose(fp);
    return SMESH_SUCCESS;
}

static int write_xdmf_surface_part(const Mesh                &mesh,
                                   const Path                &root,
                                   const std::string         &grid_name,
                                   const std::string         &rel_conn,
                                   const enum ElemType        vol_et,
                                   const idx_t *const        *elems,
                                   const element_idx_t       *parent,
                                   const i16                 *lfi,
                                   const ptrdiff_t            n_local,
                                   const int                  nnxs,
                                   const char                *topo,
                                   std::vector<XdmfSurfGrid> *surfs) {
    auto      comm     = mesh.comm();
    const i64 n_global = comm->sum((i64)n_local);
    if (n_global <= 0) {
        return SMESH_SUCCESS;
    }

    auto soa = create_host_buffer<idx_t>(nnxs, n_local);
    if (n_local > 0) {
        if (extract_surface_from_sideset(vol_et, elems, n_local, parent, lfi, soa->data()) !=
            SMESH_SUCCESS) {
            return SMESH_FAILURE;
        }
        if (comm->size() > 1 && mesh.distributed() && mesh.distributed()->node_mapping()) {
            const large_idx_t *const map   = mesh.distributed()->node_mapping()->data();
            const ptrdiff_t          n_map = mesh.n_nodes();
            idx_t *const            *d     = soa->data();
            for (int c = 0; c < nnxs; ++c) {
                idx_t *col = d[c];
                for (ptrdiff_t i = 0; i < n_local; ++i) {
                    const idx_t n = col[i];
                    if (n < 0 || (ptrdiff_t)n >= n_map) {
                        SMESH_ERROR("write_with_xdmf: surface node %ld out of range\n", (long)n);
                        return SMESH_FAILURE;
                    }
                    col[i] = (idx_t)map[n];
                }
            }
        }
    }

    const Path aos_path = root / Path(rel_conn);
    int        err      = SMESH_SUCCESS;
#ifdef SMESH_ENABLE_MPI
    if (comm->size() > 1) {
        auto          aos = create_host_buffer<idx_t>((size_t)n_local * (size_t)nnxs);
        idx_t        *out = n_local > 0 ? aos->data() : nullptr;
        idx_t *const *in  = n_local > 0 ? soa->data() : nullptr;
        for (ptrdiff_t j = 0; j < n_local; ++j) {
            for (int c = 0; c < nnxs; ++c) {
                out[j * nnxs + c] = in[c][j];
            }
        }
        err = array_write_convert_from_extension(comm->get(),
                                                 aos_path,
                                                 out,
                                                 n_local * (ptrdiff_t)nnxs,
                                                 (ptrdiff_t)n_global * (ptrdiff_t)nnxs);
    } else
#endif
    {
        err = mesh_write_soa_to_aos(aos_path, nnxs, n_local, soa->data());
    }
    if (err != SMESH_SUCCESS) {
        return err;
    }
    if (comm->rank() == 0) {
        XdmfSurfGrid g;
        g.name = grid_name;
        g.conn = rel_conn;
        g.topo = topo;
        g.n    = (ptrdiff_t)n_global;
        g.nxe  = nnxs;
        surfs->push_back(std::move(g));
    }
    return SMESH_SUCCESS;
}

static int write_xdmf_soa_grid(const Mesh                &mesh,
                               const Path                &root,
                               const std::string         &grid_name,
                               const std::string         &rel_conn,
                               idx_t *const              *soa,
                               const ptrdiff_t            n_local,
                               const int                  nnxs,
                               const char                *topo,
                               std::vector<XdmfSurfGrid> *surfs) {
    auto      comm     = mesh.comm();
    const i64 n_global = comm->sum((i64)n_local);
    if (n_global <= 0) {
        return SMESH_SUCCESS;
    }

    if (n_local > 0 && comm->size() > 1 && mesh.distributed() && mesh.distributed()->node_mapping()) {
        const large_idx_t *const map   = mesh.distributed()->node_mapping()->data();
        const ptrdiff_t          n_map = mesh.n_nodes();
        for (int c = 0; c < nnxs; ++c) {
            idx_t *col = soa[c];
            for (ptrdiff_t i = 0; i < n_local; ++i) {
                const idx_t n = col[i];
                if (n < 0 || (ptrdiff_t)n >= n_map) {
                    SMESH_ERROR("write_with_xdmf: node %ld out of range\n", (long)n);
                    return SMESH_FAILURE;
                }
                col[i] = (idx_t)map[n];
            }
        }
    }

    const Path aos_path = root / Path(rel_conn);
    int        err      = SMESH_SUCCESS;
#ifdef SMESH_ENABLE_MPI
    if (comm->size() > 1) {
        auto          aos = create_host_buffer<idx_t>((size_t)n_local * (size_t)nnxs);
        idx_t        *out = n_local > 0 ? aos->data() : nullptr;
        for (ptrdiff_t j = 0; j < n_local; ++j) {
            for (int c = 0; c < nnxs; ++c) {
                out[j * nnxs + c] = soa[c][j];
            }
        }
        err = array_write_convert_from_extension(comm->get(),
                                                 aos_path,
                                                 out,
                                                 n_local * (ptrdiff_t)nnxs,
                                                 (ptrdiff_t)n_global * (ptrdiff_t)nnxs);
    } else
#endif
    {
        err = mesh_write_soa_to_aos(aos_path, nnxs, n_local, soa);
    }
    if (err != SMESH_SUCCESS) {
        return err;
    }
    if (comm->rank() == 0) {
        XdmfSurfGrid g;
        g.name = grid_name;
        g.conn = rel_conn;
        g.topo = topo;
        g.n    = (ptrdiff_t)n_global;
        g.nxe  = nnxs;
        surfs->push_back(std::move(g));
    }
    return SMESH_SUCCESS;
}

static int write_xdmf_edgeset_grids(const Mesh                &mesh,
                                    const Path                &root,
                                    std::vector<XdmfSurfGrid> *surfs) {
    const auto &reg = mesh.edgesets();
    if (reg.empty()) {
        return SMESH_SUCCESS;
    }

    std::vector<std::string>                           names;
    std::vector<std::vector<std::shared_ptr<Edgeset>>> groups;
    for (size_t i = 0; i < reg.size(); ++i) {
        size_t gi = names.size();
        for (size_t g = 0; g < names.size(); ++g) {
            if (names[g] == reg[i].first) {
                gi = g;
                break;
            }
        }
        if (gi == names.size()) {
            names.push_back(reg[i].first);
            groups.emplace_back();
        }
        groups[gi].push_back(reg[i].second);
    }

    int err = SMESH_SUCCESS;
    for (size_t g = 0; err == SMESH_SUCCESS && g < names.size(); ++g) {
        const int multi = groups[g].size() > 1;
        for (size_t k = 0; err == SMESH_SUCCESS && k < groups[g].size(); ++k) {
            const auto &es = groups[g][k];
            if (!es) {
                continue;
            }
            auto block = mesh.block(es->block_id());
            if (!block) {
                err = SMESH_FAILURE;
                break;
            }
            const enum ElemType et = block->element_type();
            if (is_semistructured_type(et)) {
                if (mesh.comm()->rank() == 0) {
                    fprintf(stderr,
                            "write_with_xdmf: skipping SS edgeset '%s' (%s)\n",
                            names[g].c_str(),
                            type_to_string(et));
                }
                continue;
            }
            LocalEdgeTable let;
            if (let.fill(et) != SMESH_SUCCESS) {
                if (mesh.comm()->rank() == 0) {
                    LocalEdgeTable::report_unsupported("write_with_xdmf", et);
                }
                continue;
            }
            const int         nne  = let.nnxe > 0 ? let.nnxe : 2;
            const char       *topo = (nne == 3) ? "Edge_3" : "Polyline";
            std::string       leaf = std::string("edgesets/") + names[g];
            std::string       grid = names[g];
            if (multi) {
                leaf += "/" + std::to_string((long long)es->block_id());
                grid += "_" + std::to_string((long long)es->block_id());
            }
            const std::string rel = leaf + "/edge." + std::string(TypeToString<idx_t>::value());
            const ptrdiff_t   n   = es->size();
            auto              soa = create_host_buffer<idx_t>(nne, n);
            if (n > 0) {
                if (extract_edges_from_edgeset(et,
                                               block->elements()->data(),
                                               n,
                                               es->parent()->data(),
                                               es->lei()->data(),
                                               soa->data()) != SMESH_SUCCESS) {
                    err = SMESH_FAILURE;
                    break;
                }
            }
            err = write_xdmf_soa_grid(mesh, root, grid, rel, soa->data(), n, nne, topo, surfs);
        }
    }
    return err;
}

static int write_xdmf_nodeset_grids(const Mesh                &mesh,
                                    const Path                &root,
                                    std::vector<XdmfSurfGrid> *surfs) {
    const auto &reg = mesh.nodesets();
    if (reg.empty()) {
        return SMESH_SUCCESS;
    }

    std::vector<std::string>                           names;
    std::vector<std::vector<std::shared_ptr<Nodeset>>> groups;
    for (size_t i = 0; i < reg.size(); ++i) {
        size_t gi = names.size();
        for (size_t g = 0; g < names.size(); ++g) {
            if (names[g] == reg[i].first) {
                gi = g;
                break;
            }
        }
        if (gi == names.size()) {
            names.push_back(reg[i].first);
            groups.emplace_back();
        }
        groups[gi].push_back(reg[i].second);
    }

    int err = SMESH_SUCCESS;
    for (size_t g = 0; err == SMESH_SUCCESS && g < names.size(); ++g) {
        const int multi = groups[g].size() > 1;
        for (size_t k = 0; err == SMESH_SUCCESS && k < groups[g].size(); ++k) {
            const auto &ns = groups[g][k];
            if (!ns) {
                continue;
            }
            std::string leaf = std::string("nodesets/") + names[g];
            std::string grid = names[g];
            if (multi) {
                leaf += "/" + std::to_string((long long)k);
                grid += "_" + std::to_string((long long)k);
            }
            const std::string rel = leaf + "/nodes." + std::string(TypeToString<idx_t>::value());
            const ptrdiff_t   n   = ns->size();
            auto              soa = create_host_buffer<idx_t>(1, n);
            if (n > 0) {
                std::memcpy(soa->data()[0], ns->nodes()->data(), (size_t)n * sizeof(idx_t));
            }
            err = write_xdmf_soa_grid(mesh, root, grid, rel, soa->data(), n, 1, "Polyvertex", surfs);
        }
    }
    return err;
}

static int write_xdmf_sideset_surfaces(const Mesh                &mesh,
                                       const Path                &root,
                                       std::vector<XdmfSurfGrid> *surfs) {
    const auto &reg = mesh.sidesets();
    if (reg.empty()) {
        return SMESH_SUCCESS;
    }

    std::vector<std::string>                           names;
    std::vector<std::vector<std::shared_ptr<Sideset>>> groups;
    for (size_t i = 0; i < reg.size(); ++i) {
        size_t gi = names.size();
        for (size_t g = 0; g < names.size(); ++g) {
            if (names[g] == reg[i].first) {
                gi = g;
                break;
            }
        }
        if (gi == names.size()) {
            names.push_back(reg[i].first);
            groups.emplace_back();
        }
        groups[gi].push_back(reg[i].second);
    }

    int err = SMESH_SUCCESS;
    for (size_t g = 0; err == SMESH_SUCCESS && g < names.size(); ++g) {
        const int multi = groups[g].size() > 1;
        for (size_t k = 0; err == SMESH_SUCCESS && k < groups[g].size(); ++k) {
            const auto &ss = groups[g][k];
            if (!ss) {
                continue;
            }
            auto block = mesh.block(ss->block_id());
            if (!block) {
                err = SMESH_FAILURE;
                break;
            }
            const enum ElemType et = block->element_type();
            if (is_semistructured_type(et)) {
                if (mesh.comm()->rank() == 0) {
                    fprintf(stderr,
                            "write_with_xdmf: skipping SS sideset '%s' surface grid (%s)\n",
                            names[g].c_str(),
                            type_to_string(et));
                }
                continue;
            }

            std::string leaf = std::string("sidesets/") + names[g];
            if (multi) {
                leaf += "/" + std::to_string((long long)ss->block_id());
            }
            std::string grid = names[g];
            if (multi) {
                grid += "_" + std::to_string((long long)ss->block_id());
            }

            const ptrdiff_t      n_ss   = ss->size();
            const element_idx_t *parent = n_ss > 0 ? ss->parent()->data() : nullptr;
            const i16           *lfi    = n_ss > 0 ? ss->lfi()->data() : nullptr;
            idx_t *const        *elems  = block->elements()->data();

            auto emit = [&](const element_idx_t *p_use,
                            const i16           *l_use,
                            const ptrdiff_t      n_part,
                            const std::string   &stem,
                            const std::string   &gname) -> int {
                enum ElemType face_st = INVALID;
                if (n_part > 0) {
                    face_st = shell_type(side_type(et, l_use[0]));
                } else if (elem_sides_homogeneous(et)) {
                    face_st = shell_type(side_type(et));
                } else {
                    face_st = (stem.find("tri") != std::string::npos) ? TRISHELL3 : QUADSHELL4;
                }
                const char *topo = xdmf_topology_type(face_st);
                if (!topo) {
                    return SMESH_SUCCESS;
                }
                const int         nnxs = elem_num_nodes(face_st);
                const std::string rel =
                        leaf + "/" + stem + "." + std::string(TypeToString<idx_t>::value());
                return write_xdmf_surface_part(
                        mesh, root, gname, rel, et, elems, p_use, l_use, n_part, nnxs, topo, surfs);
            };

            if (elem_sides_homogeneous(et)) {
                err = emit(parent, lfi, n_ss, "surface", grid);
            } else {
                ptrdiff_t n_tri  = 0;
                ptrdiff_t n_quad = 0;
                for (ptrdiff_t i = 0; i < n_ss; ++i) {
                    if (side_type(et, lfi[i]) == TRI3) {
                        ++n_tri;
                    } else {
                        ++n_quad;
                    }
                }
                SharedBuffer<element_idx_t> p_tri, p_quad;
                SharedBuffer<i16>           l_tri, l_quad;
                if (n_tri > 0) {
                    p_tri           = create_host_buffer<element_idx_t>((size_t)n_tri);
                    l_tri           = create_host_buffer<i16>((size_t)n_tri);
                    ptrdiff_t o     = 0;
                    element_idx_t *p_tri_d = p_tri->data();
                    i16           *l_tri_d = l_tri->data();
                    for (ptrdiff_t i = 0; i < n_ss; ++i) {
                        if (side_type(et, lfi[i]) == TRI3) {
                            p_tri_d[o] = parent[i];
                            l_tri_d[o] = lfi[i];
                            ++o;
                        }
                    }
                }
                if (n_quad > 0) {
                    p_quad          = create_host_buffer<element_idx_t>((size_t)n_quad);
                    l_quad          = create_host_buffer<i16>((size_t)n_quad);
                    ptrdiff_t o     = 0;
                    element_idx_t *p_quad_d = p_quad->data();
                    i16           *l_quad_d = l_quad->data();
                    for (ptrdiff_t i = 0; i < n_ss; ++i) {
                        if (side_type(et, lfi[i]) != TRI3) {
                            p_quad_d[o] = parent[i];
                            l_quad_d[o] = lfi[i];
                            ++o;
                        }
                    }
                }
                err = emit(n_tri > 0 ? p_tri->data() : nullptr,
                           n_tri > 0 ? l_tri->data() : nullptr,
                           n_tri,
                           "surface_tri",
                           grid + "_tri");
                if (err == SMESH_SUCCESS) {
                    err = emit(n_quad > 0 ? p_quad->data() : nullptr,
                               n_quad > 0 ? l_quad->data() : nullptr,
                               n_quad,
                               "surface_quad",
                               grid + "_quad");
                }
            }
        }
    }
    return err;
}

int Mesh::write_with_xdmf(const Path &path) const {
    SMESH_TRACE_SCOPE("Mesh::write_with_xdmf");

    if (write(path) != SMESH_SUCCESS) {
        return SMESH_FAILURE;
    }

    const int n_blocks = (int)this->n_blocks();
    if (n_blocks <= 0) {
        return SMESH_FAILURE;
    }

    auto             comm = this->comm();
    std::vector<i64> n_elem(n_blocks, 0);
    for (int b = 0; b < n_blocks; ++b) {
        auto block = this->block(b);
        if (!block) {
            continue;
        }
        if (comm->size() == 1) {
            n_elem[b] = (i64)block->n_elements();
        } else {
            n_elem[b] = comm->sum((i64)block->n_elements_owned());
        }
    }

    i64 n_nodes = 0;
    if (comm->size() == 1) {
        n_nodes = (i64)this->n_nodes();
    } else if (this->distributed()) {
        n_nodes = (i64)this->distributed()->n_nodes_global();
    } else {
        n_nodes = (i64)this->n_nodes();
    }

    int                    err = SMESH_SUCCESS;
    std::vector<ptrdiff_t> n_elem_pd(n_blocks);
    for (int b = 0; b < n_blocks; ++b) {
        n_elem_pd[b] = (ptrdiff_t)n_elem[b];
    }
    if (comm->rank() == 0) {
        const int from_memory = comm->size() == 1;
        if (n_blocks == 1) {
            err = write_block_connectivity_aos(*this, path, 0, n_elem_pd[0], from_memory);
        } else {
            for (int b = 0; err == SMESH_SUCCESS && b < n_blocks; ++b) {
                auto block = this->block(b);
                if (!block) {
                    err = SMESH_FAILURE;
                    break;
                }
                err = write_block_connectivity_aos(
                        *this, path / "blocks" / block->name(), (block_idx_t)b, n_elem_pd[b], from_memory);
            }
        }
    }
    comm->broadcast(&err, 1, 0);
    if (err != SMESH_SUCCESS) {
        return err;
    }

    std::vector<XdmfSurfGrid> surfs;
    err = write_xdmf_sideset_surfaces(*this, path, &surfs);
    if (err != SMESH_SUCCESS) {
        comm->broadcast(&err, 1, 0);
        return err;
    }
    err = write_xdmf_edgeset_grids(*this, path, &surfs);
    if (err != SMESH_SUCCESS) {
        comm->broadcast(&err, 1, 0);
        return err;
    }
    err = write_xdmf_nodeset_grids(*this, path, &surfs);
    if (err != SMESH_SUCCESS) {
        comm->broadcast(&err, 1, 0);
        return err;
    }
    if (comm->rank() == 0) {
        err = write_mesh_xdmf(*this, path, n_elem_pd.data(), (ptrdiff_t)n_nodes, surfs);
    }
    comm->broadcast(&err, 1, 0);
    return err;
}

}  // namespace smesh
