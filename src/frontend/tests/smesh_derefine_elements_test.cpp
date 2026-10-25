#include <memory>
#include <vector>

#include "smesh_buffer.hpp"
#include "smesh_semistructured.hpp"
#include "smesh_ssquad4.hpp"
#include "smesh_sstet4.hpp"
#include "smesh_test.hpp"

using namespace smesh;

// Fine row r is &storage[r], so a derefined row tells which fine row it points at.
struct Rows {
    std::vector<idx_t>   storage;
    std::vector<idx_t *> rows;
    explicit Rows(const int n) : storage(static_cast<size_t>(n)), rows(static_cast<size_t>(n)) {
        for (int r = 0; r < n; ++r) rows[static_cast<size_t>(r)] = &storage[static_cast<size_t>(r)];
    }
    int index(const idx_t *const p) const { return static_cast<int>(p - storage.data()); }
};

static const enum ElemType families[] = {HEX8, TET4, WEDGE6, PYRAMID5, QUAD4, TRI3, EDGE2};

// Local vertex rows of a level-L element, in the family's vertex order.
static int corners(const enum ElemType family, const int L, int *const c) {
    int n = 0;
    if (family == TRI3) {
        c[0] = sstri_lidx(L, 0, 0);
        c[1] = sstri_lidx(L, L, 0);
        c[2] = sstri_lidx(L, 0, L);
        return 3;
    }
    if (family == EDGE2) {
        c[0] = ssedge_lidx(L, 0);
        c[1] = ssedge_lidx(L, L);
        return 2;
    }
    return ss_source_family_corners(family, L, c, &n) ? n : 0;
}

static std::vector<int> derefined(const enum ElemType family, const int from, const int to) {
    Rows                 fine(ss_family_nxe(family, from));
    std::vector<idx_t *> coarse(static_cast<size_t>(ss_family_nxe(family, to)), nullptr);
    if (derefine_element_rows(family, from, to, fine.rows.data(), coarse.data()) != SMESH_SUCCESS) return {};
    std::vector<int> ret;
    for (auto *p : coarse) ret.push_back(p ? fine.index(p) : -1);
    return ret;
}

// At level 1 only the vertices remain: each coarse vertex row is the fine element's vertex row.
static int test_level_one_keeps_the_vertices() {
    for (const auto family : families) {
        int       c1[8], cL[8];
        const int n1 = corners(family, 1, c1);
        SMESH_TEST_ASSERT(n1 > 0);
        SMESH_TEST_EQ(ss_family_nxe(family, 1), n1);
        for (const int from : {2, 3, 4, 6}) {
            const int nL = corners(family, from, cL);
            SMESH_TEST_EQ(nL, n1);
            const auto rows = derefined(family, from, 1);
            SMESH_TEST_EQ((int)rows.size(), n1);
            for (int v = 0; v < n1; ++v) {
                SMESH_TEST_EQ(rows[static_cast<size_t>(c1[v])], cL[v]);
            }
        }
    }
    return SMESH_TEST_SUCCESS;
}

// Derefining in two steps lands on the same fine rows as in one, and no fine row is taken twice.
static int test_steps_compose() {
    const int levels[][3] = {{4, 2, 1}, {6, 3, 1}, {8, 4, 2}, {6, 2, 1}};
    for (const auto family : families) {
        for (const auto &l : levels) {
            const auto direct = derefined(family, l[0], l[2]);
            const auto first  = derefined(family, l[0], l[1]);
            const auto second = derefined(family, l[1], l[2]);
            SMESH_TEST_ASSERT(!direct.empty() && !first.empty() && !second.empty());
            SMESH_TEST_EQ(direct.size(), second.size());
            std::vector<unsigned char> seen(static_cast<size_t>(ss_family_nxe(family, l[0])), 0);
            for (size_t k = 0; k < direct.size(); ++k) {
                SMESH_TEST_EQ(direct[k], first[static_cast<size_t>(second[k])]);
                SMESH_TEST_ASSERT(seen[static_cast<size_t>(direct[k])] == 0);
                seen[static_cast<size_t>(direct[k])] = 1;
            }
        }
    }
    return SMESH_TEST_SUCCESS;
}

static int test_rejects_invalid_input() {
    for (const auto family : families) {
        SMESH_TEST_ASSERT(derefined(family, 4, 3).empty());
        SMESH_TEST_ASSERT(derefined(family, 2, 4).empty());
        SMESH_TEST_ASSERT(derefined(family, 4, 0).empty());
    }
    SMESH_TEST_EQ(ss_family_nxe(TET10, 2), 0);
    Rows                 fine(4);
    std::vector<idx_t *> coarse(4);
    SMESH_TEST_ASSERT(derefine_element_rows(TET10, 2, 1, fine.rows.data(), coarse.data()) != SMESH_SUCCESS);
    return SMESH_TEST_SUCCESS;
}

// The buffer view: one column per element, rows into the fine connectivity, which it keeps alive.
static int test_view_keeps_the_fine_connectivity() {
    const int       from = 4, to = 2, nxe = ss_family_nxe(QUAD4, from);
    const ptrdiff_t ne   = 3;
    auto            fine = create_host_buffer<idx_t>(nxe, ne);
    for (int r = 0; r < nxe; ++r) {
        for (ptrdiff_t e = 0; e < ne; ++e) fine->data()[r][e] = static_cast<idx_t>(r * ne + e);
    }
    auto view = derefine_elements(QUAD4, from, to, fine);
    SMESH_TEST_ASSERT(view != nullptr);
    SMESH_TEST_EQ((int)view->extent(0), ss_family_nxe(QUAD4, to));
    SMESH_TEST_EQ((ptrdiff_t)view->extent(1), ne);
    const auto           rows = derefined(QUAD4, from, to);
    std::vector<idx_t *> fine_rows(fine->data(), fine->data() + nxe);
    fine.reset();
    for (size_t k = 0; k < rows.size(); ++k) {
        SMESH_TEST_ASSERT(view->data()[k] == fine_rows[static_cast<size_t>(rows[k])]);
        SMESH_TEST_EQ(view->data()[k][ne - 1], static_cast<idx_t>(rows[k] * ne + ne - 1));
    }
    SMESH_TEST_ASSERT(derefine_elements(QUAD4, to, from, view) == nullptr);
    return SMESH_TEST_SUCCESS;
}

int main(int argc, char *argv[]) {
    SMESH_UNIT_TEST_INIT(argc, argv);
    SMESH_RUN_TEST(test_level_one_keeps_the_vertices);
    SMESH_RUN_TEST(test_steps_compose);
    SMESH_RUN_TEST(test_rejects_invalid_input);
    SMESH_RUN_TEST(test_view_keeps_the_fine_connectivity);
    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}
