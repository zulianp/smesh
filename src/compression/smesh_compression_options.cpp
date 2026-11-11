#include "smesh_compression_options.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

namespace smesh {
namespace compression {

namespace {
bool parse_double(std::string_view s, double *v) {
    const std::string t(s);
    char             *end = nullptr;
    *v                    = std::strtod(t.c_str(), &end);
    return !t.empty() && end && *end == '\0' && std::isfinite(*v);
}
bool parse_int(std::string_view s, long long *v) {
    const std::string t(s);
    char             *end = nullptr;
    *v                    = std::strtoll(t.c_str(), &end, 10);
    return !t.empty() && end && *end == '\0';
}
}  // namespace

CompressionOptions CompressionOptions::lossless() {
    CompressionOptions o;
    o.enabled = true;
    return o;
}

int CompressionOptions::parse(std::string_view spec, CompressionOptions *out) {
    CompressionOptions o;
    size_t             start = 0;
    bool               first = true;
    while (start <= spec.size()) {
        size_t comma = spec.find(',', start);
        if (comma == std::string_view::npos) comma = spec.size();
        const std::string_view tok   = spec.substr(start, comma - start);
        const size_t           colon = tok.find(':');
        const std::string_view key   = tok.substr(0, colon);
        const std::string_view val   = colon == std::string_view::npos ? std::string_view() : tok.substr(colon + 1);
        start                        = comma + 1;
        if (tok.empty()) {
            if (comma == spec.size()) break;
            continue;
        }
        if (first) {
            first = false;
            if (key == "off" || key == "none" || key == "0") {
                o.enabled = false;
                continue;
            }
            const BoundMode modes[] = {BoundMode::Lossless, BoundMode::Abs, BoundMode::RelRange, BoundMode::PointwiseRel,
                                       BoundMode::NodeHeight};
            bool            matched = false;
            for (BoundMode m : modes) {
                if (key == compression::to_string(m) || (m == BoundMode::Lossless && (key == "on" || key == "1"))) {
                    o.enabled = true;
                    o.mode    = m;
                    matched   = true;
                    if (m != BoundMode::Lossless) {
                        if (!parse_double(val, &o.eb) || !(o.eb > 0)) {
                            fprintf(stderr, "compression: '%.*s' needs a positive bound, e.g. %s:1e-3\n", int(tok.size()),
                                    tok.data(), compression::to_string(m));
                            return SMESH_FAILURE;
                        }
                        if (m == BoundMode::PointwiseRel && !(o.eb < 1)) {
                            fprintf(stderr, "compression: pwrel needs a bound below 1\n");
                            return SMESH_FAILURE;
                        }
                        if (m == BoundMode::NodeHeight && !(o.eb < 0.25)) {
                            fprintf(stderr, "compression: h needs a bound below 0.25 (elements could invert)\n");
                            return SMESH_FAILURE;
                        }
                    }
                }
            }
            if (matched) continue;
            o.enabled = true;  // options without a mode: lossless
        }
        long long iv = 0;
        if (key == "zstd") {
            o.byte_stage = ByteStage::Zstd;
            if (!val.empty()) {
                if (!parse_int(val, &iv) || iv < -7 || iv > 22) {
                    fprintf(stderr, "compression: bad zstd level '%.*s'\n", int(val.size()), val.data());
                    return SMESH_FAILURE;
                }
                o.zstd_level = int(iv);
            }
        } else if (key == "lz4") {
            o.byte_stage = ByteStage::Lz4;
        } else if (key == "shuffle") {
            o.shuffle = true;
        } else if (key == "chunk") {
            if (!parse_int(val, &iv) || iv < 1 || iv > (i64(1) << 28)) {
                fprintf(stderr, "compression: bad chunk size '%.*s' (1 .. 2^28 values)\n", int(val.size()), val.data());
                return SMESH_FAILURE;
            }
            o.chunk_values = iv;
        } else if (key == "nochecksum") {
            o.checksums = false;
        } else if (key == "verify") {
            o.verify_write = true;
        } else {
            fprintf(stderr, "compression: unknown option '%.*s' in '%.*s'\n", int(tok.size()), tok.data(), int(spec.size()),
                    spec.data());
            return SMESH_FAILURE;
        }
    }
    *out = o;
    return SMESH_SUCCESS;
}

CompressionOptions CompressionOptions::from_env(const char *var, const CompressionOptions &dflt) {
    CompressionOptions o   = dflt;
    const char        *val = var ? std::getenv(var) : nullptr;
    if (val && *val) {
        if (parse(val, &o) != SMESH_SUCCESS) {
            SMESH_ERROR("invalid %s='%s'\n", var, val);
        }
    }
    if (const char *c = std::getenv("SMESH_COMPRESS_CHUNK"); c && *c) {
        long long iv = 0;
        if (!parse_int(c, &iv) || iv < 1 || iv > (i64(1) << 28)) {
            SMESH_ERROR("invalid SMESH_COMPRESS_CHUNK='%s'\n", c);
        }
        o.chunk_values = iv;
    }
    if (const char *v = std::getenv("SMESH_COMPRESS_VERIFY_WRITE"); v && *v) {
        o.verify_write = std::atoi(v) != 0;
    }
    return o;
}

CompressionOptions CompressionOptions::from_env(const char *var) { return from_env(var, CompressionOptions()); }

MeshWriteOptions MeshWriteOptions::all(const CompressionOptions &o) {
    MeshWriteOptions m;
    m.connectivity = o;
    m.coordinates  = o;
    m.sets         = o;
    if (o.is_lossy()) {  // a bound only makes sense for coordinates
        m.connectivity.mode = BoundMode::Lossless;
        m.sets.mode         = BoundMode::Lossless;
    }
    return m;
}

MeshWriteOptions MeshWriteOptions::from_env() {
    const CompressionOptions base = CompressionOptions::from_env("SMESH_COMPRESS");
    MeshWriteOptions         m    = all(base);
    m.connectivity                = CompressionOptions::from_env("SMESH_COMPRESS_CONNECTIVITY", m.connectivity);
    m.coordinates                 = CompressionOptions::from_env("SMESH_COMPRESS_COORDS", m.coordinates);
    m.sets                        = CompressionOptions::from_env("SMESH_COMPRESS_SETS", m.sets);
    if (m.validate() != SMESH_SUCCESS) {
        SMESH_ERROR("invalid SMESH_COMPRESS* settings\n");
    }
    return m;
}

int MeshWriteOptions::validate() const {
    if (connectivity.is_lossy() || sets.is_lossy()) {
        fprintf(stderr, "compression: connectivity and sets can only be compressed losslessly\n");
        return SMESH_FAILURE;
    }
    return SMESH_SUCCESS;
}

std::string CompressionOptions::to_string() const {
    if (!enabled) return "off";
    std::string s = compression::to_string(mode);
    if (mode != BoundMode::Lossless) {
        char buf[32];
        snprintf(buf, sizeof(buf), ":%g", eb);
        s += buf;
    }
    if (byte_stage == ByteStage::Zstd) s += ",zstd:" + std::to_string(zstd_level);
    if (byte_stage == ByteStage::Lz4) s += ",lz4";
    if (shuffle) s += ",shuffle";
    if (chunk_values != (i64(1) << 16)) s += ",chunk:" + std::to_string(chunk_values);
    if (!checksums) s += ",nochecksum";
    if (verify_write) s += ",verify";
    return s;
}

}  // namespace compression
}  // namespace smesh
