#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <queue>
#include <sstream>

using namespace Rcpp;

// ── Interval / Sequence types ─────────────────────────────────────────────────

struct Interval {
    int32_t start;
    int32_t end;
    int8_t  direction; // +1 ascending, -1 descending, 0 constant (single point)
};

typedef std::vector<Interval> Sequence;

// ── R ↔ C++ conversion helpers ───────────────────────────────────────────────

static Sequence rlist_to_sequence(const List& r_list) {
    Sequence seq;
    int n = r_list.size();
    seq.reserve(n);
    for (int i = 0; i < n; ++i) {
        List iv = r_list[i];
        Interval interval;
        interval.start     = as<int>(iv["start"]);
        interval.end       = as<int>(iv["end"]);
        interval.direction = (int8_t) as<int>(iv["direction"]);
        seq.push_back(interval);
    }
    return seq;
}

static List sequence_to_rlist(const Sequence& seq) {
    List result(seq.size());
    for (size_t i = 0; i < seq.size(); ++i) {
        result[i] = List::create(
            Named("start")     = (int) seq[i].start,
            Named("end")       = (int) seq[i].end,
            Named("direction") = (int) seq[i].direction
        );
    }
    return result;
}

// ── Internal sequence operations ──────────────────────────────────────────────

static int seq_length_impl(const Sequence& seq) {
    int L = 0;
    for (size_t i = 0; i < seq.size(); ++i) {
        L += (seq[i].direction == 0) ? 1 : (std::abs(seq[i].end - seq[i].start) + 1);
    }
    return L;
}

static Sequence reverse_sequence_impl(const Sequence& seq) {
    int n = seq.size();
    Sequence result(n);
    for (int i = 0; i < n; ++i) {
        const Interval& iv = seq[n - 1 - i];
        Interval r;
        r.start     = iv.end;
        r.end       = iv.start;
        r.direction = (int8_t)(-iv.direction);
        result[i]   = r;
    }
    return result;
}

static Sequence fuse_sequence_impl(const Sequence& seq) {
    Sequence fused = seq;
    Sequence rev   = reverse_sequence_impl(seq);
    fused.insert(fused.end(), rev.begin(), rev.end());
    return fused;
}

// Returns {left_seq, right_seq} split at cut_index (1-based count of elements).
// Elements 1..cut_index go to left; elements cut_index+1..end go to right.
static void cut_sequence_impl(const Sequence& seq, int cut_index,
                               Sequence& left_seq, Sequence& right_seq) {
    left_seq.clear();
    right_seq.clear();
    int current_length = 0;

    for (size_t i = 0; i < seq.size(); ++i) {
        const Interval& iv = seq[i];
        int iv_len = (iv.direction == 0) ? 1 : (std::abs(iv.end - iv.start) + 1);

        if (current_length >= cut_index) {
            right_seq.push_back(iv);
        } else if (current_length + iv_len <= cut_index) {
            left_seq.push_back(iv);
        } else {
            // Straddles the cut
            int cut_within = cut_index - current_length;
            if (iv.direction == 1) {
                Interval l_iv, r_iv;
                l_iv.start = iv.start; l_iv.end = iv.start + cut_within - 1; l_iv.direction = 1;
                r_iv.start = iv.start + cut_within; r_iv.end = iv.end; r_iv.direction = 1;
                left_seq.push_back(l_iv);
                right_seq.push_back(r_iv);
            } else if (iv.direction == -1) {
                Interval l_iv, r_iv;
                l_iv.start = iv.start; l_iv.end = iv.start - cut_within + 1; l_iv.direction = -1;
                r_iv.start = iv.start - cut_within; r_iv.end = iv.end; r_iv.direction = -1;
                left_seq.push_back(l_iv);
                right_seq.push_back(r_iv);
            } else {
                // direction == 0: single-element — unreachable in practice (iv_len=1
                // makes current_length+1 never straddle an integer cut_index)
                left_seq.push_back(iv);
                right_seq.push_back(iv);
            }
        }
        current_length += iv_len;
    }
}

static int hotspot_copies_impl(const Sequence& seq, int bin) {
    int count = 0;
    for (size_t i = 0; i < seq.size(); ++i) {
        int lo = std::min((int)seq[i].start, (int)seq[i].end);
        int hi = std::max((int)seq[i].start, (int)seq[i].end);
        if (bin >= lo && bin <= hi) ++count;
    }
    return count;
}

static IntegerVector seq2vec_impl(const Sequence& seq) {
    int total = seq_length_impl(seq);
    IntegerVector result(total);
    int pos = 0;
    for (size_t i = 0; i < seq.size(); ++i) {
        const Interval& iv = seq[i];
        if (iv.direction == 0) {
            result[pos++] = iv.start;
        } else if (iv.direction == 1) {
            for (int v = iv.start; v <= iv.end; ++v) result[pos++] = v;
        } else { // direction == -1
            for (int v = iv.start; v >= iv.end; --v) result[pos++] = v;
        }
    }
    return result;
}

static Sequence vec2seq_impl(const IntegerVector& vec) {
    int n = vec.size();
    if (n == 0) return Sequence();
    if (n == 1) {
        Interval iv;
        iv.start = vec[0]; iv.end = vec[0]; iv.direction = 0;
        return Sequence(1, iv);
    }

    Sequence intervals;
    intervals.reserve(n); // worst-case pre-allocation

    int32_t seg_start = vec[0];
    int32_t seg_end   = vec[0];
    int8_t  seg_dir   = 0;
    bool    dir_set   = false;

    for (int i = 1; i < n; ++i) {
        int step = vec[i] - vec[i - 1];

        if (step == 1 || step == -1) {
            int8_t step_dir = (step == 1) ? (int8_t)1 : (int8_t)-1;
            if (!dir_set) {
                seg_dir = step_dir;
                seg_end = vec[i];
                dir_set = true;
            } else if (step_dir == seg_dir) {
                seg_end = vec[i];
            } else {
                // Direction reversal — close current interval, start fresh
                Interval iv;
                iv.start = seg_start; iv.end = seg_end; iv.direction = seg_dir;
                intervals.push_back(iv);
                seg_start = vec[i];
                seg_end   = vec[i];
                dir_set   = false;
                seg_dir   = 0;
            }
        } else {
            // Non-unit step (gap) — close current interval, start fresh
            Interval iv;
            iv.start     = seg_start;
            iv.end       = seg_end;
            iv.direction = dir_set ? seg_dir : (int8_t)0;
            intervals.push_back(iv);
            seg_start = vec[i];
            seg_end   = vec[i];
            dir_set   = false;
            seg_dir   = 0;
        }
    }
    // Close final interval
    Interval iv;
    iv.start     = seg_start;
    iv.end       = seg_end;
    iv.direction = dir_set ? seg_dir : (int8_t)0;
    intervals.push_back(iv);

    return intervals;
}

// is_dup=true → duplication; is_dup=false → deletion.
// rate is the mean of the exponential event-length distribution.
static Sequence sim_amp_del_impl(const Sequence& seq, bool is_dup, double rate) {
    int n_iv = seq.size();
    if (n_iv == 0) return seq;

    // Compute per-interval lengths
    std::vector<int> iv_lengths(n_iv);
    int total_len = 0;
    for (int j = 0; j < n_iv; ++j) {
        int len = (seq[j].direction == 0) ? 1 : (std::abs(seq[j].end - seq[j].start) + 1);
        iv_lengths[j] = len;
        total_len += len;
    }
    if (total_len == 0) return seq;

    // Eligible intervals: prefer length >= 2; fall back to all
    std::vector<int> eligible;
    eligible.reserve(n_iv);
    for (int j = 0; j < n_iv; ++j) {
        if (iv_lengths[j] >= 2) eligible.push_back(j);
    }
    if (eligible.empty()) {
        for (int j = 0; j < n_iv; ++j) eligible.push_back(j);
    }

    // Sample one interval uniformly
    int pick = (int)(R::unif_rand() * (double)eligible.size());
    if (pick >= (int)eligible.size()) pick = (int)eligible.size() - 1;
    int idx    = eligible[pick];
    Interval iv = seq[idx];
    int iv_len  = iv_lengths[idx];

    // Sample event length: rexp with mean = rate, rounded, clamped to [1, iv_len]
    int event_length  = std::max(1, (int)std::round(R::rexp(rate)));
    int actual_length = std::min(event_length, iv_len);

    // Sample start offset in [0, iv_len - actual_length]
    int max_offset = iv_len - actual_length;
    int offset = (max_offset <= 0) ? 0
                 : std::min((int)(R::unif_rand() * (double)(max_offset + 1)), max_offset);

    // Derive sub-segment and surrounding fragments
    Interval before_iv, after_iv, seg;
    bool has_before = false, has_after = false;

    if (iv.direction == 1) {
        int seg_s = iv.start + offset;
        int seg_e = seg_s + actual_length - 1;
        seg.start = seg_s; seg.end = seg_e; seg.direction = 1;
        if (offset > 0) {
            before_iv.start = iv.start; before_iv.end = seg_s - 1; before_iv.direction = 1;
            has_before = true;
        }
        if (seg_e < iv.end) {
            after_iv.start = seg_e + 1; after_iv.end = iv.end; after_iv.direction = 1;
            has_after = true;
        }
    } else if (iv.direction == -1) {
        int seg_s = iv.start - offset;
        int seg_e = seg_s - actual_length + 1;
        seg.start = seg_s; seg.end = seg_e; seg.direction = -1;
        if (offset > 0) {
            before_iv.start = iv.start; before_iv.end = seg_s + 1; before_iv.direction = -1;
            has_before = true;
        }
        if (seg_e > iv.end) {
            after_iv.start = seg_e - 1; after_iv.end = iv.end; after_iv.direction = -1;
            has_after = true;
        }
    } else {
        // direction == 0
        seg = iv;
    }

    // For deletion: guard against producing an empty sequence
    if (!is_dup && idx == 0 && idx == n_iv - 1 && !has_before && !has_after) {
        return seq;
    }

    Sequence result;
    result.reserve(n_iv + 3);

    // Head (intervals before idx)
    for (int j = 0; j < idx; ++j) result.push_back(seq[j]);

    // Mid
    if (has_before) result.push_back(before_iv);
    if (is_dup) {
        result.push_back(seg);
        result.push_back(seg);
    }
    // del: seg is dropped
    if (has_after) result.push_back(after_iv);

    // Tail (intervals after idx)
    for (int j = idx + 1; j < n_iv; ++j) result.push_back(seq[j]);

    return result;
}

// WGD: concatenate sequence with itself at the interval level — O(n_intervals)
static Sequence sim_wgd_impl(const Sequence& seq) {
    Sequence result = seq;
    result.insert(result.end(), seq.begin(), seq.end());
    return result;
}

// ── Exported wrappers (Phase 1) ───────────────────────────────────────────────

//' Compute length of an interval-encoded sequence
//'
//' @param seq_list R list of intervals (each a named list with start, end, direction)
//' @return Integer total length
//' @export
// [[Rcpp::export]]
int seq_length_cpp(List seq_list) {
    return seq_length_impl(rlist_to_sequence(seq_list));
}

//' Reverse an interval-encoded sequence
//'
//' @param seq_list R list of intervals
//' @return R list of reversed intervals
//' @export
// [[Rcpp::export]]
List reverse_sequence_cpp(List seq_list) {
    return sequence_to_rlist(reverse_sequence_impl(rlist_to_sequence(seq_list)));
}

//' Fuse a sequence with its reverse (BFB palindrome creation)
//'
//' @param seq_list R list of intervals
//' @return R list of fused intervals
//' @export
// [[Rcpp::export]]
List fuse_sequence_cpp(List seq_list) {
    return sequence_to_rlist(fuse_sequence_impl(rlist_to_sequence(seq_list)));
}

//' Cut an interval-encoded sequence at a given index
//'
//' Elements 1..cut_index go to left; elements cut_index+1..end go to right.
//'
//' @param seq_list R list of intervals
//' @param cut_index Integer cut position (1-based)
//' @return Named list with \code{left_seq} and \code{right_seq}
//' @export
// [[Rcpp::export]]
List cut_sequence_cpp(List seq_list, int cut_index) {
    Sequence seq = rlist_to_sequence(seq_list);
    Sequence left_seq, right_seq;
    cut_sequence_impl(seq, cut_index, left_seq, right_seq);
    return List::create(
        Named("left_seq")  = sequence_to_rlist(left_seq),
        Named("right_seq") = sequence_to_rlist(right_seq)
    );
}

//' Count copies of a hotspot bin in an interval-encoded sequence
//'
//' @param seq_list R list of intervals
//' @param bin Integer bin position to query
//' @return Integer copy count
//' @export
// [[Rcpp::export]]
int hotspot_copies_cpp(List seq_list, int bin) {
    return hotspot_copies_impl(rlist_to_sequence(seq_list), bin);
}

//' Expand interval-encoded sequence to an integer vector
//'
//' @param seq_list R list of intervals
//' @return IntegerVector of genomic bin values
//' @export
// [[Rcpp::export]]
IntegerVector seq2vec_cpp(List seq_list) {
    return seq2vec_impl(rlist_to_sequence(seq_list));
}

//' Compress an integer vector to an interval-encoded sequence
//'
//' @param vec IntegerVector of genomic bin values
//' @return R list of intervals
//' @export
// [[Rcpp::export]]
List vec2seq_cpp(IntegerVector vec) {
    return sequence_to_rlist(vec2seq_impl(vec));
}

//' Simulate amplification or deletion on an interval-encoded sequence
//'
//' @param seq_list R list of intervals
//' @param operation \code{"dup"} for duplication or \code{"del"} for deletion
//' @param rate Mean of the exponential event-length distribution
//' @return R list of intervals
//' @export
// [[Rcpp::export]]
List sim_amp_del_cpp(List seq_list, std::string operation, double rate) {
    RNGScope rng_scope;
    bool is_dup = (operation == "dup");
    if (!is_dup && operation != "del") Rcpp::stop("operation not recognised: use 'dup' or 'del'");
    return sequence_to_rlist(sim_amp_del_impl(rlist_to_sequence(seq_list), is_dup, rate));
}

//' Simulate whole-genome duplication on an interval-encoded sequence
//'
//' Concatenates the interval list with itself — O(n_intervals), no expand/compress.
//'
//' @param seq_list R list of intervals
//' @return R list of intervals representing the doubled genome
//' @export
// [[Rcpp::export]]
List sim_wgd_cpp(List seq_list) {
    return sequence_to_rlist(sim_wgd_impl(rlist_to_sequence(seq_list)));
}

//' Simulate BFB left and right daughter sequences
//'
//' Implements the fuse-cut-reverse BFB cycle with configurable breakpoint
//' selection (uniform or beta distribution).
//'
//' @param seq_list R list of intervals
//' @param support Breakpoint distribution: \code{"uniform"} or \code{"beta"}
//'   (\code{"custom"} falls back to the R implementation)
//' @param alpha Beta distribution shape parameter (ignored unless support="beta")
//' @param beta_param Beta distribution shape parameter (ignored unless support="beta")
//' @return Named list with \code{l_seq} and \code{r_seq} (each an R interval list)
//' @export
// [[Rcpp::export]]
List sim_bfb_cpp(
    List seq_list,
    std::string support    = "uniform",
    double alpha           = NA_REAL,
    double beta_param      = NA_REAL
) {
    RNGScope rng_scope;
    Sequence seq = rlist_to_sequence(seq_list);
    int L = seq_length_impl(seq);

    // Collect fusion values: junction bins where interval[k].end == interval[k+1].start.
    // These are existing BFB fold-back points; breakpoints landing there are forbidden.
    std::vector<int> bps;
    for (size_t k = 0; k + 1 < seq.size(); ++k) {
        if (seq[k].end == seq[k + 1].start) bps.push_back(seq[k].end);
    }

    auto is_forbidden = [&](int bp) -> bool {
        if (bp == L) return true;
        for (size_t i = 0; i < bps.size(); ++i) if (bp == bps[i]) return true;
        return false;
    };

    int bp_idx      = L; // start forbidden so loop executes
    int attempts    = 0;
    const int MAX_A = 10;

    while (is_forbidden(bp_idx) && attempts < MAX_A) {
        ++attempts;
        if (support == "uniform") {
            // sample from 1..(2*L) — mirrors R: sample(1:(2*L), 1)
            bp_idx = (int)(R::unif_rand() * (double)(2 * L)) + 1;
            if (bp_idx > 2 * L) bp_idx = 2 * L;
        } else if (support == "beta") {
            if (ISNA(alpha) || ISNA(beta_param))
                Rcpp::stop("For beta distribution both alpha and beta must be provided");
            double tau = R::rbeta(alpha, beta_param);
            bp_idx = std::max(1, (int)std::round(tau * 2.0 * L));
        } else if (support == "custom") {
            // Custom breakpoints involve undefined R code in the original — defer to R
            Rcpp::stop("support='custom' is not implemented in the C++ path; "
                       "call sim_bfb_left_and_right_sequences() directly from R");
        } else {
            Rcpp::stop("Unsupported distribution type. Use 'uniform', 'beta', or 'custom'.");
        }
    }

    if (attempts >= MAX_A && is_forbidden(bp_idx)) {
        Rcpp::warning("BFB breakpoint selection failed after %d attempts "
                      "(sequence may be too fragmented). Returning unchanged sequence.",
                      MAX_A);
        return List::create(Named("l_seq") = seq_list, Named("r_seq") = seq_list);
    }

    // Fuse → cut → reverse right half
    Sequence fused = fuse_sequence_impl(seq);
    Sequence left_seq, right_seq;
    cut_sequence_impl(fused, bp_idx, left_seq, right_seq);
    Sequence r_rev = reverse_sequence_impl(right_seq);

    // Random 50/50 swap (mirrors R: if (runif(1) > .5) swap)
    bool do_swap = (R::unif_rand() > 0.5);
    return List::create(
        Named("l_seq") = sequence_to_rlist(do_swap ? r_rev   : left_seq),
        Named("r_seq") = sequence_to_rlist(do_swap ? left_seq : r_rev)
    );
}

// ── Phase 4: Full Gillespie loop in C++ ──────────────────────────────────────

// ── Phase 4 data structures ───────────────────────────────────────────────────

struct SimP {
    double birth_rate, death_rate;
    double pos_sel_rate, neg_sel_rate;
    double wgd_prob;
    int    max_cells;
    double max_time;
    int    hotspot_allele_idx;   // -1 = no hotspot
    int    hotspot_bin;
    std::string breakpoint_support;
    double bp_alpha, bp_beta;
    int    bfb_allele_idx;
    int    n_alleles;
    std::vector<std::string> allele_names;
    std::vector<double>      event_probs;  // normal, bfb, amp, del (sum to 1)
    std::vector<std::string> event_names;  // parallel labels
    double lambda, rate;
};

struct CellD {
    std::string id;
    std::string parent_id;
    std::vector<Sequence> alleles;
    double next_event_time;
    bool   hotspot_gained;
    bool   alive;
    int    generation;   // incremented on death — invalidates stale PQ entries
};

struct PQEntry4 {
    double time;
    int    pool_idx;
    int    generation;
    bool operator>(const PQEntry4& o) const { return time > o.time; }
};

typedef std::priority_queue<
    PQEntry4, std::vector<PQEntry4>, std::greater<PQEntry4>
> MinHeap4;

struct HistRec {
    std::string cell_id;
    std::string parent_id;
    bool bfb_event;
    bool wgd_event;
    std::string cn_event;    // e.g. "bfb", "none", "amp,del"
    std::string chr_allele;  // e.g. "1:A", "all", "" (= NA)
};

// ── Phase 4 helpers ───────────────────────────────────────────────────────────

static std::string join_str(const std::vector<std::string>& v, char sep) {
    std::string r;
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) r += sep;
        r += v[i];
    }
    return r;
}

// Sample one index from a normalised probability vector
static int sample_one_w(const std::vector<double>& probs) {
    double u = R::unif_rand();
    double cs = 0.0;
    for (int i = 0; i < (int)probs.size() - 1; ++i) {
        cs += probs[i];
        if (u < cs) return i;
    }
    return (int)probs.size() - 1;
}

// Internal BFB for Phase 4 (operates on Sequence directly)
static void sim_bfb_internal(const Sequence& seq, const SimP& p,
                              Sequence& l_out, Sequence& r_out) {
    int L = seq_length_impl(seq);
    std::vector<int> bps;
    for (size_t k = 0; k + 1 < seq.size(); ++k)
        if (seq[k].end == seq[k + 1].start) bps.push_back(seq[k].end);

    auto forbidden = [&](int bp) {
        if (bp == L) return true;
        for (int b : bps) if (bp == b) return true;
        return false;
    };

    int bp_idx = L, att = 0;
    while (forbidden(bp_idx) && att < 10) {
        ++att;
        if (p.breakpoint_support == "uniform") {
            bp_idx = (int)(R::unif_rand() * (double)(2 * L)) + 1;
            if (bp_idx > 2 * L) bp_idx = 2 * L;
        } else if (p.breakpoint_support == "beta") {
            double tau = R::rbeta(p.bp_alpha, p.bp_beta);
            bp_idx = std::max(1, (int)std::round(tau * 2.0 * L));
        } else {
            Rcpp::stop("support='custom' not supported in C++ sim loop; use support='uniform' or 'beta'");
        }
    }
    if (att >= 10 && forbidden(bp_idx)) {
        Rcpp::warning("BFB breakpoint selection failed; returning unchanged sequence.");
        l_out = seq; r_out = seq; return;
    }
    Sequence fused = fuse_sequence_impl(seq);
    cut_sequence_impl(fused, bp_idx, l_out, r_out);
    r_out = reverse_sequence_impl(r_out);
    if (R::unif_rand() > 0.5) std::swap(l_out, r_out);
}

// Apply a single regular event (amp/del/normal) to alleles in-place
static void apply_event(std::vector<Sequence>& alleles,
                        const std::string& ev,
                        const SimP& p,
                        std::vector<std::string>& ev_log,
                        std::vector<std::string>& chr_log) {
    if (ev == "normal") {
        ev_log.push_back("normal");
        // chr_allele is NULL/NA for normal — don't append to chr_log (mirrors R)
    } else if (ev == "amp" || ev == "del") {
        int ai = (int)(R::unif_rand() * (double)p.n_alleles);
        if (ai >= p.n_alleles) ai = p.n_alleles - 1;
        alleles[ai] = sim_amp_del_impl(alleles[ai], (ev == "amp"), p.rate);
        ev_log.push_back(ev);
        chr_log.push_back(p.allele_names[ai]);
    }
}

// ── Phase 4 birth event ───────────────────────────────────────────────────────

static void process_birth4(
    std::vector<CellD>& pool,
    MinHeap4&           pq,
    std::vector<HistRec>& hist,
    int&  alive_count,
    int&  next_id,
    int&  wgd_available,
    const SimP& p,
    int   parent_idx,
    double current_time
) {
    // Sample event counts for each daughter
    int n_left  = (int)R::rpois(p.lambda);
    int n_right = (int)R::rpois(p.lambda);
    int total_ev = n_left + n_right;

    // Copy parent alleles
    std::vector<Sequence> l_al = pool[parent_idx].alleles;
    std::vector<Sequence> r_al = pool[parent_idx].alleles;

    bool bfb_occ = false, wgd_occ = false;
    std::vector<std::string> l_ev, r_ev, l_chr, r_chr;

    bool wgd_will = (wgd_available > 0) && (R::unif_rand() < p.wgd_prob);

    std::vector<int> all_ev(total_ev);
    for (int i = 0; i < total_ev; ++i) all_ev[i] = sample_one_w(p.event_probs);

    bool special_done = false;

    if (total_ev > 0) {
        bool has_bfb = false;
        for (int e : all_ev) if (p.event_names[e] == "bfb") { has_bfb = true; break; }

        // Resolve BFB vs WGD conflict
        if (wgd_will && has_bfb) {
            if (R::unif_rand() < 0.5) has_bfb = false; else wgd_will = false;
        }

        if (has_bfb) {
            int bi = p.bfb_allele_idx;
            sim_bfb_internal(pool[parent_idx].alleles[bi], p, l_al[bi], r_al[bi]);
            bfb_occ = true; special_done = true;
            l_ev.push_back("bfb"); r_ev.push_back("bfb");
            l_chr.push_back(p.allele_names[bi]); r_chr.push_back(p.allele_names[bi]);
            // Remove BFB entries from all_ev
            std::vector<int> rem;
            for (int e : all_ev) if (p.event_names[e] != "bfb") rem.push_back(e);
            all_ev = rem;

        } else if (wgd_will) {
            for (int ai = 0; ai < p.n_alleles; ++ai) {
                l_al[ai] = sim_wgd_impl(pool[parent_idx].alleles[ai]);
                r_al[ai] = sim_wgd_impl(pool[parent_idx].alleles[ai]);
            }
            wgd_occ = true; special_done = true; --wgd_available;
            l_ev.push_back("wgd"); r_ev.push_back("wgd");
            l_chr.push_back("all"); r_chr.push_back("all");
        }

        if (special_done) {
            // Redistribute remaining events across daughters
            int nl2 = std::max(0, n_left  - 1);
            int nr2 = std::max(0, n_right - 1);
            int rem = (int)all_ev.size();

            if (rem > 0 && (nl2 + nr2) > 0) {
                int la = std::min(nl2, rem);
                int ra = std::min(nr2, rem - la);

                // Fisher-Yates partial shuffle to pick la events for left
                std::vector<int> idx(rem);
                for (int i = 0; i < rem; ++i) idx[i] = i;
                for (int i = 0; i < la; ++i) {
                    int j = i + (int)(R::unif_rand() * (double)(rem - i));
                    if (j >= rem) j = rem - 1;
                    std::swap(idx[i], idx[j]);
                }

                // Collect which original positions went to left
                std::vector<bool> used_by_left(rem, false);
                for (int i = 0; i < la; ++i) {
                    used_by_left[idx[i]] = true;
                    apply_event(l_al, p.event_names[all_ev[idx[i]]], p, l_ev, l_chr);
                }
                // Right gets the first ra unused events (in order)
                int rc = 0;
                for (int i = 0; i < rem && rc < ra; ++i) {
                    if (!used_by_left[i]) {
                        apply_event(r_al, p.event_names[all_ev[i]], p, r_ev, r_chr);
                        ++rc;
                    }
                }
            }
        } else {
            // No special event: first n_left go to left, next n_right to right
            int lc = std::min(n_left, total_ev);
            int rc = std::min(n_right, total_ev - lc);
            for (int i = 0;      i < lc;      ++i) apply_event(l_al, p.event_names[all_ev[i]], p, l_ev, l_chr);
            for (int i = lc; i < lc + rc; ++i) apply_event(r_al, p.event_names[all_ev[i]], p, r_ev, r_chr);
        }

    } else if (wgd_will) {
        for (int ai = 0; ai < p.n_alleles; ++ai) {
            l_al[ai] = sim_wgd_impl(pool[parent_idx].alleles[ai]);
            r_al[ai] = sim_wgd_impl(pool[parent_idx].alleles[ai]);
        }
        wgd_occ = true; --wgd_available;
        l_ev.push_back("wgd"); r_ev.push_back("wgd");
        l_chr.push_back("all"); r_chr.push_back("all");
    }

    // Daughter IDs
    std::string l_id = "cell_" + std::to_string(next_id++);
    std::string r_id = "cell_" + std::to_string(next_id++);

    // Hotspot and rates
    bool l_hs = (p.hotspot_allele_idx >= 0)
        && hotspot_copies_impl(l_al[p.hotspot_allele_idx], p.hotspot_bin) > 1;
    bool r_hs = (p.hotspot_allele_idx >= 0)
        && hotspot_copies_impl(r_al[p.hotspot_allele_idx], p.hotspot_bin) > 1;

    double l_cr = p.birth_rate * (1.0 + p.pos_sel_rate * l_hs)
                + p.death_rate * (1.0 + p.neg_sel_rate * l_hs);
    double r_cr = p.birth_rate * (1.0 + p.pos_sel_rate * r_hs)
                + p.death_rate * (1.0 + p.neg_sel_rate * r_hs);

    double l_t = current_time + R::rexp(1.0 / l_cr);
    double r_t = current_time + R::rexp(1.0 / r_cr);

    // Kill parent
    std::string par_id = pool[parent_idx].id;
    pool[parent_idx].alive = false;
    pool[parent_idx].generation++;
    pool[parent_idx].alleles.clear();

    // Add daughters to pool and heap
    int l_idx = (int)pool.size();
    {
        CellD c;
        c.id = l_id; c.parent_id = par_id;
        c.alleles = std::move(l_al);
        c.next_event_time = l_t; c.hotspot_gained = l_hs;
        c.alive = true; c.generation = 0;
        pool.push_back(std::move(c));
    }
    pq.push({l_t, l_idx, 0});

    int r_idx = (int)pool.size();
    {
        CellD c;
        c.id = r_id; c.parent_id = par_id;
        c.alleles = std::move(r_al);
        c.next_event_time = r_t; c.hotspot_gained = r_hs;
        c.alive = true; c.generation = 0;
        pool.push_back(std::move(c));
    }
    pq.push({r_t, r_idx, 0});

    // alive_count: -1 parent, +2 daughters = +1
    alive_count += 1;

    // History records
    std::string l_cn  = l_ev.empty()  ? "none" : join_str(l_ev,  ',');
    std::string r_cn  = r_ev.empty()  ? "none" : join_str(r_ev,  ',');
    std::string l_chr_str = l_chr.empty() ? "" : join_str(l_chr, ',');
    std::string r_chr_str = r_chr.empty() ? "" : join_str(r_chr, ',');

    hist.push_back({l_id, par_id, bfb_occ, wgd_occ, l_cn, l_chr_str});
    hist.push_back({r_id, par_id, bfb_occ, wgd_occ, r_cn, r_chr_str});
}

// ── Parameter extraction ─────────────────────────────────────────────────────

static SimP extract_simp(const Rcpp::List& sim_state, double lambda, double rate) {
    SimP p;
    List ip = sim_state["input_parameters"];

    p.birth_rate    = as<double>(ip["birth_rate"]);
    p.death_rate    = as<double>(ip["death_rate"]);
    p.pos_sel_rate  = as<double>(ip["positive_selection_rate"]);
    p.neg_sel_rate  = as<double>(ip["negative_selection_rate"]);
    p.wgd_prob      = as<double>(ip["wgd_probability"]);
    p.max_cells     = as<int>(ip["max_cells"]);
    p.max_time      = as<double>(ip["max_time"]);
    p.lambda        = lambda;
    p.rate          = rate;
    p.breakpoint_support = as<std::string>(ip["breakpoint_support"]);

    SEXP alpha_sexp = ip["alpha"];
    p.bp_alpha = (Rf_isNull(alpha_sexp) || ISNAN(Rf_asReal(alpha_sexp)))
                 ? NA_REAL : as<double>(alpha_sexp);
    SEXP beta_sexp = ip["beta"];
    p.bp_beta  = (Rf_isNull(beta_sexp)  || ISNAN(Rf_asReal(beta_sexp)))
                 ? NA_REAL : as<double>(beta_sexp);

    p.allele_names = as<std::vector<std::string>>(ip["chr_alleles"]);
    p.n_alleles    = (int)p.allele_names.size();

    std::string bfb_name = as<std::string>(ip["bfb_allele"]);
    p.bfb_allele_idx = 0;
    for (int i = 0; i < p.n_alleles; ++i)
        if (p.allele_names[i] == bfb_name) { p.bfb_allele_idx = i; break; }

    // Hotspot
    p.hotspot_allele_idx = -1;
    p.hotspot_bin = -1;
    SEXP hs_sexp = ip["hotspot"];
    if (!Rf_isNull(hs_sexp)) {
        List hs = ip["hotspot"];
        SEXP chr_s = hs["chr"]; SEXP pos_s = hs["pos"];
        if (!Rf_isNull(chr_s) && !Rf_isNull(pos_s)) {
            std::string hs_chr = as<std::string>(chr_s);
            p.hotspot_bin = as<int>(pos_s);
            for (int i = 0; i < p.n_alleles; ++i)
                if (p.allele_names[i] == hs_chr) { p.hotspot_allele_idx = i; break; }
        }
    }

    // Event probabilities
    List rates = ip["rates"];
    p.event_names = {"normal", "bfb", "amp", "del"};
    p.event_probs = {
        as<double>(rates["normal"]),
        as<double>(rates["bfb"]),
        as<double>(rates["amp"]),
        as<double>(rates["del"])
    };

    return p;
}

// ── Main exported function ────────────────────────────────────────────────────

//' Run the Gillespie simulation loop in C++
//'
//' Takes an already-initialised \code{sim_state} from R's
//' \code{initialize_simulation()} and runs the main loop, returning
//' the final state in the same list format expected by
//' \code{prepare_results()}.
//'
//' @param sim_state_r Named list returned by \code{initialize_simulation()}
//' @param lambda Poisson rate for genomic events per daughter
//' @param rate Mean of exponential event-length distribution for amp/del
//' @return Named list with same structure as \code{sim_state_r}
//' @export
// [[Rcpp::export]]
List bridge_sim_loop_cpp(List sim_state_r, double lambda, double rate) {
    RNGScope rng_scope;

    SimP p = extract_simp(sim_state_r, lambda, rate);
    List ip = sim_state_r["input_parameters"];
    int wgd_available = as<int>(ip["wgd_available"]);

    // ── Initialise C++ state from R ──────────────────────────────────────────
    std::vector<CellD> pool;
    MinHeap4 pq;
    std::vector<HistRec> hist;
    int alive_count = 0;

    // Load initial cells
    CharacterVector cell_ids_r    = sim_state_r["cell_ids"];
    List            cell_seqs_r   = sim_state_r["cell_sequences"];
    NumericVector   event_times_r = sim_state_r["cell_next_event_times"];
    LogicalVector   hs_status_r   = sim_state_r["hotspot_status"];

    pool.reserve(as<int>(ip["max_cells"]) * 4);
    for (int i = 0; i < (int)cell_ids_r.size(); ++i) {
        std::string cid = as<std::string>(cell_ids_r[i]);
        List seqs_for_cell = cell_seqs_r[cid];

        CellD c;
        c.id = cid; c.parent_id = "root";
        c.hotspot_gained = (bool)hs_status_r[i];
        c.next_event_time = event_times_r[i];
        c.alive = true; c.generation = 0;
        c.alleles.resize(p.n_alleles);
        for (int ai = 0; ai < p.n_alleles; ++ai)
            c.alleles[ai] = rlist_to_sequence(seqs_for_cell[p.allele_names[ai]]);

        int pidx = (int)pool.size();
        pool.push_back(std::move(c));
        pq.push({event_times_r[i], pidx, 0});
        ++alive_count;
    }

    // Load existing history (from initialization)
    int h_n0 = as<int>(sim_state_r["h_n"]);
    if (h_n0 > 0) {
        CharacterVector hcid = sim_state_r["h_cell_id"];
        CharacterVector hpid = sim_state_r["h_parent_id"];
        LogicalVector   hbfb = sim_state_r["h_bfb_event"];
        LogicalVector   hwgd = sim_state_r["h_wgd_event"];
        CharacterVector hcn  = sim_state_r["h_cn_event"];
        CharacterVector hchr = sim_state_r["h_chr_allele"];
        hist.reserve(h_n0 + as<int>(ip["max_cells"]) * 4);
        for (int i = 0; i < h_n0; ++i) {
            HistRec r;
            r.cell_id   = as<std::string>(hcid[i]);
            r.parent_id = as<std::string>(hpid[i]);
            r.bfb_event = (bool)hbfb[i];
            r.wgd_event = (bool)hwgd[i];
            r.cn_event  = as<std::string>(hcn[i]);
            r.chr_allele = (hchr[i] == NA_STRING) ? "" : as<std::string>(hchr[i]);
            hist.push_back(r);
        }
    }

    int next_id   = as<int>(sim_state_r["next_cell_id"]);
    double t_cur  = as<double>(sim_state_r["time"]);

    // ── Gillespie loop ────────────────────────────────────────────────────────
    while (t_cur < p.max_time && alive_count > 0 && alive_count < p.max_cells) {
        // Pop next valid event (lazy deletion)
        bool found = false;
        while (!pq.empty()) {
            PQEntry4 entry = pq.top(); pq.pop();
            if (entry.pool_idx >= (int)pool.size()) continue;
            CellD& cell = pool[entry.pool_idx];
            if (!cell.alive || cell.generation != entry.generation) continue;

            t_cur = entry.time;
            if (t_cur >= p.max_time) { found = false; break; }
            found = true;

            // Birth vs death decision
            double br = p.birth_rate * (1.0 + p.pos_sel_rate * cell.hotspot_gained);
            double dr = p.death_rate * (1.0 + p.neg_sel_rate * cell.hotspot_gained);
            bool is_birth = (R::unif_rand() < br / (br + dr));

            if (is_birth) {
                process_birth4(pool, pq, hist, alive_count,
                               next_id, wgd_available, p,
                               entry.pool_idx, t_cur);
            } else {
                cell.alive = false;
                cell.generation++;
                cell.alleles.clear();
                --alive_count;
            }
            break;
        }
        if (!found) break;
    }

    // ── Pack results ─────────────────────────────────────────────────────────
    // Collect alive cells
    std::vector<int> alive_idx;
    for (int i = 0; i < (int)pool.size(); ++i)
        if (pool[i].alive) alive_idx.push_back(i);

    int n_alive = (int)alive_idx.size();

    CharacterVector out_cell_ids(n_alive);
    List out_cell_seqs(n_alive);
    for (int i = 0; i < n_alive; ++i) {
        const CellD& c = pool[alive_idx[i]];
        out_cell_ids[i] = c.id;

        List seqs(p.n_alleles);
        CharacterVector anames(p.n_alleles);
        for (int ai = 0; ai < p.n_alleles; ++ai) {
            seqs[ai] = sequence_to_rlist(c.alleles[ai]);
            anames[ai] = p.allele_names[ai];
        }
        seqs.attr("names") = anames;
        out_cell_seqs[i] = seqs;
    }
    out_cell_seqs.attr("names") = out_cell_ids;

    // Build history vectors
    int h_n = (int)hist.size();
    CharacterVector h_cid(h_n), h_pid(h_n), h_cn(h_n), h_chr(h_n);
    LogicalVector   h_bfb(h_n), h_wgd(h_n);
    for (int i = 0; i < h_n; ++i) {
        h_cid[i] = hist[i].cell_id;
        h_pid[i] = hist[i].parent_id;
        h_bfb[i] = hist[i].bfb_event;
        h_wgd[i] = hist[i].wgd_event;
        h_cn[i]  = hist[i].cn_event;
        h_chr[i] = hist[i].chr_allele.empty() ? NA_STRING
                                               : String(hist[i].chr_allele);
    }

    // Update input_parameters with final wgd_available
    List ip_out = clone(ip);
    ip_out["wgd_available"] = wgd_available;

    return List::create(
        Named("time")                  = t_cur,
        Named("cell_ids")              = out_cell_ids,
        Named("cell_sequences")        = out_cell_seqs,
        Named("cell_next_event_times") = NumericVector(0),
        Named("hotspot_status")        = LogicalVector(0),
        Named("h_cell_id")             = h_cid,
        Named("h_parent_id")           = h_pid,
        Named("h_bfb_event")           = h_bfb,
        Named("h_wgd_event")           = h_wgd,
        Named("h_cn_event")            = h_cn,
        Named("h_chr_allele")          = h_chr,
        Named("h_n")                   = h_n,
        Named("next_cell_id")          = next_id,
        Named("input_parameters")      = ip_out
    );
}
