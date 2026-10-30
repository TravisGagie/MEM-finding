// Builds the RLBWT of a collection and its shortcuts for several window
// lengths ell and pruning thresholds theta.
//
//   build_index --text FILE --reference-length N --out DIR
//               [--bwt FILE | --rlbwt FILE] [--ell 25,50] [--theta 0,400] [--max-tier 3]
//               [--reversed] [--tag STRING]
//
// --text            the collection: one sequence per line, or FASTA (text.hpp)
// --reference-length the reference R is the first N characters of T
// --bwt             a BWT of --text with one byte per row, e.g. from Big-BWT;
//                   --text must then have one sequence per line, with no
//                   blank lines and a final newline. Without --bwt, the BWT
//                   is computed with libsais.
// --rlbwt           an RLBWT written by an earlier run, to reuse
// --tag             appended to the names of the shortcut files, e.g. _t6
// --reversed        index the reversed collection (text.hpp) instead of T,
//                   with the reverse of R as reference, for searching
//                   forward (mems.hpp); its anchors are at the same
//                   boundaries of R as those of the index of T. A --bwt file must then be a BWT of
//                   `tac FILE | rev`.
//
// Writes DIR/rlbwt.bin, DIR/ell<L>_theta<T><TAG>.shortcuts for each pair, and
// one line per pair to DIR/build.csv.
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "build.hpp"
#include "rlbwt.hpp"
#include "text.hpp"

using namespace rs;
namespace fs = std::filesystem;

namespace {

double seconds_since(std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

std::vector<int64_t> parse_list(const std::string &s) {
  std::vector<int64_t> v;
  std::stringstream ss(s);
  for (std::string item; std::getline(ss, item, ',');) v.push_back(std::stoll(item));
  return v;
}

int usage() {
  std::cerr << "usage: build_index --text FILE --reference-length N --out DIR\n"
               "                   [--bwt FILE | --rlbwt FILE] [--ell 25,50] [--theta 0,400] [--max-tier 3]\n"
               "                   [--reversed] [--tag STRING]\n";
  return 1;
}

}  // namespace

int main(int argc, char **argv) {
  std::string text_path, bwt_path, rlbwt_path, out_dir, tag;
  int64_t reference_length = -1, max_tier = 3;
  std::vector<int64_t> ells = {25}, thetas = {0};
  bool reversed = false;
  std::vector<std::string> args(argv + 1, argv + argc);
  if (std::find(args.begin(), args.end(), "--reversed") != args.end()) {
    reversed = true;
    args.erase(std::remove(args.begin(), args.end(), "--reversed"), args.end());
  }
  for (size_t i = 0; i + 1 < args.size(); i += 2) {
    const std::string &flag = args[i], &value = args[i + 1];
    if (flag == "--text") text_path = value;
    else if (flag == "--bwt") bwt_path = value;
    else if (flag == "--rlbwt") rlbwt_path = value;
    else if (flag == "--out") out_dir = value;
    else if (flag == "--reference-length") reference_length = std::stoll(value);
    else if (flag == "--max-tier") max_tier = std::stoll(value);
    else if (flag == "--ell") ells = parse_list(value);
    else if (flag == "--theta") thetas = parse_list(value);
    else if (flag == "--tag") tag = value;
    else return usage();
  }
  if (args.size() % 2 == 1 || text_path.empty() || out_dir.empty() || reference_length < 0) return usage();

  try {
    fs::create_directories(out_dir);
    auto t0 = std::chrono::steady_clock::now();
    RLBWT rlbwt;
    if (!rlbwt_path.empty()) {
      rlbwt = RLBWT::load(rlbwt_path);
    } else if (!bwt_path.empty()) {
      rlbwt = RLBWT::from_bwt_file(bwt_path);
      // Catches a BWT of a different file, or of a file with blank lines or
      // no final newline, whose text differs from the T that --text gives.
      if (rlbwt.size() != text_length(text_path))
        throw std::runtime_error("the BWT has " + std::to_string(rlbwt.size()) + " rows but T has " +
                                 std::to_string(text_length(text_path)) + " characters");
    } else {
      std::vector<uint8_t> T = read_text(text_path);
      rlbwt = RLBWT::from_text(reversed ? reverse_collection(T) : T);
    }
    const std::string rlbwt_out = (fs::path(out_dir) / "rlbwt.bin").string();
    if (rlbwt_path.empty() || fs::absolute(rlbwt_path) != fs::absolute(rlbwt_out)) rlbwt.save(rlbwt_out);
    const int64_t rlbwt_bytes = static_cast<int64_t>(fs::file_size(rlbwt_out));
    std::cerr << "RLBWT: n = " << rlbwt.size() << ", r = " << rlbwt.num_runs() << ", " << rlbwt_bytes
              << " bytes [" << seconds_since(t0) << " s]\n";

    std::vector<uint8_t> R = read_reference(text_path, reference_length);
    if (reversed) std::reverse(R.begin(), R.end());

    const std::string csv_path = (fs::path(out_dir) / "build.csv").string();
    const bool new_csv = !fs::exists(csv_path);
    std::ofstream csv(csv_path, std::ios::app);
    if (new_csv) {
      csv << "ell,theta,max_tier,n,runs,rlbwt_bytes,shortcuts_bytes,total_bytes,entrances,entrance_rows,"
             "compute_seconds,prune_seconds,anchors,dropped_unusable,discarded_colliding,dropped_by_area,"
             "backward_steps,deaths\n";
    }

    for (int64_t ell : ells) {
      auto t1 = std::chrono::steady_clock::now();
      ConstructionStats stats;
      UnprunedShortcuts unpruned = compute_shortcuts(rlbwt, R, ell, max_tier, stats, reversed);
      const double compute_seconds = seconds_since(t1);
      std::cerr << "ell = " << ell << ": " << unpruned.entrances.size() << " entrances, " << stats.deaths
                << " deaths [" << compute_seconds << " s]\n";

      for (int64_t theta : thetas) {
        auto t2 = std::chrono::steady_clock::now();
        int64_t dropped = 0;
        Shortcuts sc = prune(unpruned, theta, dropped);
        const double prune_seconds = seconds_since(t2);
        const std::string path =
            (fs::path(out_dir) / ("ell" + std::to_string(ell) + "_theta" + std::to_string(theta) + tag + ".shortcuts")).string();
        sc.save(path);
        const int64_t bytes = static_cast<int64_t>(fs::file_size(path));
        csv << ell << ',' << theta << ',' << max_tier << ',' << rlbwt.size() << ',' << rlbwt.num_runs() << ','
            << rlbwt_bytes << ',' << bytes << ',' << rlbwt_bytes + bytes << ',' << sc.num_entrances() << ','
            << sc.num_rows() << ',' << compute_seconds << ',' << prune_seconds << ',' << stats.anchors << ','
            << stats.dropped_unusable << ',' << stats.discarded_colliding << ',' << dropped << ','
            << stats.backward_steps << ',' << stats.deaths << '\n';
        std::cerr << "  theta = " << theta << ": " << sc.num_entrances() << " entrances, " << bytes << " bytes\n";
      }
    }
  } catch (const std::exception &e) {
    std::cerr << "error: " << e.what() << '\n';
    return 1;
  }
  return 0;
}
