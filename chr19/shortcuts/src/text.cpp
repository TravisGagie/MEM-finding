#include "text.hpp"

#include <algorithm>
#include <fstream>
#include <stdexcept>

namespace rs {

namespace {

// Calls emit(c) for each character of T in order, without the final $,
// until emit returns false.
template <typename Emit>
void for_each_text_char(const std::string &path, Emit emit) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("cannot open " + path);
  int first = in.peek();
  bool fasta = first == '>';
  std::string line;
  bool in_sequence = false;  // FASTA: some sequence line seen since the last header
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (fasta && !line.empty() && line[0] == '>') {
      if (in_sequence && !emit(kSeparator)) return;
      in_sequence = false;
      continue;
    }
    if (line.empty()) continue;
    for (char c : line) {
      uint8_t b = static_cast<uint8_t>(c);
      if (b == kSeparator || b == kTerminator) throw std::runtime_error(path + " contains a byte 0x00 or 0x01");
      if (!emit(b)) return;
    }
    if (fasta) {
      in_sequence = true;
    } else if (!emit(kSeparator)) {
      return;
    }
  }
  if (fasta && in_sequence) emit(kSeparator);
}

}  // namespace

std::vector<uint8_t> read_text(const std::string &path) {
  std::vector<uint8_t> text;
  for_each_text_char(path, [&](uint8_t c) {
    text.push_back(c);
    return true;
  });
  text.push_back(kTerminator);
  return text;
}

int64_t text_length(const std::string &path) {
  int64_t n = 1;  // the terminator
  for_each_text_char(path, [&](uint8_t) {
    ++n;
    return true;
  });
  return n;
}

std::vector<uint8_t> read_reference(const std::string &path, int64_t length) {
  std::vector<uint8_t> ref;
  ref.reserve(static_cast<size_t>(length));
  if (length > 0) {
    for_each_text_char(path, [&](uint8_t c) {
      ref.push_back(c);
      return static_cast<int64_t>(ref.size()) < length;
    });
  }
  if (static_cast<int64_t>(ref.size()) < length) throw std::runtime_error("the reference is longer than the text");
  return ref;
}

std::vector<uint8_t> reverse_collection(const std::vector<uint8_t> &T) {
  if (T.empty() || T.back() != kTerminator) throw std::invalid_argument("T must end with the terminator");
  // T without $ is S_1 # ... # S_k #; reversed, it is # rev(S_k) ... # rev(S_1).
  std::vector<uint8_t> rev(T.rbegin() + 1, T.rend());
  if (!rev.empty()) std::rotate(rev.begin(), rev.begin() + 1, rev.end());  // move the leading # to the end
  rev.push_back(kTerminator);
  return rev;
}

}  // namespace rs
