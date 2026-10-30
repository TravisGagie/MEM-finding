// The indexed text and the reference.
//
// A collection file holds one sequence per line, or is in FASTA format
// (if its first byte is '>'). Characters are taken as they are, so a and
// A differ; scripts/prepare_collection.py upper-cases.
// The indexed text is
//     T = S_1 # S_2 # ... # S_k # $
// where # is kSeparator, written after every sequence, and $ is
// kTerminator, which occurs once and is the smallest character. Neither
// can occur in a pattern, so no match crosses from one sequence to the
// next.
//
// The reference R is a prefix of T (paper, Section 4): the first
// `reference_length` characters. If R spans several sequences, the
// separators between them are not A/C/G/T, so no window or tier string
// crosses from one to the next.
#pragma once
#include <cstdint>
#include <string>
#include <vector>

namespace rs {

constexpr uint8_t kTerminator = 0x00;
constexpr uint8_t kSeparator = 0x01;

inline bool is_acgt(uint8_t c) { return c == 'A' || c == 'C' || c == 'G' || c == 'T'; }

// Reads a collection file and returns T.
std::vector<uint8_t> read_text(const std::string &path);

// The length of T, without building it.
int64_t text_length(const std::string &path);

// Returns the first `length` characters of T, reading only as much of
// the collection file as needed.
std::vector<uint8_t> read_reference(const std::string &path, int64_t length);

// The reversed collection, for searching forward (mems.hpp): for
// T = S_1 # ... # S_k # $, returns
//     rev(S_k) # ... # rev(S_1) # $,
// the text of the collection file with its lines in reverse order and
// each line reversed (as `tac FILE | rev` writes it). A string occurs in
// T if and only if its reverse occurs in the reversed collection, and the
// reverse of the reference R is a substring of it.
std::vector<uint8_t> reverse_collection(const std::vector<uint8_t> &T);

}  // namespace rs
