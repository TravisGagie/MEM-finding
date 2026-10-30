// Binary reading and writing of integers and vectors, for save/load.
#pragma once
#include <cstdint>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <vector>

namespace rs::io {

template <typename T>
void write(std::ostream &out, T value) {
  out.write(reinterpret_cast<const char *>(&value), sizeof(T));
}

template <typename T>
T read(std::istream &in) {
  T value{};
  in.read(reinterpret_cast<char *>(&value), sizeof(T));
  if (!in) throw std::runtime_error("unexpected end of file");
  return value;
}

template <typename T>
void write_vector(std::ostream &out, const std::vector<T> &v) {
  write(out, static_cast<int64_t>(v.size()));
  if (!v.empty()) out.write(reinterpret_cast<const char *>(v.data()), static_cast<std::streamsize>(v.size() * sizeof(T)));
}

template <typename T>
std::vector<T> read_vector(std::istream &in) {
  std::vector<T> v(static_cast<size_t>(read<int64_t>(in)));
  if (!v.empty()) in.read(reinterpret_cast<char *>(v.data()), static_cast<std::streamsize>(v.size() * sizeof(T)));
  if (!in) throw std::runtime_error("unexpected end of file");
  return v;
}

}  // namespace rs::io
