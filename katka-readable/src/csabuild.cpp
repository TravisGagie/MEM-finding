// katka-csabuild: build the RLCSA (rlcsa.hpp) from a BWT file, for katka-classify -C.
// usage: katka-csabuild bwt out.csa [explicit|ef]     (default: explicit run lists, as in the paper)
#include "rlcsa.hpp"

int main(int argc, char **argv) {
    if (argc < 3) katka::die("usage: katka-csabuild bwt out.csa [explicit|ef]");
    bool ef = argc > 3 && std::string(argv[3]) == "ef";
    double t0 = katka::now();
    katka::rlcsa X;
    X.build(argv[1], ef);
    X.save(argv[2]);
    printf("n=%lu r=%lu %s RLCSA bytes=%lu time=%.1fs\n", X.n, X.R, ef ? "Elias-Fano" : "explicit", X.bytes(), katka::now() - t0);
    return 0;
}
