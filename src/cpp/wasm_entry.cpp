#include <cstddef>

int run_analysis(int argc, const char* const* argv);

extern "C" int reax_run_analysis(int argc, const char** argv) {
    return run_analysis(argc, argv);
}
