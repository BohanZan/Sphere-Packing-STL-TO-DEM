#pragma once
#include <exception>
#include <string>

namespace sp {
// Intel LLVM 2025.2 on Windows corrupts ASan-instrumented catch parameters
// (upstream LLVM PR 159618). Isolate just exception-message extraction; the
// packing, geometry, I/O, and callers retain full ASan instrumentation.
// Call only while handling an active exception. Do not inline this boundary.
#if defined(_WIN32) && defined(__INTEL_LLVM_COMPILER)
__declspec(noinline) __attribute__((no_sanitize("address")))
#endif
inline std::string currentExceptionMessage() {
  try {
    throw;
  } catch (const std::exception &error) {
    return error.what();
  } catch (...) {
    return "Unknown exception";
  }
}
} // namespace sp
