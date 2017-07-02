#if defined(_MSC_VER)
#   pragma warning(push)
#   pragma warning(disable: 4275 4251 4231)
#elif defined(__clang__)
#   pragma clang diagnostic push
#elif defined(__GNUC__)
#   if (__GNUC__ >= 4) && (__GNUC_MINOR__ >= 6)
#       pragma GCC diagnostic push
#       if (__GNUC__ >= 4) && (__GNUC_MINOR__ >= 8)
#           pragma GCC diagnostic ignored "-Wpedantic"
#       endif // (__GNUC__ >= 4) && (__GNUC_MINOR__ >= 8)
#   endif // (__GNUC__ >= 4) && (__GNUC_MINOR__ >= 6)
#endif // defined(_MSC_VER)
