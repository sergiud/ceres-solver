#if defined(_MSC_VER)
#   pragma warning(pop)
#elif defined(__clang__)
#   pragma clang diagnostic pop
#elif defined(__GNUC__)
#   if (__GNUC__ >= 4) && (__GNUC_MINOR__ >= 6)
#       pragma GCC diagnostic pop
#   endif // (__GNUC__ >= 4) && (__GNUC_MINOR__ >= 6)
#endif // defined(_MSC_VER)
