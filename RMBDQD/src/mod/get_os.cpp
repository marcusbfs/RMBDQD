extern"C" void get_os_(int &i) {

#if defined(__linux) || defined(__unix__) || defined(__unix) || defined(unix)
    i = 1;
#elif defined(_WIN32)
    i = 2;
#else
    i = 3;
#endif
}
