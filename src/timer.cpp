
#ifdef _WIN32
#include <windows.h>
double timer() {
	return (double)GetTickCount()/1000.0;
}
#else
#include <sys/time.h>
#include <sys/resource.h>
double timer() {
  struct rusage r;
  getrusage(0, &r);
  return (double)(r.ru_utime.tv_sec+r.ru_utime.tv_usec/1000000.0);
}
#endif
