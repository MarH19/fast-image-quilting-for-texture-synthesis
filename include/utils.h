#define REL_TOL 1e-9
#define ABS_TOL 0.0
#define ABS(a) (((a) < 0) ? -1 * (a) : (a))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define IS_CLOSE(a, b) (ABS(a - b) <= MAX(REL_TOL * MAX(ABS(a), ABS(b)), ABS_TOL))

/* This way we can make things still static while making them accessible for unit testing. */
#if UNITTEST
#define ustatic 
#else
#define ustatic static
#endif