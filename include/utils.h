/* 
* This way we can make things still static while making them accessible for unit testing.
*/
#if UNITTEST
#define ustatic 
#else
#define ustatic static
#endif