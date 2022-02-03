#ifdef MPI 
#define io_begin IF (NODE_ME==IONODE) THEN
#define io_end   ENDIF
#define do_io    IF (NODE_ME==IONODE)
#else
#define io_begin
#define io_end
#define do_io
#endif
