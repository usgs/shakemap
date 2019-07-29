typedef struct cpoint {
        double   x;        /* x coordinate the contour intersects the edge */
        double   y;        /* y coordinate the contour intersects the edge */
        int     ti;        /* triangle index -- temporary use only */
        int     ei;        /* local edge index of the intersected edge */
        int     gei;       /* global edge index of the intersected edge */
        int     order;     /* indicating the segment's handedness */
        struct cpoint *next;        /* next point in the segment */
} CPoint;

typedef struct cseg {
        double xmin;
        double xmax;
        double ymin;
        double ymax;
        int    ci;
        int    nesting;
        int    neg_depth;
        int    pos_depth;
        struct cseg *neg_holes;
        struct cseg *neg_next;
        struct cseg *pos_holes;
        struct cseg *pos_next;
        struct cseg *tmp_hole;
        struct cseg *next;
        CPoint *first;
        CPoint *last;
} CSeg;

typedef struct cinfo {
        double cv;        /* value of the contour level */
        CSeg *open;        /* chain of open contours at this level */
        CSeg *closed;        /* chain of closed contours at this level */
} CInfo;

typedef struct cresult {
    size_t ncarray;
    CInfo *contour0;
    CInfo *carray;
} CResult;

CResult contour_grid(double *grid, size_t gnx, size_t gny, double gdx,          
                     double gdy, double ul_x, double ul_y,                      
                     double *contour_levels, size_t ncont,                      
                     int outstyle, int verb);
