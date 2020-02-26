#if !defined(__has_include)
#define __has_include(x) 0
#endif
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <strings.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>

#include "contour.h"

#ifndef DBL_MAX
#define DBL_MAX         1.7976931348623157E+308 /* max decimal value of a */
                                                /* "double" */
#endif /* DBL_MAX */

#ifndef DBL_MIN
#define DBL_MIN         -DBL_MAX
#endif /* DBL_MIN */

#define NEG_HOLES        0
#define POS_HOLES        1

#define NGRAB 1000

/*
 * In this contouring scheme, the grid is divided into a series of
 * boxes, one grid point on each corner (p0 - p3).  A fifth point (p4)
 * is placed at the center of the box, and is the average of the
 * the other four points (this is an approximation, and could
 * probably be improved by fitting a surface).
 *
 * The box is then divided into four triangular regions (t0 - t3)
 * with eight total edges (e0 - e7), as illustrated below:

       p0        e1        p1
        +------------------+
        |\                /|
        | \      t1      / |
        |  \            /  |
        |   \e4      e5/   |
        |    \        /    |
        |     \      /     |
        |      \    /      |
        |       \  /       |
     e0 |  t0    \/    t2  | e2
        |        /\p4      |
        |       /  \       |
        |      /    \      |
        |     /      \     |
        |    /        \    |
        |   /e7      e6\   |
        |  /            \  |
        | /     t3       \ |
        |/                \|
        +------------------+
       p3        e3         p2

 * Contouring takes place by determining if a contour of a particular
 * interval passes through one of the four outside edges of the box
 * (e0 - e3), and then is traced through the box until it reaches
 * another outside edge.  For example, if a contour passes through
 * e0, it must then intersect one of e4 or e7.  If it intersects e4,
 * then it must intersect e1 or e5, and so on.  We use linear 
 * interpolation, so no edge will be penetrated by a contour of a
 * given level more than once.
 *
 * The four outside edges can each be given a unique global identifier
 * (see get_global_edge(), below).  By contouring each box, then
 * connecting the the ends of the contours by means of the global
 * identifiers, complete contours may be traced through the grid.
 */

/*
 * These are the edges that make up each triangle (by index)
 */
int tri_egs[4][3] = {
        { 0, 4, 7 },
        { 1, 4, 5 },
        { 2, 5, 6 },
        { 3, 6, 7 } };

/*
 * These are the points that make up each edge (by index)
 * Since we would ultimately like to create "right-handed" polygons
 * (i.e. when moving from point to point the interior of the polygon
 * is to the right), we want the points for the four outer edges to
 * be in a known order, so we choose them so that when standing on
 * the edge and facing inward, the first coordinate is on the right,
 * and the second is on the left.  The point order of the interior
 * edges doesn't matter.
 */
struct {
        int p1;
        int p2;
} eg_pts[] = {
        { 3, 0 },
        { 0, 1 },
        { 1, 2 },
        { 2, 3 },
        { 0, 4 },
        { 1, 4 },
        { 2, 4 },
        { 3, 4 } };

/*
 * Each interior edge (4-7) belongs to two triangles.  Given
 * one of those triangles, we want to know the other (so that
 * we may trace the contour from triangle to triangle).  This
 * array is set up so that for a given edge index, and a given
 * triangle index, the dereferenced value is the other triangle's 
 * index.  Got that?
 */
int eg_tri[8][4] = {
        { -1, -1, -1, -1 },
        { -1, -1, -1, -1 },
        { -1, -1, -1, -1 },
        { -1, -1, -1, -1 },
        {  1,  0, -1, -1 },
        { -1,  2,  1, -1 },
        { -1, -1,  3,  2 },
        {  3, -1, -1,  0 } };

CInfo *carray = NULL;
int ncarray = 0;
int ncpoint = NGRAB;
CPoint *cpbuf;
CInfo contour0;
int ncontours = 0;
double avecontint = 0;
double contour_delta;
double xorig, yorig;
double dx, dy;
int nx, ny;
int verbose;

typedef struct {
        double        x;
        double        y;
        double        z;
} XYZPt;

struct {
        double xmin;
        double xmax;
        double ymin;
        double ymax;
} gedge;

void contour_box(int ix, int iy, double *grid);
void get_p2(XYZPt pt[], double cv, CPoint *p1, CPoint **p2);
void fill_cpoint(XYZPt pt[], int ei, int ti, double cv, int order, CPoint **p);
void store_segment(CSeg *newseg, int ci);
int get_global_edge(int ix, int iy, int ei);
CPoint *get_new_cpoint(void);
CPoint *dup_cpoint(CPoint *p);
CSeg *get_new_cseg(void);
void reverse_cseg(CSeg *seg);
CSeg *dup_cseg(CSeg *s);
CSeg *get_outer_boundry_polygon(void);
void fixval(int i, int j, double *g);
double epps(void);
void set_depressions(void);
int is_contained(CSeg *out, CSeg *in, CPoint *cp);
int is_depression(CSeg *out);
void set_xyminmax(CSeg *s);
void add_contour(double cval);
void finish_contours(void);
void finish_polygons(void);
void connect_segments(CInfo *ca);
int grid_edge(double x, double y);
void clean_up_polygons(void);
void nest_polygons(void);
int almost_equal(double a, double b, int maxUlps);
int64_t abs64(int64_t a);

CResult contour_grid(double *grid, size_t gnx, size_t gny, double gdx, 
                     double gdy, double ul_x, double ul_y, 
                     double *contour_levels, size_t ncont,
                     int outstyle, int verb) {

  int i, j, k; 
  size_t ui;
  double cv;
  CResult cresult;

  carray = NULL;
  ncarray = 0;
  contour0.cv = 0;
  contour0.open = NULL;
  contour0.closed = NULL;
  ncontours = 0;
  ncpoint = NGRAB;
  cpbuf = NULL;

  /* 
   * Using a fixed seed means we get the same random numbers
   * every time, but that gives us the same results every time
   * (which is good) and I don't think the systematic errors it 
   * introduces are a problem.
   */
  srand48(444);
  /*
  srand48(getpid());
  */

  /*
   * Set some useful variables...
   */
  dx = gdx;
  dy = -1 * gdy;
  nx = gnx;
  ny = gny;
  xorig = ul_x;
  yorig = ul_y;

  /*
   * These values are at the borders of the grid
   * Note: the strange calculation (i.e. (nx-2)*dx then later
   * adding dx again) is to duplicate the way this value will
   * be calculated later on, and since we compare this value to
   * another floating point value for equality, we don't want to 
   * introduce even the smallest round-off difference.
   */
  if( outstyle == 2 || outstyle == 3 || outstyle == 4 ) {
    gedge.xmin = xorig;
    gedge.xmax = (nx - 2) * dx + xorig;
    gedge.xmax += dx;
    gedge.ymax = yorig;
    gedge.ymin = (ny - 2) * dy + yorig;
    gedge.ymin += dy;
  }

  /*
   * look for NaNs, set them to something sane
   */
  for(j = 0; j < ny; j++) {
    for(i = 0; i < nx; i++) {
      if (isnan(grid[i + j*nx]) )
        fixval(i, j, grid);
    }
  }

  /*
   * get contour info
   */
  for(ui = 0; ui < ncont; ui++) {
      add_contour(contour_levels[ui]);
  }

  /*
   * Compute the average contour interval
   */
  if( ncontours == 1 ) {
    avecontint = carray[0].cv;
  } else {
    for( i = 1; i < ncontours; i++ ) {
      avecontint += fabs(carray[i].cv - carray[i - 1].cv);
    }
    avecontint /= (ncontours - 1);
  }
  if( avecontint == 0 )
    avecontint = 1.0;
  contour_delta = avecontint * 1.e-4;

  /*
   * We can't handle grid points that fall exactly on contours,
   * so we nudge them off, just a bit...
   */
  for(k = 0; k < ncontours; k++) {
    cv = carray[k].cv;
    for(j = 0; j < ny; j++) {
      for(i = 0; i < nx; i++) {
        if( grid[i + j*nx] > cv - contour_delta 
         && grid[i + j*nx] < cv + contour_delta ) {
          grid[i + j*nx] += epps();
          /*
           * We want to recheck this point to make sure we changed
           * the value within the resolution of the machine's
           * floating point representation.
           */
          i--;
        }
      }
    }
  }
  if( verbose ) fprintf(stderr, "Moved nodes off contour values.\n");

  /*
   * Do the actual contouring
   */
  for(j = 0; j < ny - 1; j++) {
    for(i = 0; i < nx - 1; i++) {
      contour_box(i, j, grid);
    }
  }
  if( verbose ) fprintf(stderr, "Did contouring.\n");

  /*
   * output the results
   */
  cresult.ncarray = ncontours;
  if( outstyle == 1 ) {
    cresult.carray = carray;
    cresult.contour0 = &contour0;
    return cresult;
  } else if( outstyle == 2 ) {
    finish_contours();
    cresult.carray = carray;
    cresult.contour0 = &contour0;
    return cresult;
  } else {
    finish_polygons();
    cresult.carray = carray;
    cresult.contour0 = &contour0;
    return cresult;
  }
}

/*
 * ci => contour index
 * cv => contour value
 * ti => triangle index
 * ei => edge index
 */

void contour_box(int ix, int iy, double *grid) {

  int lp1, lp2, ci, ei, order, tmp;
  int done[4];
  double cv;
  CSeg *seg;
  XYZPt pt[5];
  CPoint *p1, *p2, *last;
  double x = ix * dx + xorig;
  double y = iy * dy + yorig;

  pt[0].x = x;
  pt[1].x = x + dx;
  pt[2].x = x + dx;
  pt[3].x = x;
  pt[4].x = x + dx / 2.0;
  pt[0].y = y;
  pt[1].y = y;
  pt[2].y = y + dy;
  pt[3].y = y + dy;
  pt[4].y = y + dy / 2.0;
  pt[0].z = grid[iy * nx + ix];
  pt[1].z = grid[iy * nx + ix + 1];
  pt[2].z = grid[(iy + 1) * nx + ix + 1];
  pt[3].z = grid[(iy + 1) * nx + ix];
  pt[4].z = (pt[0].z + pt[1].z + pt[2].z + pt[3].z) / 4.0;

  /*
   * Loop over the contour values
   */
  for(ci = 0; ci < ncontours; ci++) {
    cv = carray[ci].cv;

    if( pt[4].z > cv - contour_delta && pt[4].z < cv + contour_delta ) {
      pt[4].z += epps();
      ci--;
      continue;
    }

    done[0] = done[1] = done[2] = done[3] = 0;
    /*
     * Loop over the four outside edges of this block
     */
    for(ei = 0; ei < 4; ei++) {
      if( done[ei] )
        continue;
      lp1 = eg_pts[ei].p1;
      lp2 = eg_pts[ei].p2;
      if( pt[lp1].z >= pt[lp2].z ) {
        /* normal order */
        order = 1;
      } else {
        /* reverse order */
        tmp = lp1;
        lp1 = lp2;
        lp2 = tmp;
        order = -1;
      }
      if( cv > pt[lp1].z || cv < pt[lp2].z )
        continue;
      /*
       * Fill the first point with initial info; note that this
       * is an outside edge, so ei == ti
       */
      fill_cpoint(pt, ei, ei, cv, order, &p1);
      last = p1;
      /*
       * global edges only matter on the anchor points of a segment
       */
      p1->gei = get_global_edge(ix, iy, ei);

      /* 
       *while p2 not on edge 0 -> 3, continue through the current box 
       */
      do {

        get_p2(pt, cv, last, &p2);

        /* add p2 to current segment */
        last->next = p2;
        last = p2;
      } while (last->ei > 3);
      last->gei = get_global_edge(ix, iy, last->ei);
      /* 
       * segment is now complete, store it 
       */
      done[last->ei] = 1;
      seg = get_new_cseg();
      seg->first = p1;
      seg->last  = last;
      seg->next = carray[ci].open;
      carray[ci].open = seg;
      store_segment(seg, ci);
    }
  }
  return;
}

void get_p2(XYZPt pt[], double cv, CPoint *p1, CPoint **p2) {

  int i, ei, lp1, lp2, tmp, order;

  for( i = 0; i < 3; i++ ) {
    if((ei = tri_egs[p1->ti][i]) == p1->ei)
      continue;
    lp1 = eg_pts[ei].p1;
    lp2 = eg_pts[ei].p2;
    if( pt[lp1].z >= pt[lp2].z ) {
      /* normal order */
      order = 1;
    } else {
      /* reverse order */
      tmp = lp1;
      lp1 = lp2;
      lp2 = tmp;
      order = -1;
    }
    if( cv > pt[lp1].z || cv < pt[lp2].z )
      continue;
    fill_cpoint(pt, ei, p1->ti, cv, order, p2);
    break;
  }
  if( i >= 3 ) {
    fprintf( stderr, "get_p2: couldn't find another crossing point\n");
    exit(1);
  }
  /*
   * He we set up the next iteration: if this edge is an interior
   * edge, we will need to do the next computation on the adjacent
   * triangle, rather than the current one.
   */
  if( (*p2)->ei > 3 )
    (*p2)->ti = eg_tri[(*p2)->ei][(*p2)->ti];
  return;
}

void fill_cpoint(XYZPt pt[], int ei, int ti, double cv, int order, CPoint **p) {

  int lp1, lp2;
  double f;

  *p = get_new_cpoint();
  (*p)->order = order;
  (*p)->ti    = ti;
  (*p)->ei    = ei;
  if(order > 0) {
    lp1 = eg_pts[ei].p1;
    lp2 = eg_pts[ei].p2;
  } else {
    lp1 = eg_pts[ei].p2;
    lp2 = eg_pts[ei].p1;
  }
  f = (pt[lp1].z - cv) / (pt[lp1].z - pt[lp2].z);
  if(ei == 0 || ei == 2) {
    /* vertical edge */
    (*p)->x = pt[lp1].x;
    (*p)->y = pt[lp2].y * f + pt[lp1].y * (1.0 - f);
  } else if(ei == 1 || ei == 3) {
    /* horizontal edge */
    (*p)->x = pt[lp2].x * f + pt[lp1].x * (1.0 - f);
    (*p)->y = pt[lp1].y;
  } else {
    /* diagonal edge */
    (*p)->x = pt[lp2].x * f + pt[lp1].x * (1.0 - f);
    (*p)->y = pt[lp2].y * f + pt[lp1].y * (1.0 - f);
  }
  return;
}

/*
 * store_segment tries to add a segment to an existing segment, if
 * no current segment matches, it starts a new segment.  If, by
 * adding the new segment, a segment becomes closed, the segment
 * is moved to the closed segment list.
 */
void store_segment(CSeg *newseg, int ci) {

  CSeg *seg, *prev;
  int didmerge = 1;

  if( newseg->first->order < 0 )
    reverse_cseg(newseg);

  /*
   * Try to find a matching segment on the open segment list
   */
  while( didmerge ) {
    for( prev = 0, seg = carray[ci].open; seg; prev = seg, seg = seg->next ) {
      if( seg == newseg )
        continue;
      if( newseg->first->gei == seg->last->gei ) {
        /* 
         * add the existing segment to the front of the new one
         * skip the first point in the new segment since it is a duplicate
         * WARNING: memory leak -- lose track of 1 CPoint
         */
        seg->last->next = newseg->first->next;
        newseg->first = seg->first;
        break;
      } else if( newseg->last->gei == seg->first->gei ) {
        /* 
         * add new segment to the front of the existing one 
         * skip the first point in the old segment since it is a duplicate
         * WARNING: memory leak -- lose track of 1 CPoint
         */
        newseg->last->next = seg->first->next;
        newseg->last = seg->last;
        break;
      }
    }
    if( seg ) {
      /*
       * We did a merge, so first we want to get rid of the old
       * segment, which is 'seg'
       */
      if( prev ) 
        prev->next = seg->next;
      else
        carray[ci].open = seg->next;
      free(seg);
      /*
       * Now see if the merged segment is closed
       */
      if( newseg->first->gei == newseg->last->gei ) {
        /* 
         * pop the segment out of the list and add it to the closed list 
         */
        carray[ci].open = newseg->next;
        newseg->next = carray[ci].closed;
        carray[ci].closed = newseg;
        /*
         * Yes, we did merge, but that segment is done, so there's
         * nothing left to do here
         */
        didmerge = 0;
      }
    } else {
      didmerge = 0;
    }
  }
  return;
}

/*
 * egmap remaps the edge indices so that they work within the
 * get_global_edge function
 */
int egmap[] = { 0, 2, 1, 3 };

/*
 * f1eg and f2eg are functions that are used by get_global_edge
 * to compute global edge numbers
 */
int f1eg[]  = { 0, 0, 0, 1 };
int f2eg[]  = { 0, 0, 1, 1 };

/*
 * get_global_edge(ix, iy, ei) returns a unique numeric identifier for 
 * each edge in the grid.  For edges that overlap (i.e. the right edge
 * of a box overlaps the left edge of the box to its right), the 
 * identifier will be the same.  'ix', 'iy' are the grid indices of the
 * upper left corner of the box, and ei is the index of the edge in
 * question.  If you know a better way to do this, let me know.
 */

int get_global_edge(int ix, int iy, int ei) {

  int mi = egmap[ei];

  return 2 * (mi + ix + iy * nx + (nx - 1) * f1eg[mi]) - 3 * f2eg[mi];
}

CPoint *get_new_cpoint(void) {

  if(ncpoint == NGRAB) {
    cpbuf = (CPoint *)malloc(NGRAB * sizeof(CPoint));
    if (cpbuf == NULL) {
      fprintf(stderr, "Out of memory in get_new_cpoint\n");
      exit(1);
    }
    memset(cpbuf, 0, NGRAB * sizeof(CPoint));
    ncpoint = 0;
  }
  return &cpbuf[ncpoint++];
}

CPoint *dup_cpoint(CPoint *p) {

  CPoint *np = get_new_cpoint();

  memcpy(np, p, sizeof(CPoint));
  return np;
}

CSeg *get_new_cseg() {

  CSeg *s = (CSeg *)malloc(sizeof(CSeg));

  if(s == NULL) {
    fprintf(stderr, "Out of memory in get_new_cseg\n");
    exit(1);
  }
  s->nesting = 0;
  s->neg_depth = INT_MAX;
  s->pos_depth = INT_MAX;
  s->neg_holes = NULL;
  s->neg_next  = NULL;
  s->pos_holes = NULL;
  s->pos_next  = NULL;
  s->next  = NULL;
  s->first = NULL;
  s->last  = NULL;
  return s;
}

void reverse_cseg(CSeg *seg) {

  CPoint *s = seg->first;
  CPoint *tmp, *top;

  tmp = s->next;
  s->next = NULL;
  top = s;
  s = tmp;
  while( s->next ) {
    tmp = s->next;
    s->next = top;
    top = s;
    s = tmp;
  }
  s->next = top;
  seg->last = seg->first;
  seg->first = s;
  return;
}

CSeg *dup_cseg(CSeg *s) {

  CSeg *new = get_new_cseg();
  CPoint *ncp = dup_cpoint(s->first);
  CPoint *cp = s->first->next;

  new->first = ncp;

  while( cp ) {
    ncp->next = dup_cpoint(cp);
    ncp = ncp->next;
    cp = cp->next;
  }

  new->last = ncp;

  new->nesting = s->nesting;
  new->xmin = s->xmin;
  new->xmax = s->xmax;
  new->ymin = s->ymin;
  new->ymax = s->ymax;

  return new;
}

CSeg *get_outer_boundry_polygon() {

  CSeg *s = get_new_cseg();
  CPoint *p;

  p = get_new_cpoint();
  s->first = p;
  p->x = gedge.xmin;
  p->y = gedge.ymax;
  p->next = get_new_cpoint();
  p = p->next;
  p->x = gedge.xmax;
  p->y = gedge.ymax;
  p->next = get_new_cpoint();
  p = p->next;
  p->x = gedge.xmax;
  p->y = gedge.ymin;
  p->next = get_new_cpoint();
  p = p->next;
  p->x = gedge.xmin;
  p->y = gedge.ymin;
  p->next = get_new_cpoint();
  p = p->next;
  p->x = gedge.xmin;
  p->y = gedge.ymax;
  s->last  = p;

  return s;
}

void fixval(int i, int j, double *g) {

  int cnt = 0;
  int val = 0;
  double t;

  if( i == 0 ) {
    if( j == 0 ) {
      val += isnan(t = g[i + 1 +      j *nx]) ? 0 : (cnt++, t);
      val += isnan(t = g[i     + (j + 1)*nx]) ? 0 : (cnt++, t);
      val += isnan(t = g[i + 1 + (j + 1)*nx]) ? 0 : (cnt++, t);
    } else if( j == ny - 1 ) {
      val += isnan(t = g[i     + (j - 1)*nx]) ? 0 : (cnt++, t);
      val += isnan(t = g[i + 1 + (j - 1)*nx]) ? 0 : (cnt++, t);
      val += isnan(t = g[i + 1 +      j *nx]) ? 0 : (cnt++, t);
    } else {
      val += isnan(t = g[i     + (j - 1)*nx]) ? 0 : (cnt++, t);
      val += isnan(t = g[i + 1 +      j *nx]) ? 0 : (cnt++, t);
      val += isnan(t = g[i     + (j + 1)*nx]) ? 0 : (cnt++, t);
    }
  } else if( i == nx - 1 ) {
    if( j == 0 ) {
      val += isnan(t = g[i     + (j + 1)*nx]) ? 0 : (cnt++, t);
      val += isnan(t = g[i - 1 + (j + 1)*nx]) ? 0 : (cnt++, t);
      val += isnan(t = g[i - 1 +      j *nx]) ? 0 : (cnt++, t);
    } else if( j == ny - 1 ) {
      val += isnan(t = g[i - 1 +      j *nx]) ? 0 : (cnt++, t);
      val += isnan(t = g[i - 1 + (j - 1)*nx]) ? 0 : (cnt++, t);
      val += isnan(t = g[i     + (j - 1)*nx]) ? 0 : (cnt++, t);
    } else {
      val += isnan(t = g[i     + (j + 1)*nx]) ? 0 : (cnt++, t);
      val += isnan(t = g[i - 1 +      j *nx]) ? 0 : (cnt++, t);
      val += isnan(t = g[i     + (j - 1)*nx]) ? 0 : (cnt++, t);
    }
  } else if( j == 0 ) {
    val += isnan(t = g[i + 1 +      j *nx]) ? 0 : (cnt++, t);
    val += isnan(t = g[i     + (j + 1)*nx]) ? 0 : (cnt++, t);
    val += isnan(t = g[i - 1 +      j *nx]) ? 0 : (cnt++, t);
  } else if( j == ny - 1 ) {
    val += isnan(t = g[i - 1 +      j *nx]) ? 0 : (cnt++, t);
    val += isnan(t = g[i     + (j - 1)*nx]) ? 0 : (cnt++, t);
    val += isnan(t = g[i + 1 +      j *nx]) ? 0 : (cnt++, t);
  } else {
    val += isnan(t = g[i     + (j - 1)*nx]) ? 0 : (cnt++, t);
    val += isnan(t = g[i     + (j + 1)*nx]) ? 0 : (cnt++, t);
    val += isnan(t = g[i + 1 +      j *nx]) ? 0 : (cnt++, t);
    val += isnan(t = g[i - 1 +      j *nx]) ? 0 : (cnt++, t);
  }
  if( cnt > 0 )
    g[i + j*nx] = val / cnt;
  else
    g[i + j*nx] = 0.0;
  return;
}

double epps(void) {

  return (drand48() - 0.5) * 20.0 * contour_delta;
}

void clean_up_polygons(void) {

  int ci;
  int shorts = 0;
  int fixes  = 0;
  CInfo *clist;
  CSeg *s, *sprev;
  CPoint *first;

  for( ci = -1; ci < ncontours; ci++ ) {

    if( ci == -1 )
      clist = &contour0;
    else
      clist = &(carray[ci]);

    for( s = clist->closed; s; s = s->next ) {
      first = s->first;
      if( first->x != s->last->x || first->y != s->last->y ) {
        if(verbose)
          fprintf(stderr, "Points not identical in contour %d (%g %g)\n", 
                ci, first->x - s->last->x, first->y - s->last->y);
        s->last->x = first->x;
        s->last->y = first->y;
        fixes++;
      }
    }
    /*
     * Get rid of short contours (i.e. < 4 points)
     */
    for( s = clist->closed, sprev = NULL; s; ) {
      if( !s->first || !s->first->next || !s->first->next->next 
       || !s->first->next->next->next ) {
        if( sprev ) {
          sprev->next = s->next;
          free(s);
          s = sprev->next;
        } else {
          clist->closed = s->next;
          free(s);
          s = clist->closed;
        }
        shorts++;
      } else {
        sprev = s;
        s = s->next;
      }
    }
  }
  if( verbose )
    fprintf(stderr, "Fixed %d mismatches, discarded %d short polygons\n",
            fixes, shorts);
  return;
}

void set_xyminmax(CSeg *s) {

  CPoint *cp;

  while( s ) {
    cp = s->first;
    s->xmin = s->xmax = cp->x;
    s->ymin = s->ymax = cp->y;
    for( cp = cp->next; cp; cp = cp->next ) {
      if( cp->x < s->xmin )
        s->xmin = cp->x;
      if( cp->x > s->xmax )
        s->xmax = cp->x;
      if( cp->y < s->ymin )
        s->ymin = cp->y;
      if( cp->y > s->ymax )
        s->ymax = cp->y;
    }
    s = s->next;
  }
  return;
}

int is_contained(CSeg *out, CSeg *in, CPoint *cp) {

  CPoint *p1, *p2;
  int inside = 0;

  /*
   * A polygon can't be contained within another one if any point
   * is outside a box bounding the container
   */
  if( cp->x < out->xmin || cp->x > out->xmax ||
      cp->y < out->ymin || cp->y > out->ymax ) {
    return 0;
  }
  /*
   * Find a point on the polygon to be tested that does not lie
   * atop a point in the container -- a polygon can contain another
   * if they share a vertex, but not if they share an edge (this
   * is a rule for ArcView shapefiles).
   */
  for( p1 = out->first, p2 = out->last; p1 != out->last; p2 = p1, p1 = p1->next ) {
    if( p1->x == cp->x && p1->y == cp->y ) {
      cp = cp->next;
      if( !cp )
        return 0;
      if( ( p2->x       == cp->x && p2->y       == cp->y )
       || ( p1->next->x == cp->x && p1->next->y == cp->y )) {
        cp = 0;
        break;
      }
      p1 = out->first;
      p2 = out->last;
    }
  }
  if( !cp )
    return 0;

  /*
   * Might have a new point, so check the box again (anything
   * to short-circuit this routine...)
   */
  if( cp->x < out->xmin || cp->x > out->xmax ||
      cp->y < out->ymin || cp->y > out->ymax ) {
    return 0;
  }

  /*
   * If any vertex of a polygon is inside the other, then the whole
   * polygon is contained in the other.
   *
   * Count the edges intersecting a ray shot from the point: an odd
   * number of intersections means the point is inside the polygon,
   * even means it is outside.  I found this algorithm in some damn
   * book or another.
   */
  for( p1 = out->first, p2 = out->last; p1; p2 = p1, p1 = p1->next ) {
    if( 
        ( 
          ( p1->y <= cp->y && cp->y < p2->y ) ||
          ( p2->y <= cp->y && cp->y < p1->y ) 
        )
        &&
        ( 
          cp->x < (p2->x - p1->x) * (cp->y - p1->y) / (p2->y - p1->y) + p1->x 
        ) 
      )
      inside = !inside;
  }
  return inside ? 1 : 0;
}

int is_depression(CSeg *s) {

  CPoint *p1, *p2, *pmin;
  double y1, y2, x1, x2;
  double x, xmin, yval;

  /*
   * Find the leftmost edge of the polygon (x == xmin) at some
   * fixed value of y.  If, along that edge, the interior of the
   * polygon is to the left, the polygon is a depression, otherwise
   * it is a normal polygon.
   */

  /*
   * First, find the first two points that aren't horizontal and use
   * their midpoint as the reference x and y values
   */
  for( p1 = s->first, p2 = p1->next; p2; p1 = p2, p2 = p1->next ) {
    if( almost_equal(p1->y, p2->y, 1) ) 
      continue;
    yval = (p1->y + p2->y) / 2.0;
    xmin = (p1->x + p2->x) / 2.0;
    pmin = p1;
    break;
  }
  if( !p2 ) {
    fprintf(stderr, "Error in is_depression: can't find a point\n");
    return 0;
  }
  /* 
   * Now go through the rest of the points and look for the leftmost
   * one at yval.
   */
  for( p1 = p2, p2 = p1->next; p2; p1 = p2, p2 = p1->next ) {
    if( almost_equal(p1->y, p2->y, 1) ) {
      continue;
    } else if( p1->y < p2->y ) {
      y1 = p1->y;
      x1 = p1->x;
      y2 = p2->y;
      x2 = p2->x;
    } else {
      y2 = p1->y;
      x2 = p1->x;
      y1 = p2->y;
      x1 = p2->x;
    }
    if( yval >= y2 || yval < y1 )
      continue;
    x = x2 + (x1 - x2) * ((y2 - yval) / (y2 - y1));
    if( x < xmin ) {
      pmin = p1;
      xmin = x;
    }
  }
  if( pmin->y < pmin->next->y )
    return 0;
  else
    return 1;
}

void add_contour(double cval) {

  if( ncontours >= ncarray ) {
    if( ncarray == 0 )
      ncarray = 50;
    ncarray *= 2;
    carray = (CInfo *)realloc(carray, ncarray * sizeof(CInfo));
    if( carray == NULL ) {
      fprintf( stderr, "Can't reallocate memory for carray\n");
      exit(1);
    }
  }
  carray[ncontours].cv     = cval;
  carray[ncontours].open   = NULL;
  carray[ncontours].closed = NULL;
  ncontours++;
  return;
}

void finish_contours(void) {

  int ci;

  for( ci = 0; ci < ncontours; ci++ ) {
    connect_segments( &(carray[ci]) );
  }
  return;
}

void finish_polygons(void) {

  int ci, ciopen, ciclosed;
  CSeg *nseg, *s;

  /*
   * Find the lowest levels with open and with closed contours
   */
  ciopen   = -1;
  ciclosed = -1;
  for( ci = 0; ci < ncontours - 1; ci++ ) {
    if( ciopen == -1 && carray[ci].open != NULL )
      ciopen = ci;
    if( ciclosed == -1 && carray[ci].closed != NULL )
      ciclosed = ci;
  }

  if( ciopen == -1 ) {
    /*
     * If there are no open contours then the outside border of the map 
     * won't be included.  So we dummy up a set of points to create the 
     * enclosing polygon.
     */
    s = get_outer_boundry_polygon();
    if( ciclosed <= 0 ) {
      s->next = contour0.closed;
      contour0.closed = s;
    } else {
      s->next = carray[ciclosed - 1].closed;
      carray[ciclosed - 1].closed = s;
    }
  }

  /* 
   * Close the polygons -- this is similar to closing the contours 
   * (i.e. finish_contours()) except that contours of the next larger 
   * level can close off a polygon, so we duplicate and reverse them, 
   * (which is, essentially, the same as going to the left of the
   * *first* point instead of the last) and put them on the open list 
   * for the current level.
   */

  /*
   * The lowest level segments go to contour0...
   */
  for( s = carray[0].open; s; s = s->next ) {
    nseg = dup_cseg(s);
    reverse_cseg(nseg);
    nseg->next = contour0.open;
    contour0.open = nseg;
  }
  connect_segments(&contour0);
  /*
   * Need to set contour0.cv  -- we know it is less than the
   * first contour, but we have no way of knowing by how
   * much, so we set it to the value of the first contour
   */
  contour0.cv = carray[0].cv;

  for( ci = 0; ci < ncontours - 1; ci++ ) {
    for( s = carray[ci+1].open; s; s = s->next ) {
      nseg = dup_cseg(s);
      reverse_cseg(nseg);
      nseg->next = carray[ci].open;
      carray[ci].open = nseg;
    }
    connect_segments(&(carray[ci]));
    /*
     * Values inside a polygon are assumed to be halfway between
     * the inner and outer edges (contours)
     */
    carray[ci].cv = (carray[ci].cv + carray[ci+1].cv) / 2.0;
  }

  /*
   * Simply close the highest level contours.
   * We leave the contour value alone, because, as with the lowest
   * contour, we know the enclosed area has a higher value, but we 
   * have no way of knowing how much higher.
   */
  connect_segments(&(carray[ncontours - 1]));

  clean_up_polygons();

  set_xyminmax(contour0.closed);
  for( ci = 0; ci < ncontours; ci++ )
    set_xyminmax(carray[ci].closed);

  set_depressions();
  nest_polygons();

  return;
}

void connect_segments(CInfo *ca) {

  int fei, lei;
  double mmpt, curval;
  CSeg *this, *fseg, *prev;
  CSeg *curseg, *curprv;
  CPoint *newpt;

  while( ca->open ) {
    fseg = ca->open;
    if( (lei = grid_edge(fseg->last->x, fseg->last->y)) < 0 ) {
      /*
       * This can only happen if an open contour does not reach the
       * edge of the grid, which is an unrecoverable error and indicates
       * a bug in the program...
       */
      fprintf(stderr, "Error contouring val %f: gei1 %d "
                      "gei2 %d\n", ca->cv, fseg->first->gei,
                      fseg->last->gei);
      fprintf(stderr, "UNKNOWN edge x=%f y=%f xmin=%f xmax=%f ymin=%f "
                      "ymax=%f\n", fseg->last->x, fseg->last->y, 
                      gedge.xmin, gedge.xmax, gedge.ymin, gedge.ymax);
      exit(1);
    }
    if( lei == 0 || lei == 2 )
      mmpt = fseg->last->y;
    else
      mmpt = fseg->last->x;

    curval = DBL_MAX;
    curseg = NULL;
    curprv = NULL;
    /*
     * Inner loop over possible mates
     */
    for( this = fseg, prev = NULL; this; prev = this, this = this->next ) {
      if( (fei = grid_edge(this->first->x, this->first->y)) < 0 ) {
        fprintf(stderr, "Error contouring val %f: gei1 %d "
                        "gei2 %d\n", ca->cv, this->first->gei,
                        this->last->gei);
        fprintf(stderr, "first edge %d, last edge %d\n", fei, lei);
        fprintf(stderr, "UNKNOWN edge x=%f y=%f xmin=%f xmax=%f ymin=%f "
                        "ymax=%f\n", this->first->x, this->first->y, 
                        gedge.xmin, gedge.xmax, gedge.ymin, gedge.ymax);
        exit(1);
      }
      if( fei != lei )
        continue;
      if( lei == 0 ) {                 /* want smallest y > y0 */
        if( this->first->y < mmpt )
          continue;
        if( this->first->y < curval || curseg == NULL ) {
          curval = this->first->y;
          curseg = this;
          curprv = prev;
        }
      } else if( lei == 1 ) {                 /* want smallest x > x0 */
        if( this->first->x < mmpt )
          continue;
        if( this->first->x < curval || curseg == NULL ) {
          curval = this->first->x;
          curseg = this;
          curprv = prev;
        }
      } else if( lei == 2 ) {         /* want greatest y < y0 */
        if( this->first->y > mmpt )
          continue;
        if( this->first->y > curval || curseg == NULL ) {
          curval = this->first->y;
          curseg = this;
          curprv = prev;
        }
      } else {                         /* want greatest x < x0 */
        if( this->first->x > mmpt )
          continue;
        if( this->first->x > curval || curseg == NULL ) {
          curval = this->first->x;
          curseg = this;
          curprv = prev;
        }
      }
    }
    if( curseg ) {
      /*
       * found something?
       */
      if( curseg == fseg ) {
        /*
         * connect last point to first point, then put on closed list
         */
        newpt = dup_cpoint(fseg->first);
        fseg->last->next = newpt;
        fseg->last = newpt;
        newpt->next = NULL;

        ca->open = fseg->next;
        fseg->next = ca->closed;
        ca->closed = fseg;
      } else {
        /*
         * Connect the matching segment to the first segment,
         * remove it from the open list, then go back to the loop
         */
        fseg->last->next = curseg->first;
        fseg->last = curseg->last;
        curprv->next = curseg->next;
        free(curseg);
      }
    } else {
      /*
       * No point to connect to on this edge; make a new point 
       * located at the leftmost point on the edge and add it
       * to the segment; (the next loop iteration will view
       * this point as being on the next sequential edge -- see
       * the grid_edge() function to see how this works
       */
      newpt = get_new_cpoint();
      if( lei == 0 ) {
        newpt->x = gedge.xmin;
        newpt->y = gedge.ymax;
      } else if( lei == 1 ) {
        newpt->x = gedge.xmax;
        newpt->y = gedge.ymax;
      } else if( lei == 2 ) {
        newpt->x = gedge.xmax;
        newpt->y = gedge.ymin;
      } else {
        newpt->x = gedge.xmin;
        newpt->y = gedge.ymin;
      }
      fseg->last->next = newpt;
      fseg->last = newpt;
      newpt->next = NULL;
    }
  }
  return;
}

void set_depressions(void) {

  int ci;
  CSeg *s, *prev;

  for( s = contour0.closed; s; s = s->next ) {
    if( is_depression(s) ) {
      fprintf(stderr, "Found a depression in the lowest contour...\n");
      exit(1);
    }
  }
  for( ci = 0; ci < ncontours; ci++ ) {
    for( s = carray[ci].closed, prev = NULL; s; ) {
      if( is_depression(s) ) {
        /*
         * Depressions in this level are simply polygons at the level below...
         * Reverse the polygon:
         */
        reverse_cseg(s);
        /*
         * Remove the polygon from the current list
         */
        if( prev ) 
          prev->next = s->next;
        else
          carray[ci].closed = s->next;
        /*
         * Put the polygon on the next-lower list
         */
        if( ci == 0 ) {
          s->next = contour0.closed;
          contour0.closed = s;
        } else {
          s->next = carray[ci-1].closed;
          carray[ci-1].closed = s;
        }
        if( prev ) 
          s = prev->next;
        else
          s = carray[ci].closed;
      } else {
        prev = s;
        s = s->next;
      }
    }
  }
  return;
}

/*
 * Given the coordinates of a point on the edge of the grid, returns
 * the index of the edge upon which the point lies
 */
int grid_edge(double x, double y) {

  /*
   * Normally, you are on the x/y max/min of your coordinate, but
   * the leftmost point (looking inward from outside of the box)
   * should be considered on the following edge.
   */
  if( x == gedge.xmin ) {
    if( y == gedge.ymax )
      return 1;
    return 0;
  } else if( y == gedge.ymax ) {
    if( x == gedge.xmax )
      return 2;
    return 1;
  } else if( x == gedge.xmax ) {
    if( y == gedge.ymin )
      return 3;
    return 2;
  } else if( y == gedge.ymin ) {
    /* x == gedge.xmin has already been handled */
    return 3;
  } else {
    return -1;
  }
}

void nest_polygons(void) {

  int i, j, best;
  CInfo *clist, *clist1, *clist2, *clist3;
  CInfo dummy;
  CSeg *s, *t, *potholes;

  if( verbose ) fprintf(stderr,"Starting to compute nesting...");

  /*
   * Compute the self-nesting...
   */
  for( i = -1; i < ncontours; i++ ) {
    if( i == -1 )
      clist = &contour0;
    else
      clist = &(carray[i]);
    for( s = clist->closed; s; s = s->next ) {
      s->nesting = 0;
      for( t = clist->closed; t; t = t->next ) {
        if( t == s )
          continue;
        if( is_contained(t, s, s->first) )
          s->nesting++;
      }
    }
  }

  /*
   * Holes in any given polygon can be either polygons of one level
   * lower or one level higher.
   */
  dummy.closed = NULL;
  for( i = -1; i < ncontours; i++ ) {
    if( i == -1 ) {
      clist1 = &dummy;
      clist2 = &contour0;
    } else if( i == 0 ) {
      clist1 = &contour0;
      clist2 = &(carray[i]);
    } else {
      clist1 = &(carray[i-1]);
      clist2 = &(carray[i]);
    }
    if( i == ncontours - 1 )
      clist3 = &dummy;
    else
      clist3 = &(carray[i+1]);
    /*
     * Find the holes in the polygons.  There can be negative holes (i.e. 
     * depressions) and positive holes (i.e. the next-larger polygons).
     * This is a bit tricky: for each potential hole, we go through
     * the list of polygons we could be a hole in, looking for the
     * most deeply nested one inside of which we reside, then we
     * go through the holes, and place them on the (neg/pos)_hole list
     * of the polygon in which they form a hole (keeping only the 
     * least-nested ones.
     */
    for( j = NEG_HOLES; j <= POS_HOLES; j++ ) {
      if( j == NEG_HOLES ) {
        potholes = clist1->closed;
      } else {
        potholes = clist3->closed;
      }
      for( s = potholes; s; s = s->next ) {
        best = -1;
        s->tmp_hole = NULL;
        for( t = clist2->closed; t; t = t->next ) {
          /*
           * If I know I'm a hole inside a polygon of nesting level 2,
           * for instance, then there is no point in checking other
           * polygons at level 2 or greater...
           */
          if( best >= t->nesting )
            continue;
          if( is_contained(t, s, s->first) ) {
            best = t->nesting;
            s->tmp_hole = t;
          }
        }
      }
      for( s = potholes; s; s = s->next ) {
        if( s->tmp_hole != NULL ) {
          /*
           * Add the temp variable 't' for clarity...
           * (s->tmp_hole points to the polygon within which s is a hole)
           */
          t = s->tmp_hole;
          /*
           * If s is a lower nesting level than the current holes on the
           * list (or the list is empty) start the list with s, and set
           * the nesting level; if s is at the same nesting level, push 
           * it onto t's list of holes...
           */
          if( j == NEG_HOLES ) {
            if( s->nesting < t->neg_depth ) {
              s->neg_next  = NULL;
              t->neg_holes = s;
              t->neg_depth = s->nesting;
            } else if( s->nesting == t->neg_depth ) {
              s->neg_next = t->neg_holes;
              t->neg_holes = s;
            }
          } else {
            if( s->nesting < t->pos_depth ) {
              s->pos_next  = NULL;
              t->pos_holes = s;
              t->pos_depth = s->nesting;
            } else if( s->nesting == t->pos_depth ) {
              s->pos_next = t->pos_holes;
              t->pos_holes = s;
            }
          }
        }
      }
    }
  }
  if( verbose ) fprintf(stderr,"done.\n");
  return;
}

/*
 * Compare two doubles to determine if thay are within a
 * specified number of ULPs of one another
 */
int almost_equal(double a, double b, int maxUlps) {
  /*
   * Make sure maxUlps is non-negative and small enough that the
   * default NAN won't compare as equal to anything.
   */
  assert(sizeof(double) == sizeof(int64_t));
  assert(maxUlps > 0 && maxUlps < 4 * 1024 * 1024);
  int64_t aInt; /* = *(int64_t*)&a; */
  memcpy(&aInt, &a, sizeof(int64_t));
  /* Make aInt lexicographically ordered as a twos-complement int */
  if (aInt < 0)
    aInt = 0x8000000000000000 - aInt;
  /* Make bInt lexicographically ordered as a twos-complement int */
  int64_t bInt; /* = *(int64_t*)&b; */
  memcpy(&bInt, &b, sizeof(int64_t));
  if (bInt < 0)
    bInt = 0x8000000000000000 - bInt;
  int64_t intDiff = abs64(aInt - bInt);
  if (intDiff <= maxUlps)
    return 1;
  return 0;
}

/*
 * 64-bit integer absolute value
 */
int64_t abs64(int64_t a) {

  if (a < 0) {
    return a * (int64_t)-1;
  } else {
    return a;
  }
}

