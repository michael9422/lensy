/*
 * spectrograph.c - v1.5 (codename FlamingMarshmallow)
 *
 * A model of echelle spectrograph optics using the lensy library.
 *
 *-----------------------------------------------------------------------
 * Copyright (C) 2011 - 2015 FSF
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *-----------------------------------------------------------------------
 *
 * Compile this program with:
 *
 * 	# gcc spectrograph.c -Wall -lm -lSDL2 -o spectrograph
 *
 * Run this program with:
 *
 *	# ./spectrograph
 *
 *-----------------------------------------------------------------------
 *
 * This program is written for focusing tests of a specific camera lens.
 * For another optical system, the program must be modified and recompiled.
 *
 * The program works by first specifying what surfaces the light rays will
 * pass through. They can be flat, spherical, parabolic, cylindrical, or
 * hyperbolic. For each surface, a structure of that type is declared and the
 * surface parameters are specified (i.e. the vertex position vector, the
 * aperture diameter, and the vector from the vertex to the center of
 * curvature, or to the focus for a parabola or hyperbola).
 *
 * Then, a light ray structure is initialized, and it is traced through the
 * optics by determining where it intersects the surfaces, in order, and
 * adjusting the ray trajectory at each surface. This is done by calling
 * one of the functions:
 *
 *    lensy_intersect_paraboloid()
 *    lensy_intersect_sphere()
 *    lensy_intersect_cylinder()
 *    lensy_intersect_plane()
 *    lensy_intersect_hyperboloid()
 *
 * followed by calling one of the functions:
 *
 *    lensy_redirect_reflect()
 *    lensy_redirect_refract()
 *    lensy_redirect_diffract()
 *    lensy_redirect_impact()
 *
 * The calculations are done for three dimensional rays and surfaces. The
 * sequence that rays intersect surfaces is determined by the sequence that
 * the functions are called in the program code, and the program does not
 * attempt to determine if interference would prevent that sequence of
 * intersections.
 *
 * The program uses the SDL library for graphic display, so it must be
 * installed, including the development files.
 *
 * The program creates a image of the focal plane in a file "lensy.fits".
 * Viewing the FITS file can be done with the display program 'ds9', for
 * example.
 *
 * Send questions or contributions to me at <michael.h.williamson@gmail.com>.
 *
 *
 *
 * Change log:
 *
 *   2015-07-01  Modified for use with separate library functions & SDL2.
 *   2011-03-13  Added functions ray_cone() and ray_beam().
 *   2010-10-28  Added a line list for 3-D rendering in the graphic window.
 *
 * TO DO:
 *
 *   - Add a command line interface.
 *   - Store the parameters in a data file(s) (i.e., the surface parameters,
 *     the material index of refraction, diffraction grating parameters, etc..)
 *     This would require significant program changes.
 *   - Check malloc return value everywhere (xmalloc?).
 *   - Make the display window zoom-in, zoom-out, and 3-D.
 *   - Handle non-circular, off-center apertures.
 *   - Create a function for drawing the optic elements cross-section in the
 *     graphic window.
 */


#include <arpa/inet.h>
#include <stdio.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <search.h>
#include <stdbool.h>
#include <sys/errno.h>
#include <sys/fcntl.h>
#include <sys/mman.h>
#include <sys/un.h>
#include <sys/file.h>

#include <lensy.h>


#include <SDL2/SDL.h>

SDL_Window *optic_win, *focus_win;
SDL_Renderer *optic_ren, *focus_ren;
float optic_scale = 0.005;		// meter/pixel for display
float focus_scale = 80.0e-3/640;	// meter/pixel for display


#define	PI		3.14159265
#define	DEG2RAD		0.01745329252
#define	RAD2DEG		57.29577951

/*
 * The following structures are used for specifying the rays to be traced
 * through the optics.
 */

double in_air		= 1.000293;	// index of refraction for air
double in_vacuum	= 1.000;	// index of refraction for vacuum

//---------------------------- point source list
struct ptsource_struct {
	double p[3];			// <x, y, z> position vector (meters).
	double d[3];			// direction vector.
	double cone_dia;		// diameter of a cone of rays, (degrees)
	double cone_step;		// angular spacing of rays in the cone (degrees)
	double wavelength;		// wavelength (meters)
	char red, green, blue;		// colors for graphics.
	bool draw;
};


#define NMAX_PTS	10000

struct ptsource_struct pts[NMAX_PTS] = {
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 2.0, 490e-9,    0,   0, 255,  true  },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 0.5, 495e-9,    0,   0, 255,  false },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 0.5, 500e-9,    0,   0, 255,  false },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 0.5, 505e-9,    0,   0, 255,  false },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 0.5, 510e-9,    0,   0, 255,  false },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 0.5, 515e-9,    0,   0, 255,  false },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 0.5, 525e-9,    0,   0, 255,  false },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 0.5, 530e-9,    0,   0, 255,  false },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 0.5, 540e-9,    0,   0, 255,  false },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 0.5, 545e-9,    0,   0, 255,  false },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 0.5, 550e-9,    0, 255, 255,  false },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 0.5, 555e-9,    0, 255, 255,  false },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 0.5, 565e-9,    0, 255, 255,  false },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 0.5, 570e-9,    0, 255, 255,  false },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 0.5, 575e-9,    0, 255, 255,  false },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 0.5, 585e-9,    0, 255, 255,  false },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 0.5, 590e-9,    0, 255, 255,  false },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 0.5, 595e-9,    0, 255, 255,  false },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 0.5, 605e-9,    0, 255,   0,  false },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 0.5, 610e-9,    0, 255,   0,  false },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 0.5, 615e-9,    0, 255,   0,  false },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 0.5, 625e-9,    0, 255,   0,  false },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 0.5, 630e-9,    0, 255,   0,  false },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 0.5, 635e-9,    0, 255,   0,  false },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 0.5, 645e-9,    0, 255,   0,  false },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 0.5, 650e-9,    0, 255,   0,  false },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 0.5, 655e-9,    0, 255,   0,  false },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 0.5, 665e-9,    0, 255,   0,  false },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 0.5, 670e-9,  255, 255,   0,  false },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 0.5, 675e-9,  255, 255,   0,  false },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 0.5, 685e-9,  255, 255,   0,  false },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 0.5, 690e-9,  255, 255,   0,  false },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 0.5, 695e-9,  255, 255,   0,  false },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 0.5, 705e-9,  255, 255,   0,  false },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 }, 10.0, 0.5, 710e-9,  255,   0,   0,  false },
   { { -0.611, +0.054, 0.000 }, { +0.7790, -1.1973, 0.0000 },  0.0, 0.5, 600e-9,    0, 255,   0,  false },
   { {      0,      0,     0 }, {       0,       0,      0 },  NAN, NAN,    NAN,    0,   0,   0,  false }
};

int32_t n_pts;


LIST_HEAD(raylist);

#define	NMAX_ASS	1000

double ass[3][NMAX_ASS];	// array of spot sizes
int32_t n_ass;


LIST_HEAD(spotsize_list);

struct spotsize_struct {
	double p[3];
	int32_t n;
	double rms, rms_v[3];
	struct list_head spotsize_list;
} spotsize, *pspotsize;


ENTRY he, *phe;


LIST_HEAD(line_list);		// list for rendering in 3-D in graphic window

struct line_struct {
   double p0[3], p1[3];		// begin and end point 3-D vectors
   char red, green, blue;
   struct list_head line_list;
} le, *ple;

bool make_ll_picture;		// add to line list



/*---------------------------------------------------- plot
 * Plot one pixel on the graphics window.
 */
void plot(double x, double y, char red, char green, char blue)
{
	int x0, y0;

	x0 =  x / optic_scale + 640 / 2;
	y0 = -y / optic_scale + 480 / 2;

	SDL_SetRenderDrawColor(optic_ren, red, green, blue, 255);
	SDL_RenderDrawPoint(optic_ren, x0, y0);
}


/*---------------------------------------------------- line
 * Plot a line on the graphics window.
 *
 *   p0, p1	 	- the line segment endpoints
 *   red, green, blue	- the 8-bit color intensities
 *
 */
void line(double p0[3], double p1[3], char red, char green, char blue)
{
	int x0, y0, x1, y1;

	if (make_ll_picture) {
		ple = (struct line_struct *) malloc(sizeof(struct line_struct));
		if (ple != NULL) {
			ple->p0[0] = p0[0];
			ple->p0[1] = p0[1];
			ple->p0[2] = p0[2];
			ple->p1[0] = p1[0];
			ple->p1[1] = p1[1];
			ple->p1[2] = p1[2];
			ple->red = red;
			ple->green = green;
			ple->blue = blue;
			list_add(&(ple->line_list), &line_list);
		}
	}

	x0 =  p0[0] / optic_scale + 640 / 2;
	y0 = -p0[1] / optic_scale + 480 / 2;
	x1 =  p1[0] / optic_scale + 640 / 2;
	y1 = -p1[1] / optic_scale + 480 / 2;

	SDL_SetRenderDrawColor(optic_ren, red, green, blue, 255);
	SDL_RenderDrawLine(optic_ren, x0, y0, x1, y1);
}



//---------------------------------------------------- main
int main(int argc, char *argv[])
{
	int32_t i, j, k, i0;
	FILE *fp;
	struct lensy_ray_struct ray, *pray;
	struct list_head *pos, *pos0;
	char s80[80], s100[100], s200[200];
//	char red, green, blue;
	bool draw;
//	double wavelength, best_focus;
	double d0, d1, d2;
	double dd0, dd1, dd2;
	double w0[3], w1[3], w2[3], u0[3], u1[3], u2[3];
	struct lensy_paraboloid_struct pm;
	struct lensy_sphere_struct sp[12];
	struct lensy_cylinder_struct cyl;
	struct lensy_plane_struct pl;
	int32_t nr_hdr;
	char hdr[180][80], zeros[2880];
	int x0, y0;


	/*--------------------------------- optic elements
	 * The radii for lenses are from measurements supplied by the
	 * lens manufacturer.
	 * The spacings were from measurments of schematic diagrams, but have been
	 * modified for focus tests, and probably are now wrong in places. (The orientation
	 * of the cylindrical lens may be rotated by 90 degrees, too. I don't remember.)
	 */
	struct lensy_paraboloid_struct collimator1 = {
		{ +0.29657, -1.04348, 0.00 }, { -0.91404, +1.08932, 0.00 }, 0.6096
	};

	struct lensy_plane_struct echelleg = {
		{ -0.75515, -0.04119, 0.00 }, {  0.2927321, -0.3367498, -0.8949344 }, 0.456
	};

	struct lensy_plane_struct foldm = {
		{ -0.57391, +0.03844, 0.00 }, {  0.6427876, -0.7660444,  0.0000000 }, 0.160
	};

	struct lensy_paraboloid_struct collimator2 = {
		{ +0.35973,  -1.11213, 0.00 }, { -0.91404, +1.08932, 0.00 }, 0.6096
	};

	/*
	 * The cross dispersion grating normal vector is <cos(-19.5), sin(-19), 0>.
	 */
	struct lensy_plane_struct crossdisp = {
		{ -0.30893, +0.00000, 0.00 }, {  0.9426415, -0.3338069, 0.0000000 }, 0.260
	};

	/*
	 * The camera lens optical surfaces
	 */
	struct lensy_sphere_struct sp1[12] = {
		{ {    +0.0e-3, 0.0, 0.0 },  { +310.085e-3, 0.0, 0.0 }, 256.0e-3 },
		{ {  +37.19e-3, 0.0, 0.0 },  {  +3010.0e-3, 0.0, 0.0 }, 256.0e-3 },
		{ { +216.54e-3, 0.0, 0.0 },  { +294.167e-3, 0.0, 0.0 }, 212.0e-3 },
		{ { +224.44e-3, 0.0, 0.0 },  { +137.589e-3, 0.0, 0.0 }, 196.0e-3 },
		{ { +292.64e-3, 0.0, 0.0 },  { -279.363e-3, 0.0, 0.0 }, 196.0e-3 },
		{ { +303.61e-3, 0.0, 0.0 },  { +774.610e-3, 0.0, 0.0 }, 188.0e-3 },
		{ { +603.79e-3, 0.0, 0.0 },  { +175.180e-3, 0.0, 0.0 }, 173.0e-3 },
		{ { +663.37e-3, 0.0, 0.0 },  { -153.651e-3, 0.0, 0.0 }, 173.0e-3 },
		{ { +670.79e-3, 0.0, 0.0 },  { -348.256e-3, 0.0, 0.0 }, 173.0e-3 },
		{ { +755.55e-3, 0.0, 0.0 },  { -196.175e-3, 0.0, 0.0 },  82.0e-3 },
		{ { +760.10e-3, 0.0, 0.0 },  { -769.560e-3, 0.0, 0.0 },  92.0e-3 },
		{ { +767.36e-3, 0.0, 0.0 },  { -144.410e-3, 0.0, 0.0 },  78.0e-3 }
	};

	struct lensy_cylinder_struct cyl1 = {
		{ +776.48e-3, 0.0, 0.0 }, { -280.0e-3, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, 73.9e-3
	};

	/*
	 * The CCD
	 */
	struct lensy_ccd_struct ccd1;
	ccd1.v[0] = 783.48e-3 - 2.0e-3;
	ccd1.v[1] = 0.0;
	ccd1.v[2] = 0.0;

	ccd1.vx[0] = 0.0;
	ccd1.vx[1] = 15.0e-6;
	ccd1.vx[2] = 0.0;

	ccd1.vy[0] = 0.0;
	ccd1.vy[1] = 0.0;
	ccd1.vy[2] = 15.0e-6;

	ccd1.x_nmax = 4096;
	ccd1.y_nmax = 4096;
	lensy_init_ccd(&ccd1);
	memset(ccd1.b, 0, ccd1.b_size);


	memset(zeros, 0, sizeof(zeros));
	strcpy(s80, "");
	strcpy(s100, "");
	strcpy(s200, "");

	for (n_pts = 0; !isnan(pts[n_pts].wavelength); n_pts++);


	//------------------- Init SDL
	if (SDL_Init(SDL_INIT_EVERYTHING) < 0) {
		fprintf(stderr, "initialize SDL failed: %s\n", SDL_GetError());
		exit(-1);
	}
	atexit(SDL_Quit);

	//------------------- create a display window for optics
	optic_win = SDL_CreateWindow("Optics",
					SDL_WINDOWPOS_UNDEFINED,
					SDL_WINDOWPOS_UNDEFINED,
					640, 480, SDL_WINDOW_INPUT_GRABBED);

	optic_ren = SDL_CreateRenderer(optic_win, -1, SDL_RENDERER_ACCELERATED);
	if ((optic_win == NULL) || (optic_ren == NULL)) {
		fprintf(stderr, "create SDL window failed\n");
		exit(-1);
	}

	SDL_SetHint(SDL_HINT_RENDER_SCALE_QUALITY, "linear");
	SDL_RenderSetLogicalSize(optic_ren, 640, 480);
	SDL_SetRenderDrawColor(optic_ren, 0, 0, 0, 255);
	SDL_RenderClear(optic_ren);
	SDL_RenderPresent(optic_ren);

	//------------------- create a display window for focal plane
	focus_win = SDL_CreateWindow("Focal Plane",
					0,	// SDL_WINDOWPOS_UNDEFINED,
					0,	// SDL_WINDOWPOS_UNDEFINED,
					640, 480, 0);

	focus_ren = SDL_CreateRenderer(focus_win, -1, SDL_RENDERER_ACCELERATED);
	if ((focus_win == NULL) || (focus_ren == NULL)) {
		fprintf(stderr, "create SDL window failed\n");
		exit(-1);
	}

	SDL_SetHint(SDL_HINT_RENDER_SCALE_QUALITY, "linear");
	SDL_RenderSetLogicalSize(focus_ren, 640, 480);
	SDL_SetRenderDrawColor(focus_ren, 0, 0, 0, 255);
	SDL_RenderClear(focus_ren);
	SDL_RenderPresent(focus_ren);

	make_ll_picture = true;
	draw = true;
	i0 = 0;

	//--- focus adjusment parameters
//	best_focus = 1e6;
	dd0 =  0.0e-4;
	dd1 =  0.0e-3;
	dd2 =  0.0e-3;

ray_trace_loop:
	SDL_SetRenderDrawColor(optic_ren, 0, 0, 0, 255);
	SDL_RenderClear(optic_ren);

	SDL_SetRenderDrawColor(focus_ren, 0, 0, 0, 255);
	SDL_RenderClear(focus_ren);

	memcpy(&pm, &collimator2, sizeof(pm));
	memcpy(&sp, &sp1, sizeof(sp));
	memcpy(&cyl, &cyl1, sizeof(cyl));
	memcpy(&pl, &ccd1.p, sizeof(pl));

	sp[9].v[0] 	+= dd0;
	sp[10].v[0] 	+= dd0;
	sp[11].v[0] 	+= dd0 + dd1;
	if (sp[11].v[0] - sp[10].v[0] < 4.0e-3) return (-1);	// minimum spacing

	cyl.v[0] 	+= dd0 + dd1;
	pl.v[0] 	+= dd0 + dd1 + dd2;
	if (pl.v[0] - cyl.v[0] < 4.0e-3) return (-1);		// minimum spacing

#if 0
	/*------------------------ realignment test
	 * Test tilting the collimator2 slightly. The cross-dispersion grating needs
	 * to be rotated as well, for this test.
	 */
	w0[0] = pm.f[0];
	w0[1] = pm.f[1];
	w0[2] = pm.f[2];

	d0 = lensy_mag3(pm.f);
	d1 = 0.01;
	w0[0] = d1*w0[0]/d0;
	w0[1] = d1*w0[1]/d0;
	w0[2] = d1*w0[2]/d0;

	pm.f[0] +=  w0[1];
	pm.f[1] += -w0[0];
	pm.f[2] +=  0.0;
#endif

	n_ass = 0;

	/*
	 * Create a cone of rays, in a linked list.
	 * The raylist is generated and erased for each iteration of the ray
	 * trace loop.
	 *
	 * The source for the spectrograph is a polished optic fiber. To model
	 * that accurately here would require multiple cones with slightly
	 * different positions (ray.p[]) to account for the fiber diameter.
	 * A simplification is to use a point source.
	 */
	j = 0;

	for (i = 0; i < n_pts; i++) {
		ray.p[0] = pts[i].p[0];
		ray.p[1] = pts[i].p[1];
		ray.p[2] = pts[i].p[2];

		ray.d[0] = pts[i].d[0];
		ray.d[1] = pts[i].d[1];
		ray.d[2] = pts[i].d[2];

		ray.wavelength	= pts[i].wavelength;
		ray.red		= pts[i].red;
		ray.green	= pts[i].green;
		ray.blue	= pts[i].blue;

		j += lensy_cone(&raylist, &ray, pts[i].cone_dia, pts[i].cone_step);
	}
	printf("rays %d\n", j);

	//--------------- collimator1
	list_for_each_safe(pos, pos0, &raylist) {
		pray = list_entry(pos, struct lensy_ray_struct, raylist);

		i = lensy_intersect_paraboloid(pray, &collimator1, w0, w1);
		if (draw)
			line(pray->p, w0, pray->red, pray->green, pray->blue);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}
		lensy_redirect_reflect(pray, w0, w1);
	}

	/*---------------- echelle grating
	 * For the diffraction calculation, the grating ruling direction is
	 * determined from both the normal vector to the intersect plane, and
	 * the normal vector to the rulings, i.e., the vector given here for
	 * 'diffract' is not, and is not required to be in the surface plane
	 * of the grating.
	 */
	w2[0] = +0.65606;
	w2[1] = -0.75471;
	w2[2] =     0.00;

	d0 = lensy_mag3(w2);
	d1 = 1.901141e-5;		// the ruling spacing

	w2[0] = d1*w2[0]/d0;
	w2[1] = d1*w2[1]/d0;
	w2[2] = d1*w2[2]/d0;

	list_for_each_safe(pos, pos0, &raylist) {
		pray = list_entry(pos, struct lensy_ray_struct, raylist);

		i = lensy_intersect_plane(pray, &echelleg, w0, w1);
		if (draw)
			line(pray->p, w0, pray->red, pray->green, pray->blue);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}

		/*--------------------------------
		 * The original ray is deleted and here a bunch of new rays are created
		 * and added to the raylist for a range of reflected orders off of the
		 * grating.
		 */
		memcpy(&ray, pray, sizeof(struct lensy_ray_struct));
		list_del(pos);
		free(pray);

		for (j = 40; j < 100; j++) {
			pray = (struct lensy_ray_struct *) malloc (sizeof(struct lensy_ray_struct));
			if (pray == NULL) continue;

			memcpy(pray, &ray, sizeof(struct lensy_ray_struct));
			sprintf (s100, "%d", j);
			i = sizeof(pray->pathkey) - strlen(pray->pathkey) - 1;
			strncat (pray->pathkey, s100, i);

			i = lensy_redirect_diffract(pray, w0, w1, w2, pray->wavelength, pray->wavelength, j);
			if (i < 0) {
				free(pray);
				continue;
			}
			list_add(&(pray->raylist), &raylist);
		}
	}

	//--------------- collimator1 again
	list_for_each_safe(pos, pos0, &raylist) {
		pray = list_entry(pos, struct lensy_ray_struct, raylist);

		i = lensy_intersect_paraboloid(pray, &collimator1, w0, w1);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}
		if (draw)
			line(pray->p, w0, pray->red, pray->green, pray->blue);
		lensy_redirect_reflect(pray, w0, w1);
	}

	//--------------- fold mirror
	list_for_each_safe(pos, pos0, &raylist) {
		pray = list_entry(pos, struct lensy_ray_struct, raylist);

		i = lensy_intersect_plane(pray, &foldm, w0, w1);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}
		if (draw)
			line(pray->p, w0, pray->red, pray->green, pray->blue);
		lensy_redirect_reflect(pray, w0, w1);
	}

	//--------------- collimator2
	list_for_each_safe(pos, pos0, &raylist) {
		pray = list_entry(pos, struct lensy_ray_struct, raylist);

		i = lensy_intersect_paraboloid(pray, &pm, w0, w1);	// &collimator2
		if (draw)
			line(pray->p, w0, pray->red, pray->green, pray->blue);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}
		lensy_redirect_reflect(pray, w0, w1);
	}

	//---------------- cross-dispersion grating
	w2[0] = +0.00000;
	w2[1] = -1.00000;
	w2[2] =     0.00;

	d0 = lensy_mag3(w2);
	d1 = 4.0e-6;		// the ruling spacing

	w2[0] = d1*w2[0]/d0;
	w2[1] = d1*w2[1]/d0;
	w2[2] = d1*w2[2]/d0;

	list_for_each_safe(pos, pos0, &raylist) {
		pray = list_entry(pos, struct lensy_ray_struct, raylist);

		i = lensy_intersect_plane(pray, &crossdisp, w0, w1);
		if (draw)
			line(pray->p, w0, pray->red, pray->green, pray->blue);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}
		i = lensy_redirect_diffract(pray, w0, w1, w2, pray->wavelength, pray->wavelength, +1);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}
	}

	//------------- camera lens
	list_for_each_safe(pos, pos0, &raylist) {
		pray = list_entry(pos, struct lensy_ray_struct, raylist);

		i = lensy_intersect_sphere(pray, &sp[0], w0, w1);
		if (draw)
			line(pray->p, w0, pray->red, pray->green, pray->blue);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}
		d0 = in_air / lensy_index_of_refraction(pray->wavelength, CaF2);
		i = lensy_redirect_refract(pray, w0, w1, d0);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}

		i = lensy_intersect_sphere(pray, &sp[1], w0, w1);
		if (draw)
			line(pray->p, w0, pray->red, pray->green, pray->blue);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}
		d0 = lensy_index_of_refraction(pray->wavelength, CaF2) / in_air;
		i = lensy_redirect_refract(pray, w0, w1, d0);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}

		i = lensy_intersect_sphere(pray, &sp[2], w0, w1);
		if (draw) line(pray->p, w0, pray->red, pray->green, pray->blue);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}
		d0 = in_air / lensy_index_of_refraction(pray->wavelength, tsu2);
		i = lensy_redirect_refract(pray, w0, w1, d0);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}

		i = lensy_intersect_sphere(pray, &sp[3], w0, w1);
		if (draw)
			line(pray->p, w0, pray->red, pray->green, pray->blue);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}
		d0 = lensy_index_of_refraction(pray->wavelength, tsu2) /
			lensy_index_of_refraction(pray->wavelength, CaF2);
		i = lensy_redirect_refract(pray, w0, w1, d0);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}

		i = lensy_intersect_sphere(pray, &sp[4], w0, w1);
		if (draw)
			line(pray->p, w0, pray->red, pray->green, pray->blue);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}
		d0 = lensy_index_of_refraction(pray->wavelength, CaF2) /
			lensy_index_of_refraction(pray->wavelength, tsu4);
		i = lensy_redirect_refract(pray, w0, w1, d0);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}

		i = lensy_intersect_sphere(pray, &sp[5], w0, w1);
		if (draw)
			line(pray->p, w0, pray->red, pray->green, pray->blue);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}
		d0 = lensy_index_of_refraction(pray->wavelength, tsu4) / in_air;
		i = lensy_redirect_refract(pray, w0, w1, d0);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}

		i = lensy_intersect_sphere(pray, &sp[6], w0, w1);
		if (draw)
			line(pray->p, w0, pray->red, pray->green, pray->blue);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}
		d0 = in_air / lensy_index_of_refraction(pray->wavelength, tsu5);
		i = lensy_redirect_refract(pray, w0, w1, d0);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}

		i = lensy_intersect_sphere(pray, &sp[7], w0, w1);
		if (draw)
			line(pray->p, w0, pray->red, pray->green, pray->blue);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}
		d0 = lensy_index_of_refraction(pray->wavelength, tsu5) /
			lensy_index_of_refraction(pray->wavelength, tsu6);
		i = lensy_redirect_refract(pray, w0, w1, d0);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}

		i = lensy_intersect_sphere(pray, &sp[8], w0, w1);
		if (draw)
			line(pray->p, w0, pray->red, pray->green, pray->blue);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}
		d0 = lensy_index_of_refraction(pray->wavelength, tsu6)/in_air;
		i = lensy_redirect_refract(pray, w0, w1, d0);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}

		i = lensy_intersect_sphere(pray, &sp[9], w0, w1);
		if (draw)
			line(pray->p, w0, pray->red, pray->green, pray->blue);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}
		d0 = in_air / lensy_index_of_refraction(pray->wavelength, tsu7);
		i = lensy_redirect_refract(pray, w0, w1, d0);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}

		i = lensy_intersect_sphere(pray, &sp[10], w0, w1);
		if (draw)
			line(pray->p, w0, pray->red, pray->green, pray->blue);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}
		d0 = lensy_index_of_refraction(pray->wavelength, tsu7) / in_air;
		i = lensy_redirect_refract(pray, w0, w1, d0);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}

		i = lensy_intersect_sphere(pray, &sp[11], w0, w1);
		if (draw)
			line(pray->p, w0, pray->red, pray->green, pray->blue);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}
		d0 = in_air / lensy_index_of_refraction(pray->wavelength, fsilica);
		i = lensy_redirect_refract(pray, w0, w1, d0);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}

		i = lensy_intersect_cylinder(pray, &cyl, w0, w1);
		if (draw)
			line(pray->p, w0, pray->red, pray->green, pray->blue);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}
		d0 = lensy_index_of_refraction(pray->wavelength, fsilica) / in_vacuum;
		i = lensy_redirect_refract(pray, w0, w1, d0);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}

		i = lensy_intersect_plane(pray, &ccd1.p, w0, w1);
		if (draw)
			line(pray->p, w0, pray->red, pray->green, pray->blue);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}
		lensy_redirect_impact (pray, w0, w1);
	}

	//------- add the impact positions to the focal plane picture
	list_for_each_safe(pos, pos0, &raylist) {
		pray = list_entry(pos, struct lensy_ray_struct, raylist);

		w1[0] = pray->p[0] - ccd1.v[0];
		w1[1] = pray->p[1] - ccd1.v[1];
		w1[2] = pray->p[2] - ccd1.v[2];

		i = floor(lensy_inner3(w1, ccd1.vx) / lensy_inner3(ccd1.vx, ccd1.vx));
		j = floor(lensy_inner3(w1, ccd1.vy) / lensy_inner3(ccd1.vy, ccd1.vy));

		i += ccd1.x_nmax/2;
		j += ccd1.y_nmax/2;

		if ((i >= 0) && (i < ccd1.x_nmax) &&
		    (j >= 0) && (j < ccd1.y_nmax)) {
			k = ccd1.b[j * ccd1.x_nmax + i];
			k += (k < 65000) ? 100 : 0;
			ccd1.b[j*ccd1.x_nmax + i] = k;

			x0 =  i * lensy_mag3(ccd1.vx) / focus_scale;
			y0 =  j * lensy_mag3(ccd1.vy) / focus_scale;

			SDL_SetRenderDrawColor(focus_ren, pray->red, pray->green, pray->blue, 255);
			SDL_RenderDrawPoint(focus_ren, x0, y0);
		}
	}


	/*---------------- perform the average spot size calculation
	 * The idea behind using the hash table and 'pathkey' is to identify
	 * rays with identical parameters (start position, wavelength, grating
	 * reflection order), for the purpose of calculating a spot size (for
	 * example, for a spectrograph).
	 */
	hcreate(100);

	list_for_each_safe(pos, pos0, &raylist) {
		pray = list_entry(pos, struct lensy_ray_struct, raylist);

		he.key = (char *) pray->pathkey;
		phe = hsearch(he, FIND);
		if (phe == NULL) {
			pspotsize = (struct spotsize_struct *) malloc(sizeof(struct spotsize_struct));
			pspotsize->p[0] 	= pray->p[0];
			pspotsize->p[1] 	= pray->p[1];
			pspotsize->p[2] 	= pray->p[2];
			pspotsize->n 		= 1;
			pspotsize->rms_v[0] 	= 0.0;
			pspotsize->rms_v[1] 	= 0.0;
			pspotsize->rms_v[2] 	= 0.0;
			pspotsize->rms 		= 0.0;
			list_add(&(pspotsize->spotsize_list), &spotsize_list);
			he.key = (char *) pray->pathkey;
			he.data = (void *) pspotsize;
			if (hsearch (he, ENTER) == NULL) {
				fprintf (stderr, "hash action ENTER failed: \"%s\"\n", he.key);
				exit(-1);
			}
		} else {
			pspotsize = (struct spotsize_struct *) phe->data;
			pspotsize->p[0] += pray->p[0];
			pspotsize->p[1] += pray->p[1];
			pspotsize->p[2] += pray->p[2];
			pspotsize->n++;
		}
	}

	//---------------- find the centroid of the spots
	list_for_each_safe(pos, pos0, &spotsize_list) {
		pspotsize = list_entry(pos, struct spotsize_struct, spotsize_list);
		pspotsize->p[0] /= pspotsize->n;
		pspotsize->p[1] /= pspotsize->n;
		pspotsize->p[2] /= pspotsize->n;
	}

	//------------- calculate the sum of squares
	list_for_each_safe(pos, pos0, &raylist) {
		pray = list_entry(pos, struct lensy_ray_struct, raylist);

		he.key = (char *) pray->pathkey;
		phe = hsearch(he, FIND);
		if (phe == NULL) {
			fprintf(stderr, "hash action FIND failed: \"%s\"\n", he.key);
//			exit(-1);
		} else {
			pspotsize = (struct spotsize_struct *) phe->data;

			w0[0] = pray->p[0] - pspotsize->p[0];
			w0[1] = pray->p[1] - pspotsize->p[1];
			w0[2] = pray->p[2] - pspotsize->p[2];

			pspotsize->rms_v[0] += w0[0] * w0[0];
			pspotsize->rms_v[1] += w0[1] * w0[1];
			pspotsize->rms_v[2] += w0[2] * w0[2];

			pspotsize->rms += lensy_inner3(w0, w0);
		}
	}

	//---------------- calculate the RMS
	list_for_each_safe(pos, pos0, &spotsize_list) {
		pspotsize = list_entry(pos, struct spotsize_struct, spotsize_list);

		pspotsize->rms_v[0] = sqrt(pspotsize->rms_v[0] / pspotsize->n);
		pspotsize->rms_v[1] = sqrt(pspotsize->rms_v[1] / pspotsize->n);
		pspotsize->rms_v[2] = sqrt(pspotsize->rms_v[2] / pspotsize->n);

		pspotsize->rms = sqrt(pspotsize->rms / pspotsize->n);
	}
	hdestroy ();

	//------------- store in an array and free the spotsize_list
	list_for_each_safe(pos, pos0, &spotsize_list) {
		pspotsize = list_entry(pos, struct spotsize_struct, spotsize_list);

		if (pspotsize->n <= 1) continue;

		if (n_ass < NMAX_ASS) {
			ass[0][n_ass] = pspotsize->rms_v[0];
			ass[1][n_ass] = pspotsize->rms_v[1];
			ass[2][n_ass] = pspotsize->rms_v[2];
			n_ass++;
		}
		// pspotsize->rms is not used?

		list_del(pos);
		free(pspotsize);
	}

	//---------------- free the raylist
	list_for_each_safe(pos, pos0, &raylist) {
		pray = list_entry(pos, struct lensy_ray_struct, raylist);
		list_del(pos);
		free(pray);
	}

	if (make_ll_picture) {
		//------------- write a FITS file showing the focal plane
		fp = fopen("lensy.fits", "w");
		if (fp == NULL) {
			fprintf(stderr, "fopen FITS file failed\n");
			exit(0);
		}

		i = 0;
		memset (hdr, ' ', sizeof(hdr));
		snprintf (hdr[i++], 80, "SIMPLE  = %20s%-50s", "T", "");
		snprintf (hdr[i++], 80, "BITPIX  = %20d%-50s", 16, "");
		snprintf (hdr[i++], 80, "NAXIS   = %20d%-50s", 2, "");
		snprintf (hdr[i++], 80, "NAXIS1  = %20d%-50s", ccd1.x_nmax, "");
		snprintf (hdr[i++], 80, "NAXIS2  = %20d%-50s", ccd1.y_nmax, "");
		snprintf (hdr[i++], 80, "ORIGIN  = %-70s", "'lensy'");
		snprintf (hdr[i++], 80, "BZERO   = %20.0f%-50s", 32768.0, "");
		snprintf (hdr[i++], 80, "BSCALE  = %20.0f%-50s", 1.0, "");
		snprintf (hdr[i++], 80, "%-80s", "END");

		nr_hdr = (i-1)/36 + 1;
		while (--i >= 0) hdr[i][79] = ' ';
		fwrite (hdr, 1, nr_hdr*36*80, fp);

		// swap bytes/sign for FITS file format (big-endian)
		for (i = 0; i < (ccd1.x_nmax * ccd1.y_nmax); i++)
			ccd1.b[i] = htons(ccd1.b[i] ^ 0x8000);
		fwrite (ccd1.b, 1, ccd1.b_size, fp);
		fwrite (zeros, 1, (2880 - ccd1.b_size%2880) % 2880, fp);
		fclose (fp);

		// undo byte swap/sign
		for (i = 0; i < (ccd1.x_nmax * ccd1.y_nmax); i++)
			ccd1.b[i] = ntohs(ccd1.b[i]) ^ 0x8000;
	}

	//------ draw x, y axes
	w0[0] = -1.0;
	w0[1] =  0.0;
	w0[2] =  0.0;

	w1[0] = +1.0;
	w1[1] =  0.0;
	w1[2] =  0.0;

	line(w0, w1, 255, 0, 255);

	w0[0] =  0.0;
	w0[1] = -1.0;
	w0[2] =  0.0;

	w1[0] =  0.0;
	w1[1] = +1.0;
	w1[2] =  0.0;

	line(w0, w1, 255, 255, 0);

	w0[0] =  0.0;
	w0[1] =  0.0;
	w0[2] = -1.0;

	w1[0] =  0.0;
	w1[1] =  0.0;
	w1[2] = +1.0;

	line(w0, w1, 0, 255, 255);



   //------------- draw the surfaces of the optics in the graphic window
   for (d0 = -1.5; d0 < +0.0; d0 += 0.2) {
      ray.p[0] 	=  0.0;		// ray starting position vector in meters
      ray.p[1] 	=  d0;
      ray.p[2] 	=  0.0;

      ray.d[0]	=  1.0;		// ray direction vector
      ray.d[1] 	= -1.0;
      ray.d[2] 	=  0.0;

      i = lensy_intersect_paraboloid(&ray, &collimator1, w0, w1);
      j = (i == 0) ? 255 : 100;
      line(w0, w0, j, j, j);

      i = lensy_intersect_paraboloid(&ray, &pm, w0, w1);		// &collimator2
      j = (i == 0) ? 255 : 100;
      line(w0, w0, j, j, j);


      ray.d[0]	= -1.0;		// ray direction vector
      ray.d[1] 	=  1.0;
      ray.d[2] 	=  0.0;

      i = lensy_intersect_plane(&ray, &echelleg, w0, w1);
      j = (i == 0) ? 255 : 100;
      line(w0, w0, j, j, j);

      i = lensy_intersect_plane(&ray, &foldm, w0, w1);
      j = (i == 0) ? 255 : 100;
      line(w0, w0, j, j, j);

      i = lensy_intersect_plane(&ray, &crossdisp, w0, w1);
      j = (i == 0) ? 255 : 100;
      line(w0, w0, j, j, j);
   }

   for (d0 = -30.0e-2; d0 < +30.0e-2; d0 += 0.2) {
      ray.p[0] 	=  -1.0;	// ray starting position vector in meters
      ray.p[1] 	=  d0;
      ray.p[2] 	=  0.0;

      ray.d[0]	=  1.0;		// ray direction vector
      ray.d[1] 	=  0.0;
      ray.d[2] 	=  0.0;
      for (k = 0; k < 12; k++) {
         i = lensy_intersect_sphere(&ray, &sp[k], w0, w1);
         j = (i < 0) ? 64 : 255;
         if (i != -2) line(w0, w0, j, j, j);
      }

      i = lensy_intersect_cylinder(&ray, &cyl, w0, w1);
      j = (i < 0) ? 64 : 255;
      if (i != -2) line(w0, w0, j, j, j);

      i = lensy_intersect_plane(&ray, &pl, w0, w1);
      j = (i < 0) ? 64 : 255;
      if (i != -2) line(w0, w0, j, j, j);
   }


	//-------------- show the ray trace picture
	SDL_RenderPresent(optic_ren);
	SDL_RenderPresent(focus_ren);


	w0[0] = 0.0;
	w0[1] = 0.0;
	w0[2] = 0.0;
	for (i = 0; i < n_ass; i++) {
		w0[0] += ass[0][i];
		w0[1] += ass[1][i];
		w0[2] += ass[2][i];
	}
	if (n_ass > 0) {
		w0[0] /= n_ass;
		w0[1] /= n_ass;
		w0[2] /= n_ass;
	}

	printf ("spotsize x = %5.0lfum, y = %5.0lfum, z = %5.0lfum\n",
		 2 * w0[0] * 1e6, 2 * w0[1] * 1e6, 2 * w0[2] * 1e6);
	SDL_PumpEvents();
	SDL_Delay(1);

	if (i0++ < 5) {
		dd0 += 0.1e-3;
		// do not add to list subsequent calls to the line() function
		make_ll_picture = false;
		goto ray_trace_loop;
	}
	SDL_Delay(500);
/*
	printf ("press enter...");
	fflush (stdout);
	fgets (s100, 99, stdin);
*/

	//---------------------------- show the 3-D line_list picture

	/*
	 * Create two angles d0, and d1 for the primed coordinates, which will
	 * be used to project the lines from the line list onto. Also multiply
	 * by a shrinking scale factor, d2.
	 */
	d2 = 1.0;
	for (d0 = 0.0; d0 < PI / 2; d0 += 3.0 * PI / 180.0) {
		d1 =  1.5 * d0;
//		d2 *= 0.99;
		optic_scale *= 0.98;

		SDL_SetRenderDrawColor(optic_ren, 0, 0, 0, 255);
		SDL_RenderClear(optic_ren);

		//------ create unit vectors for a basis for primed coordinates
		u0[0] = cos(d1) * cos(d0);
		u0[1] = cos(d1) * sin(d0);
		u0[2] = sin(d1);

		u1[0] = -sin(d0);
		u1[1] = cos(d0);
		u1[2] = 0.0;

		lensy_cross3(u0, u1, u2);

		//--------- draw the list of lines projected onto two of the primed axes
		list_for_each_safe(pos, pos0, &line_list) {
			ple = list_entry(pos, struct line_struct, line_list);

			w0[0] = d2 * lensy_inner3(ple->p0, u0);
			w0[1] = d2 * lensy_inner3(ple->p0, u1);
			w0[2] = d2 * lensy_inner3(ple->p0, u2);

			w1[0] = d2 * lensy_inner3(ple->p1, u0);
			w1[1] = d2 * lensy_inner3(ple->p1, u1);
			w1[2] = d2 * lensy_inner3(ple->p1, u2);

			line(w0, w1, ple->red, ple->green, ple->blue);
		}

		SDL_RenderPresent(optic_ren);
		SDL_PumpEvents();
		SDL_Delay(2);
	}

	//-------------------------- free the line_list
	list_for_each_safe(pos, pos0, &line_list) {
		ple = list_entry(pos, struct line_struct, line_list);
		list_del(pos);
		free(ple);
	}

	exit(0);
}
