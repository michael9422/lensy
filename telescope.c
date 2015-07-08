/*
 * telescope.c - v1.5 (codename FlamingMarshmallow)
 *
 * A model of telescope optics using the lensy library functions.
 *
 *-----------------------------------------------------------------------
 * Copyright (C) 2013 - 2015 FSF
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
 * 	# gcc telescope.c -Wall -llensy -lm -lSDL2 -o telescope
 *
 * Run this program with:
 *
 *	# ./telescope
 *
 *-----------------------------------------------------------------------
 *
 * This program is written for focusing tests of a specific telescope.
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
 * optics by calculating where it intersects the surfaces, in order, and
 * adjusting the ray direction at each surface. This is done by calling
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
 * The program creates an image of the focal plane in a file "lensy.fits".
 * Viewing the FITS file can be done with the display program 'ds9', for
 * example.
 *
 * Send questions or contributions to me at <michael.h.williamson@gmail.com>.
 *
 * History:
 *
 *  2015-06-11  Change to use SDL2 library
 *
 */


#include <arpa/inet.h>
#include <ctype.h>
#include <fcntl.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <search.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/errno.h>
#include <sys/fcntl.h>
#include <sys/file.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/un.h>
#include <time.h>
#include <unistd.h>

#include <lensy.h>


#include <SDL2/SDL.h>

SDL_Window *optic_win, *focus_win;
SDL_Renderer *optic_ren, *focus_ren;
float optic_scale = 0.007;		// meter/pixel for display
float focus_scale = 4.0e-3/640;		// meter/pixel for display

/*
 * The following structures are used for specifying the rays to be traced
 * through the optics.
 */

double in_air		= 1.000293;	// index of refraction for air
double in_vacuum	= 1.000;	// index of refraction for vacuum


LIST_HEAD(raylist);


#define	NMAX_ASS	1000

double ass[3][NMAX_ASS];		// array of spot sizes for x, y, z axes
int32_t n_ass;

LIST_HEAD(spotsize_list);

struct spotsize_struct {
	double p[3];
	int32_t n;
	double rms, rms_v[3];
	struct list_head spotsize_list;
} spotsize, *pspotsize;


ENTRY he, *phe;


LIST_HEAD(line_list);			// list for showing ray paths

struct line_struct {
	double p0[3], p1[3];		// begin and end point 3-D vectors
	char red, green, blue;
	struct list_head line_list;
} le, *ple;

bool make_ll_picture;			// make a line list picture




/*---------------------------------------------------- plot
 * Plot one pixel on the graphics window.
 */
void plot(double x, double y, char red, char green, char blue)
{
	int x0, y0;

	x0 =  x / optic_scale + 640 / 2 - 200;
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

	x0 =  p0[0] / optic_scale + 640 / 2 - 200;
	y0 = -p0[1] / optic_scale + 480 / 2;
	x1 =  p1[0] / optic_scale + 640 / 2 - 200;
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
//	double wavelength;
	double d0, d1, d2;
	double w0[3], w1[3], w2[3];
	double u0[3], u1[3], u2[3];
	int32_t nr_hdr;
	char hdr[180][80], zeros[2880];
	int x0, y0;

	/*----------------------------------------- surface parameters
	 * Define the optical elements of the system.
	 */
	struct lensy_paraboloid_struct primary = {
		{  0, 0, 0 }, {  +3.0432, 0, 0 }, 2.0
	};
	struct lensy_hyperboloid_struct secondary = {
		{  2.6314 + 0.3e-3, 0,  0 }, {  -0.9007, 0, 0 }, 1.4577, 0.279
	};
	struct lensy_plane_struct flat1 = {
		{ 0.420 + 66.0e-3, 0, 0 }, { +1.0, 0, 0 }, 50.0e-3
	};
	struct lensy_sphere_struct sphere1 = {
		{ 0.420 + 63.0e-3, 0, 0 }, { -100.0e-3, 0, 0 }, 50.0e-3
	};
	struct lensy_plane_struct cube0 = {
		{ 0.420 + 15.0e-3 + 30e-3, 0, 0 }, { +1.0, 0, 0 }, 30.0e-3
	};
	struct lensy_plane_struct cube1 = {
		{ 0.420 + 15.0e-3, 0, 0 }, { +1.0, 0, 0 }, 30.0e-3
	};

	struct lensy_ccd_struct ccd1;

	ccd1.v[0] = 0.420;
	ccd1.v[1] = 0.0;
	ccd1.v[2] = 0.0;

	ccd1.vx[0] = 0.0;
	ccd1.vx[1] = 0.0;
	ccd1.vx[2] = -4.0e-6;

	ccd1.vy[0] = 0.0;
	ccd1.vy[1] = 4.0e-6;
	ccd1.vy[2] = 0.0;

	ccd1.x_nmax = 1000;
	ccd1.y_nmax = 1000;
	lensy_init_ccd(&ccd1);
	memset(ccd1.b, 0, ccd1.b_size);

	memset(zeros, 0, sizeof(zeros));
	strcpy(s80, "");
	strcpy(s100, "");
	strcpy(s200, "");

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

ray_trace_loop:
	SDL_SetRenderDrawColor(optic_ren, 0, 0, 0, 255);
	SDL_RenderClear(optic_ren);

	SDL_SetRenderDrawColor(focus_ren, 0, 0, 0, 255);
	SDL_RenderClear(focus_ren);

	n_ass = 0;

	/*
	 * Create a beam of rays, in a linked list.
	 * The raylist is generated and erased for each iteration of the ray
	 * trace loop.
	 */
	ray.p[0] =  1.0;
	ray.p[1] =  0.0;
	ray.p[2] =  0.0;
	ray.d[0] = -1.0;
	ray.d[1] =  0.0; // + 2.0e-4;
	ray.d[2] =  0.0;

	i = 0;

	ray.wavelength	= 800e-9;
	ray.red		= 200;
	ray.green	= 40;
	ray.blue	= 0;
	i += lensy_beam(&raylist, &ray, 2.1, 0.07);

	ray.wavelength	= 600e-9;
	ray.red		= 40;
	ray.green	= 200;
	ray.blue	= 0;
	i += lensy_beam(&raylist, &ray, 2.1, 0.07);

	ray.wavelength	= 400e-9;
	ray.red		= 0;
	ray.green	= 40;
	ray.blue	= 200;
	i += lensy_beam(&raylist, &ray, 2.1, 0.07);

	printf("rays in the beam %d\n", i);

	//---- eliminate rays in the center (for the central hole)
	list_for_each_safe(pos, pos0, &raylist) {
		pray = list_entry(pos, struct lensy_ray_struct, raylist);

		d0 = hypot(pray->p[1], pray->p[2]);
		if (d0 < 0.254) {
			list_del(pos);
			free(pray);
			continue;
		}
	}

	//--------------- primary mirror
	list_for_each_safe(pos, pos0, &raylist) {
		pray = list_entry(pos, struct lensy_ray_struct, raylist);

		i = lensy_intersect_paraboloid(pray, &primary, w0, w1);
		if (draw)
			line(pray->p, w0, pray->red, pray->green, pray->blue);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}
		lensy_redirect_reflect(pray, w0, w1);
	}

	//--------------- secondary mirror
	list_for_each_safe(pos, pos0, &raylist) {
		pray = list_entry(pos, struct lensy_ray_struct, raylist);

		i = lensy_intersect_hyperboloid(pray, &secondary, w0, w1);
		if (draw)
			line(pray->p, w0, pray->red, pray->green, pray->blue);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}
		lensy_redirect_reflect(pray, w0, w1);
	}

	//----------------- flat1
	list_for_each_safe(pos, pos0, &raylist) {
		pray = list_entry(pos, struct lensy_ray_struct, raylist);

		i = lensy_intersect_plane(pray, &flat1, w0, w1);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}
		if (draw)
			line(pray->p, w0, pray->red, pray->green, pray->blue);
		d0 = in_air / lensy_index_sellmeier(pray->wavelength, &N_BK7);
		i = lensy_redirect_refract(pray, w0, w1, d0);
	}

	//----------------- sphere1
	list_for_each_safe(pos, pos0, &raylist) {
		pray = list_entry(pos, struct lensy_ray_struct, raylist);

		i = lensy_intersect_sphere(pray, &sphere1, w0, w1);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}
		if (draw)
			line(pray->p, w0, pray->red, pray->green, pray->blue);
		d0 = lensy_index_sellmeier(pray->wavelength, &N_BK7) / in_air;
		i = lensy_redirect_refract(pray, w0, w1, d0);
	}

	//----------------- cube0
	list_for_each_safe(pos, pos0, &raylist) {
		pray = list_entry(pos, struct lensy_ray_struct, raylist);

		i = lensy_intersect_plane(pray, &cube0, w0, w1);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}
		if (draw)
			line(pray->p, w0, pray->red, pray->green, pray->blue);
		d0 = in_air / lensy_index_sellmeier(pray->wavelength, &N_BK7);
		i = lensy_redirect_refract (pray, w0, w1, d0);
	}

	//----------------- cube1
	list_for_each_safe(pos, pos0, &raylist) {
		pray = list_entry(pos, struct lensy_ray_struct, raylist);

		i = lensy_intersect_plane(pray, &cube1, w0, w1);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}
		if (draw)
			line(pray->p, w0, pray->red, pray->green, pray->blue);
		d0 = lensy_index_sellmeier(pray->wavelength, &N_BK7) / in_air;
		i = lensy_redirect_refract(pray, w0, w1, d0);
	}

	//---------------- focal plane
	list_for_each_safe(pos, pos0, &raylist) {
		pray = list_entry(pos, struct lensy_ray_struct, raylist);

		i = lensy_intersect_plane(pray, &ccd1.p, w0, w1);
		if (draw)
			line(pray->p, w0, pray->red, pray->green, pray->blue);
		if (i < 0) {
			list_del(pos);
			free(pray);
			continue;
		}
		lensy_redirect_impact(pray, w0, w1);
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
		free (pray);
	}

	if (make_ll_picture) {
		//------------- write a FITS file showing the focal plane
		fp = fopen("lensy.fits", "w");
		if (fp == NULL) {
			fprintf(stderr, "fopen FITS file failed\n");
			exit(0);
		}

		i = 0;
		memset (hdr, ' ', sizeof (hdr));
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


	//--- draw sections of the surfaces of the optics in the graphic window
	for (d0 = -1.5; d0 < +1.5; d0 += optic_scale/2) {
		ray.p[0] =  1.0;	// ray starting position vector in meters
		ray.p[1] =  d0;
		ray.p[2] =  0.0;

		ray.d[0] = -1.0;	// ray direction vector
		ray.d[1] =  0.0;
		ray.d[2] =  0.0;

		i = lensy_intersect_paraboloid(&ray, &primary, w0, w2);
		ray.p[1] = d0 + optic_scale;
		i = lensy_intersect_paraboloid(&ray, &primary, w1, w2);
		j = (i == 0) ? 255 : 100;
		line(w0, w1, j, j, j);

		ray.p[1] = d0;
		i = lensy_intersect_plane(&ray, &ccd1.p, w0, w2);	// &collimator2
		ray.p[1] =  d0 + optic_scale;
		i = lensy_intersect_plane(&ray, &ccd1.p, w1, w2);	// &collimator2
		j = (i == 0) ? 255 : 100;
		line(w0, w1, j, j, j);

		ray.d[0] = +1.0;	// ray direction vector
		ray.d[1] =  0.0;
		ray.d[2] =  0.0;

		ray.p[1] = d0;
		i = lensy_intersect_hyperboloid (&ray, &secondary, w0, w2);
		ray.p[1] = d0 + optic_scale;
		i = lensy_intersect_hyperboloid (&ray, &secondary, w1, w2);
		j = (i == 0) ? 255 : 100;
		line(w0, w1, j, j, j);
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
	SDL_Delay(1);
	SDL_PumpEvents();

	if (i0++ < 25) {
		secondary.v[0] += 0.015e-3;
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
		SDL_Delay(10);
	}

	//-------------------------- free the line_list
	list_for_each_safe(pos, pos0, &line_list) {
		ple = list_entry(pos, struct line_struct, line_list);
		list_del (pos);
		free (ple);
	}

	exit(0);
}
