/*
 * lensy.h - v1.4 (codename FlamingMarshmallow)
 *
 * This is the header file for the lensy optics library functions.
 *
 *-----------------------------------------------------------------------
 * Copyright (C) 2015 FSF
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
 * The calculations are done for three dimensional rays and surfaces.
 *
 * Send questions or contributions to me at <michael.h.williamson@gmail.com>.
 *
 *-----------------------------------------------------------------------
 * History:
 *
 *   2015-06-09 This file is created.
 *
 */

#ifndef LENSY_H
#define LENSY_H		1

#include <ctype.h>
#include <float.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <search.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/errno.h>
#include <sys/types.h>
#include <sys/un.h>
#include <unistd.h>

#include <list.h>

#ifndef PI
 #define PI		3.14159265
#endif

#define DEG2RAD		0.01745329252
#define RAD2DEG		57.29577951



/*----------------------------- surface structures
 * Arrays have elements of
 *
 *	 <[0], [1], [2]> = <x, y, z>.
 *
 * All values are in meters.
 */

struct lensy_paraboloid_struct {
	double v[3];		// vertex position
	double f[3];		// vector from vertex to focus
	double aperture;	// circular aperture diameter
};

struct lensy_sphere_struct {
	double v[3];		// vertex position
	double vr[3];		// vector from vertex to center of the sphere
	double aperture;	// circular aperture diameter
};

struct lensy_plane_struct {
	double v[3];		// vertex position
	double n[3];		// normal vector to the plane
	double aperture;	// circular aperture diameter
};

struct lensy_ccd_struct {
	double v[3];		// vertex position
	double vx[3], vy[3];	// pixel axis vectors for a CCD -
				// 	vector lengths determine pixel size -
				// 	and must be perpendicular vectors
	int x_nmax, y_nmax;	// detector dimensions
	uint16_t *b;		// image buffer
	int b_size;		// buffer size in bytes
	struct lensy_plane_struct p;	// plane structure dervied from above values
};

struct lensy_cylinder_struct {
	double v[3];		// vertex position
	double va[3];		// vector from vertex to the cylinder axis
	double a[3];		// vector parallel to the cylinder axis
	double aperture;	// circular aperture diameter
};

struct lensy_hyperboloid_struct {
	double v[3];		// vertex position
	double a[3];		// vector from vertex to the center
	double e;		// eccentricity (e > 1).
	double aperture;	// circular aperture diameter
};


//------------------------------ light ray structure
struct lensy_ray_struct {
	double p[3];		// 3-D vector position of the light ray
	double d[3];		// ray direction vector
	double wavelength;	// wavelength in vacuum (meters)
	struct list_head raylist;
	char red, green, blue;
	char pathkey[80];	// ray path history, for spot size calculation
};


//--------------------------------------------------- vector functions
double lensy_inner3(double a[3], double b[3]);

double lensy_mag3(double a[3]);

void lensy_cross3(double a[3], double b[3], double r[3]);


/*--------------------------------------------- lensy_intersect_paraboloid
 * Calculate where a ray intersects a paraboloid.
 *
 * The ray structure 'r' and the paraboloid structure 'p' passed to this
 * function must be initialized with values such that the ray will intersect
 * the paraboloid, starting from the start position 'p', by multiplying it
 * (the ray) by a positive value.
 *
 * The return vector q[3] is the calculated intersect point.
 * The return vector n[3] is the unit normal vector to the surface at the
 * intersect point.
 *
 * A return value of zero means OK.
 * A return value of -1 means that the ray is outside the aperture.
 * A return value of -2 means that there is no intersect point.
 *
 */
int32_t lensy_intersect_paraboloid(struct lensy_ray_struct *r,
					   struct lensy_paraboloid_struct *p,
					   double q[3], double n[3]);


/*------------------------------------------------ lensy_intersect_sphere
 * Calculate where a ray intersects a sphere.
 *
 * The ray structure 'r' and the sphere structure 's' passed to this function
 * must be initialized with values such that the ray will intersect the sphere,
 * starting from the start position 'p', by multiplying it (the ray) by a
 * positive value.
 *
 * The return vector q[3] is the calculated intersect point.
 * The return vector n[3] is the unit (outward) normal vector to the surface
 * at the intersect point.
 *
 * A return value of zero means OK.
 * A return value of -1 means that the ray is outside the aperture.
 * A return value of -2 means that there is no intersect point.
 *
 */
int32_t lensy_intersect_sphere(struct lensy_ray_struct *r,
						struct lensy_sphere_struct *s,
						double q[3], double n[3]);


/*---------------------------------------------- lensy_intersect_cylinder
 * Calculate where a ray intersects a cylinder.
 *
 * The ray structure 'r' and the cylinder structure 'c' passed to this function
 * must be initialized with values such that the ray will intersect the
 * cylinder, starting from the start position 'p', by multiplying it (the ray)
 * by a positive value.
 *
 * The return vector q[3] is the calculated intersect point.
 * The return vector n[3] is the unit (outward) normal vector to the surface
 * at the intersect point.
 *
 * A return value of zero means OK.
 * A return value of -1 means that the ray is outside the aperture.
 * A return value of -2 means that there is no intersect point.
 *
 */
int32_t lensy_intersect_cylinder(struct lensy_ray_struct *r,
					struct lensy_cylinder_struct *c,
					double q[3], double n[3]);


/*---------------------------------------------- lensy_intersect_plane
 * Calculate where a ray intersects a plane.
 *
 * The ray structure 'r' and the plane structure 'p' passed to this function
 * must be initialized with values such that the ray will intersect the plane,
 * starting from the start position 'p', by multiplying it (the ray) by a
 * positive value.
 *
 * The return vector q[3] is the calculated intersect point.
 * The return vector n[3] is the unit normal vector to the surface.
 *
 * A return value of zero means OK.
 * A return value of -1 means that the ray is outside the aperture.
 * A return value of -2 means that there is no intersect point.
 *
 */
int32_t lensy_intersect_plane(struct lensy_ray_struct *r,
						struct lensy_plane_struct *p,
						double q[3], double n[3]);


/*------------------------------------------- lensy_intersect_hyperboloid
 * Calculate where a ray intersects a hyperboloid.
 *
 * The ray structure 'r' and the hyperboloid structure 'h' passed to this
 * function must be initialized with values such that the ray will intersect
 * the hyperboloid, starting from the start position 'p', by multiplying it
 * (the ray) by a positive value.
 *
 * The return vector q[3] is the calculated intersect point.
 * The return vector n[3] is the unit (outward) normal vector to the surface
 * at the intersect point.
 *
 * A return value of zero means OK.
 * A return value of -1 means that the ray is outside the aperture.
 * A return value of -2 means that there is no intersect point.
 *
 */
int32_t lensy_intersect_hyperboloid(struct lensy_ray_struct *r,
					struct lensy_hyperboloid_struct *h,
					double q[3], double n[3]);


/*-------------------------------------------------- lensy_redirect_reflect
 * 'r' is the ray to be reflected.
 * 'q' is the 3-D intersect point.
 * 'n' is the unit normal vector to the surface.
 */
void lensy_redirect_reflect(struct lensy_ray_struct *r,
						 double q[3], double n[3]);


/*--------------------------------------------------- lensy_redirect_refract
 * 'r'  is the ray to be refracted.
 * 'q'  is the 3-D intersect point.
 * 'n'  is the unit (outward) normal vector to the surface.
 * 'm'  is the ratio of the index of refraction for the incident medium divided
 *      by the index of refraction for the transmission medium.
 *
 * A return value of zero means OK.
 * A return value of -1 means 'total internal reflection'.
 */
int32_t lensy_redirect_refract(struct lensy_ray_struct *r,
					double q[3], double n[3], double m);


/*-------------------------------------------------- lensy_redirect_diffract
 * Diffract 'ray' from a grating. If the ray direction inner product with
 * the surface normal is negative, 'ray' is reflected. Otherwise, it is
 * transmitted.
 *
 * 'r'  is the ray to be diffracted.
 * 'q'  is the 3-D intersect point.
 * 'n'  is the unit normal vector to the surface.
 * 'a'  is the 3-D vector perpendicular to the diffraction grating rulings
 *      and the surface normal vector, whose length is equal to the spacing
 *      between adjacent rulings.
 * 'wli' is the wavelength of the incident light in meters.
 * 'wlt' is the wavelength of the reflected or transmitted light (it should
 *       equal 'wli' for reflection).
 * 'm'  is the order that should be used (..., -2, -1, 0, +1, +2, ...) for
 *      reflecting.
 *
 * A return value of zero means OK.
 * A return value of -1 means 'invalid'.
 * A return value of -2 means 'invalid'.
 */
int32_t lensy_redirect_diffract(struct lensy_ray_struct *r,
				double q[3], double n[3], double a[3],
				double wli, double wlt, int32_t m);


/*---------------------------------------------------- lensy_redirect_impact
 * The ray reaches the intersect point (focal plane).
 *
 * 'r' is the ray to be impacted.
 * 'q' is the 3-D intersect point.
 * 'n' is the unit normal vector to the surface.
 */
void lensy_redirect_impact(struct lensy_ray_struct *r,
						double q[3], double n[3]);



//---------- coefficients for calculating the index of refraction
extern double CaF2[6];
extern double tsu2[6];
extern double tsu4[6];
extern double tsu5[6];
extern double tsu6[6];
extern double tsu7[6];
extern double fsilica[6];


/*-------------------------------------------- lensy_index_of_refraction
 * Return the index of refraction for the wavelength wl in meters.
 */
double lensy_index_of_refraction(double wl, double a[6]);



//----- values for calculating index of refraction using sellmeier formula
struct lensy_sellmeier_struct {
   double b1, b2, b3, c1, c2, c3;
};


extern struct lensy_sellmeier_struct N_BAF10;
extern struct lensy_sellmeier_struct N_SF6;
extern struct lensy_sellmeier_struct N_BK7;
extern struct lensy_sellmeier_struct SF2;


/*---------------------------------------------------- lensy_index_sellmeier
 * Return the index of refraction using the Sellmeier formula.
 * Wavelength wl is in meters.
 */
double lensy_index_sellmeier(double wl, struct lensy_sellmeier_struct *a);



/*----------------------------------------------------- lensy_cone
 * Create a cone of rays in a list 'pl', centered around the ray 'pr'.
 *
 * Input parameters:
 *
 *	cone_dia	- The diameter of a cone of rays, in degrees.
 *	cone_step	- The angular spacing of rays in the cone, in degrees.
 *
 *
 */
int32_t lensy_cone(struct list_head *pl, struct lensy_ray_struct *pr,
				double cone_dia, double cone_step);


/*----------------------------------------------------- lensy_beam
 * Create a circular beam of parallel rays in a list 'pl', centered around
 * the ray 'pr'.
 *
 * Input parameters:
 *
 *	beam_dia	- The beam diameter in meters.
 *	beam_step	- The spacing between rays in the beam, in meters.
 *
 */
int32_t lensy_beam(struct list_head *pl, struct lensy_ray_struct *pr,
				double beam_dia, double beam_step);


/*----------------------------------------------------- lensy_init_ccd
 * This function fills in the values for the associated plane structure,
 * and allocates space for an image buffer.
 */
void lensy_init_ccd(struct lensy_ccd_struct *ccd);

#endif
