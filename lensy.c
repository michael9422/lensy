/*
 * lensy.c - v1.4 (codename FlamingMarshmallow)
 *
 * Library functions for modelling geometric light ray paths with lenses
 * and mirrors that have simple surfaces.
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
 * These functions work by first specifying what surfaces the light rays will
 * pass through. They can be flat, spherical, parabolic, cylindrical, or
 * hyperbolic. For each surface, a structure of that type is declared and the
 * surface parameters are specified (i.e. the vertex position vector, the
 * aperture diameter, and the vector from the vertex to the center of
 * curvature, or to the focus for a parabola or hyperbola).
 *
 * Then, a light ray structure is initialized, and it is traced through the
 * optics by calculating where it intersects the surfaces, in order, and
 * adjusting the the ray direction at each surface. This is done by calling
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
 * Send questions or contributions to me at <michael.h.williamson@gmail.com>.
 *
 * History:
 *
 *  2015-06-09  Reorganize into library functions
 *  2013-11-14  Added a lensy_ccd_struct type for (flat) CCD type detector.
 *  2011-03-13  Added functions lensy_cone() and lensy_beam().
 *  2010-10-28  Added a line list for 3-D rendering in the graphic window.
 *
 * To do:
 *
 *  - Add a command line interface.
 *  - Store the parameters in a data files (i.e., the surface parameters,
 *    the material index of refraction, diffraction grating parameters, etc.)
 *    This would require significant program changes.
 *  - Check malloc return value everywhere (xmalloc?).
 *  - Make the display window zoom-in, zoom-out, and 3-D. (Perhaps use openGL
 *    instead of SDL.)
 *  - Handle non-circular, off-center apertures.
 *  - Create a function for drawing the optic elements cross-section in the
 *    graphic window.
 *  - The calculations could be done faster. Also, ray tracing could use
 *    multiple threads for increased speed.
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


//--------------------------------------------------- lensy_inner3
double lensy_inner3(double a[3], double b[3])
{
	return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}


//--------------------------------------------------- lensy_mag3
double lensy_mag3(double a[3])
{
	return (sqrt(lensy_inner3(a, a)));
}


//--------------------------------------------------- lensy_cross3
void lensy_cross3(double a[3], double b[3], double r[3])
{
	r[0] =   a[1]*b[2] - b[1]*a[2];
	r[1] = -(a[0]*b[2] - b[0]*a[2]);
	r[2] =   a[0]*b[1] - b[0]*a[1];
}


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
					double q[3], double n[3])
{
	double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10;
	double w0[3], w1[3], w2[3];

	q[0] = r->p[0];
	q[1] = r->p[1];
	q[2] = r->p[2];

	//--------- w0 (unit vector parallel to p->f)
	d0 = lensy_mag3(p->f);
	w0[0] = p->f[0] / d0;
	w0[1] = p->f[1] / d0;
	w0[2] = p->f[2] / d0;

	//--------- w1
	w1[0] = r->p[0] - p->v[0] - p->f[0];
	w1[1] = r->p[1] - p->v[1] - p->f[1];
	w1[2] = r->p[2] - p->v[2] - p->f[2];

	//------ solve a*x^2 + b*x + c = 0 for x
	d1 = lensy_inner3(r->d, r->d) - pow(lensy_inner3(r->d, w0), 2.0); // a
	d4 = 2 * lensy_mag3(p->f) + lensy_inner3(w1, w0);
	d2 = 2 * lensy_inner3(r->d, w1) - 2 * lensy_inner3(r->d, w0) * d4; // b
	d3 = lensy_inner3(w1, w1) - d4 * d4;			// c
	if (d1 == 0.0) {
		d5 = -d3 / d2;
//	} else if (fabs(d1) < 1e-6) {		// is this needed?
//		return -2;
	} else {
		d10 = (d2 * d2) - (4 * d1 * d3);
		if (d10 < 0.0) return -2;
		d5 = (-d2 + sqrt(d10)) / (2 * d1);
		d9 = (-d2 - sqrt(d10)) / (2 * d1);
		if ((d5 < 0.0) || ((d9 > 0.0) && (d9 < d5))) d5 = d9;
		if (d5 < 0.0) return -2;
	}
	q[0] = r->p[0] + d5 * r->d[0];
	q[1] = r->p[1] + d5 * r->d[1];
	q[2] = r->p[2] + d5 * r->d[2];

	w2[0] = q[0] - p->v[0];
	w2[1] = q[1] - p->v[1];
	w2[2] = q[2] - p->v[2];

	d6 = lensy_inner3(w2, w0);
	w2[0] += -d6 * w0[0];
	w2[1] += -d6 * w0[1];
	w2[2] += -d6 * w0[2];

	d7 = lensy_mag3(w2);
	if ((d7 == 0.0) /* || (fabs (d7) < 1e-6) */ ) {
		n[0] = w0[0];
		n[1] = w0[1];
		n[2] = w0[2];
	} else {
		w2[0] /= d7;
		w2[1] /= d7;
		w2[2] /= d7;

		n[0] = -d7/(2*d0) * w2[0] + w0[0];
		n[1] = -d7/(2*d0) * w2[1] + w0[1];
		n[2] = -d7/(2*d0) * w2[2] + w0[2];

		d8 = lensy_mag3(n);
		n[0] /= d8;
		n[1] /= d8;
		n[2] /= d8;
	}

	if (d7 > p->aperture / 2.0) return -1;
	return 0;
}


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
 *
 */
int32_t lensy_intersect_sphere(struct lensy_ray_struct *r,
						struct lensy_sphere_struct *s,
						double q[3], double n[3])
{
	double d0, d1, d2, d3, d4, d5, d7, d8, d9;
	double w0[3], w1[3], w2[3], w3[3], q1[3], w4[3];


	q[0] = r->p[0];
	q[1] = r->p[1];
	q[2] = r->p[2];

	//--------- w1 = s->v + s->vr (sphere center vector)
	w1[0] = s->v[0] + s->vr[0];
	w1[1] = s->v[1] + s->vr[1];
	w1[2] = s->v[2] + s->vr[2];

	//--------- w0 = r->p - s->c
	w0[0] = r->p[0] - w1[0];
	w0[1] = r->p[1] - w1[1];
	w0[2] = r->p[2] - w1[2];

	//------ solve a*x^2 + b*x + c = 0 for x
	d1 = lensy_inner3(r->d, r->d);			// a
	d2 = 2 * lensy_inner3(r->d, w0);			// b
	d3 = lensy_inner3(w0, w0) - lensy_inner3(s->vr, s->vr);	// c
	if (d1 == 0.0) {
		d4 = -d3 / d2;
//	} else if (fabs(d1) < 1e-6) {	// is this needed?
//		return -2;
	} else {
		/*
		 * Select the correct intersect point, from the two possible
		 * intersect points, by choosing the one on the 'vertex' side
		 * of the center of the sphere.
		 */
		d5 = (d2 * d2) - (4 * d1 * d3);
		if (d5 < 0.0) return -2;

		d4 = (-d2 + sqrt(d5)) / (2 * d1);
		q1[0] = r->p[0] + d4 * r->d[0];
		q1[1] = r->p[1] + d4 * r->d[1];
		q1[2] = r->p[2] + d4 * r->d[2];

		w4[0] = q1[0] - w1[0];
		w4[1] = q1[1] - w1[1];
		w4[2] = q1[2] - w1[2];

		if (lensy_inner3(w4, s->vr) >= 0.0)
			d4 = (-d2 - sqrt(d5)) / (2 * d1);
	}
	q1[0] = r->p[0] + d4 * r->d[0];
	q1[1] = r->p[1] + d4 * r->d[1];
	q1[2] = r->p[2] + d4 * r->d[2];

	w4[0] = q1[0] - w1[0];
	w4[1] = q1[1] - w1[1];
	w4[2] = q1[2] - w1[2];

	if (lensy_inner3(w4, s->vr) >= 0.0) return -2;

	q[0] = q1[0];
	q[1] = q1[1];
	q[2] = q1[2];

	n[0] = q[0] - w1[0];
	n[1] = q[1] - w1[1];
	n[2] = q[2] - w1[2];

	d0 = lensy_mag3(n);
	if (d0 > 0.0) {
		n[0] /= d0;
		n[1] /= d0;
		n[2] /= d0;
	}

	w2[0] = q[0] - s->v[0];
	w2[1] = q[1] - s->v[1];
	w2[2] = q[2] - s->v[2];

	d7 = lensy_mag3(s->vr);
	w3[0] = s->vr[0] / d7;
	w3[1] = s->vr[1] / d7;
	w3[2] = s->vr[2] / d7;

	d8 = lensy_inner3(w2, w3);
	w2[0] += -d8 * w3[0];
	w2[1] += -d8 * w3[1];
	w2[2] += -d8 * w3[2];

	d9 = lensy_mag3(w2);
	if (d9 > s->aperture / 2.0) return -1;
	return 0;
}


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
					double q[3], double n[3])
{
	double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, da, db, dc;
	double a[3], w0[3], w1[3], w2[3], w3[3], q1[3], w4[3], w5[3], w6[3];

	q[0] = r->p[0];
	q[1] = r->p[1];
	q[2] = r->p[2];

	/*------- Make the vector 'a' parallel to the cylinder axis and ensure
	 * it is perpendicular to c->va and has unit length.
	 */
	d0 = lensy_mag3(c->va);
	if (d0 == 0.0) return -2;
	w0[0] = c->va[0] / d0;
	w0[1] = c->va[1] / d0;
	w0[2] = c->va[2] / d0;

	d1 = lensy_inner3(c->a, w0);
	a[0] =  c->a[0] - d1 * w0[0];
	a[1] =  c->a[1] - d1 * w0[1];
	a[2] =  c->a[2] - d1 * w0[2];
	d2 = lensy_mag3(a);
	if (d2 == 0.0) return -2;
	a[0] /= d2;
	a[1] /= d2;
	a[2] /= d2;

	//--------- w1 = c->v + c->va (cylinder center vector)
	w1[0] = c->v[0] + c->va[0];
	w1[1] = c->v[1] + c->va[1];
	w1[2] = c->v[2] + c->va[2];

	//--------- w2 = r->p - w1
	w2[0] = r->p[0] - w1[0];
	w2[1] = r->p[1] - w1[1];
	w2[2] = r->p[2] - w1[2];

	//------ solve a*x^2 + b*x + c = 0 for x
	d3 = lensy_inner3(w2, a);
	d4 = lensy_inner3(r->d, a);
	w3[0] = r->d[0] - d4 * a[0];
	w3[1] = r->d[1] - d4 * a[1];
	w3[2] = r->d[2] - d4 * a[2];

	w4[0] = w2[0] - d3 * a[0];
	w4[1] = w2[1] - d3 * a[1];
	w4[2] = w2[2] - d3 * a[2];

	da = lensy_inner3(w3, w3);
	db = 2 * lensy_inner3(w3, w4);
	dc = lensy_inner3(w4, w4) - lensy_inner3(c->va, c->va);

	if (da == 0.0) {
		d5 = -dc / db;
//	} else if (fabs(da) < 1e-6) {	// is this needed?
//		return -2;
	} else {
		/*
		 * Select the correct intersect point, from the two possible
		 * intersect points, by choosing the one on the 'vertex' side
		 * of the center of the cylinder.
		 */
		d6 = (db * db) - (4 * da * dc);
		if (d6 < 0.0) return -2;

		d5 = (-db + sqrt(d6)) / (2 * da);
		q1[0] = r->p[0] + d5 * r->d[0];
		q1[1] = r->p[1] + d5 * r->d[1];
		q1[2] = r->p[2] + d5 * r->d[2];

		w5[0] = q1[0] - w1[0];
		w5[1] = q1[1] - w1[1];
		w5[2] = q1[2] - w1[2];

		if (lensy_inner3(w5, c->va) >= 0.0)
			d5 = (-db - sqrt(d6)) / (2 * da);
	}
	q1[0] = r->p[0] + d5 * r->d[0];
	q1[1] = r->p[1] + d5 * r->d[1];
	q1[2] = r->p[2] + d5 * r->d[2];

	w5[0] = q1[0] - w1[0];
	w5[1] = q1[1] - w1[1];
	w5[2] = q1[2] - w1[2];

	if (lensy_inner3(w5, c->va) >= 0.0) return -2;

	q[0] = q1[0];
	q[1] = q1[1];
	q[2] = q1[2];

	//--------------------- compute the surface normal
	d7 = lensy_inner3(w5, a);
	n[0] = w5[0] - d7 * a[0];
	n[1] = w5[1] - d7 * a[1];
	n[2] = w5[2] - d7 * a[2];

	d8 = lensy_mag3(n);
	if (d8 == 0.0) return -2;
	n[0] /= d8;
	n[1] /= d8;
	n[2] /= d8;

	//------- determine if the intersect point is inside the aperture
	d9 = lensy_inner3(w5, w0);
	w6[0] = w5[0] - d9 * w0[0];
	w6[1] = w5[1] - d9 * w0[1];
	w6[2] = w5[2] - d9 * w0[2];

	d10 = lensy_mag3(w6);
	if (d10 > c->aperture / 2.0) return -1;
	return 0;
}


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
						double q[3], double n[3])
{
	double d0, d1, d2, d3;
	double w[3];

	q[0] = r->p[0];
	q[1] = r->p[1];
	q[2] = r->p[2];

	d0 = lensy_inner3(r->d, p->n);
	if ((d0 == 0.0) /* || (fabs (d0) < 1e-6) */) return -2;

	d1 = (lensy_inner3(p->v, p->n) - lensy_inner3(r->p, p->n)) / d0;
	if (d1 < 0.0) return -2;

	q[0] = r->p[0] + d1 * r->d[0];
	q[1] = r->p[1] + d1 * r->d[1];
	q[2] = r->p[2] + d1 * r->d[2];

	d3 = lensy_mag3(p->n);
	if (d3 == 0.0) return -2;
	n[0] = p->n[0] / d3;
	n[1] = p->n[1] / d3;
	n[2] = p->n[2] / d3;

	w[0] = q[0] - p->v[0];
	w[1] = q[1] - p->v[1];
	w[2] = q[2] - p->v[2];
	d2 = lensy_mag3(w);

	if (d2 > p->aperture / 2.0) return -1;

	return 0;
}


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
					double q[3], double n[3])
{
	double d0, d1, d2, d3, d4, d5, d7, d8, d9;
	double f[3], w0[3], w1[3], w2[3], w3[3], w4[3], q1[3];


	q[0] = r->p[0];
	q[1] = r->p[1];
	q[2] = r->p[2];

	//--------- w1 = h->v + h->a (hyperboloid center)
	w1[0] = h->v[0] + h->a[0];
	w1[1] = h->v[1] + h->a[1];
	w1[2] = h->v[2] + h->a[2];

	//--------- focus
	f[0] = w1[0] + (h->e * -h->a[0]);
	f[1] = w1[1] + (h->e * -h->a[1]);
	f[2] = w1[2] + (h->e * -h->a[2]);

	//---------
	d8 = lensy_mag3(h->a);
	w2[0] = -h->a[0] / d8;
	w2[1] = -h->a[1] / d8;
	w2[2] = -h->a[2] / d8;

	w3[0] = r->p[0] - w1[0] + h->a[0] / h->e;
	w3[1] = r->p[1] - w1[1] + h->a[1] / h->e;
	w3[2] = r->p[2] - w1[2] + h->a[2] / h->e;

	w4[0] = r->p[0] - f[0];
	w4[1] = r->p[1] - f[1];
	w4[2] = r->p[2] - f[2];

	//------ solve a*x^2 + b*x + c = 0 for x
	d0 = h->e * h->e;
	d4 = lensy_inner3(w2, r->d);
	d5 = lensy_inner3(w2, w3);

	d1 = lensy_inner3(r->d, r->d) - (d0 * d4 * d4);		// a
	d2 = 2 * (lensy_inner3(r->d, w4) - (d0 * d4 * d5));	// b
	d3 = lensy_inner3(w4, w4) - (d0 * d5 * d5);		// c

	if (d1 == 0.0) {
		d4 = -d3 / d2;
//	} else if (fabs (d1) < 1e-6) {	// is this needed?
//		return -2;
	} else {
		/*
		 * Select the correct intersect point, from the two possible
		 * intersect points.
		 */
		d5 = (d2 * d2) - (4 * d1 * d3);
		if (d5 < 0.0) return -2;

		d4 = (-d2 + sqrt(d5)) / (2 * d1);
		q1[0] = r->p[0] + d4 * r->d[0];
		q1[1] = r->p[1] + d4 * r->d[1];
		q1[2] = r->p[2] + d4 * r->d[2];

		w4[0] = q1[0] - w1[0];
		w4[1] = q1[1] - w1[1];
		w4[2] = q1[2] - w1[2];

		if (lensy_inner3(w4, h->a) >= 0.0)
			d4 = (-d2 - sqrt(d5))/(2 * d1);
	}
	q1[0] = r->p[0] + d4 * r->d[0];
	q1[1] = r->p[1] + d4 * r->d[1];
	q1[2] = r->p[2] + d4 * r->d[2];

	w4[0] = q1[0] - w1[0];
	w4[1] = q1[1] - w1[1];
	w4[2] = q1[2] - w1[2];

	if (lensy_inner3(w4, h->a) >= 0.0) return -2;

	q[0] = q1[0];
	q[1] = q1[1];
	q[2] = q1[2];

	w0[0] = q[0] - h->v[0];
	w0[1] = q[1] - h->v[1];
	w0[2] = q[2] - h->v[2];

	d0 = lensy_inner3(w0, w2);
	w0[0] += -d0 * w2[0];
	w0[1] += -d0 * w2[1];
	w0[2] += -d0 * w2[2];

	d9 = lensy_mag3(w0);
	if (d0 == 0.0) {
		n[0] = w2[0];
		n[1] = w2[1];
		n[2] = w2[2];
	} else {
		w0[0] /= d9;
		w0[1] /= d9;
		w0[2] /= d9;

		d0 = sqrt((d8 * d8) * (h->e * h->e - 1));
		d1 = (d8 / d0) * (d9 / sqrt(d0 * d0 + d9 * d9));

		n[0] = w2[0] - d1 * w0[0];
		n[1] = w2[1] - d1 * w0[1];
		n[2] = w2[2] - d1 * w0[2];
	}

	d7 = lensy_mag3(n);
	if (d7 > 0.0) {
		n[0] /= d7;
		n[1] /= d7;
		n[2] /= d7;
	} else {
		return -2;
	}

	if (d9 > h->aperture / 2.0) return -1;
	return 0;
}


/*----------------------------------------------------- lensy_redirect_reflect
 * 'r' is the ray to be reflected.
 * 'q' is the 3-D intersect point.
 * 'n' is the unit normal vector to the surface.
 */
void lensy_redirect_reflect(struct lensy_ray_struct *r,
						double q[3], double n[3])
{
	double d0;

	r->p[0] = q[0];
	r->p[1] = q[1];
	r->p[2]	= q[2];

	d0 = lensy_inner3(r->d, n);
	r->d[0] += -2 * d0 * n[0];
	r->d[1] += -2 * d0 * n[1];
	r->d[2] += -2 * d0 * n[2];
}



/*----------------------------------------------------- lensy_redirect_refract
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
					double q[3], double n[3], double m)
{
	double d0, d1, d2, d3;
	double n1[3], u[3], w[3], w1[3], v[3];

	r->p[0]	=  q[0];
	r->p[1]	=  q[1];
	r->p[2]	=  q[2];

	d0 = lensy_mag3(r->d);
	if (d0 == 0.0) return -2;

	u[0] = -r->d[0] / d0;
	u[1] = -r->d[1] / d0;
	u[2] = -r->d[2] / d0;

	n1[0] = n[0];
	n1[1] = n[1];
	n1[2] = n[2];
//	if (fabs(lensy_inner3(n1, n1) - 1.0) > 1e-6) return -2;

	if (lensy_inner3(u, n1) < 0.0) {
		n1[0] = -n1[0];
		n1[1] = -n1[1];
		n1[2] = -n1[2];
	}
	lensy_cross3(u, n1, w);
	d1 = lensy_mag3(w);
	d2 = m * d1;
	if (fabs(d2) >= 1.0) return -1;
	d3 = asin(d2);	// angle of transmission w.r.t. the surface normal

	if (d1 > 0.0) {
		w1[0] = w[0] / d1;
		w1[1] = w[1] / d1;
		w1[2] = w[2] / d1;

		lensy_cross3(w1, n1, v);
		r->d[0] = d0 * (cos(d3) * (-n1[0]) + sin(d3) * v[0]);
		r->d[1] = d0 * (cos(d3) * (-n1[1]) + sin(d3) * v[1]);
		r->d[2] = d0 * (cos(d3) * (-n1[2]) + sin(d3) * v[2]);
	} else {
		r->d[0] = d0 * (-n1[0]);
		r->d[1] = d0 * (-n1[1]);
		r->d[2] = d0 * (-n1[2]);
	}
	return 0;
}


/*--------------------------------------------------- lensy_redirect_diffract
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
				double wli, double wlt, int32_t m)
{
	double wli1, wlt1;
	double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10;
	double n1[3], a1[3], t1[3], w0[3], w1[3];


	r->p[0]	=  q[0];
	r->p[1]	=  q[1];
	r->p[2]	=  q[2];

	d3 = lensy_mag3(n);
	if (d3 == 0.0) return -2;
	n1[0] = n[0] / d3;
	n1[1] = n[1] / d3;
	n1[2] = n[2] / d3;

	d0 = lensy_mag3(r->d);
	if (d0 == 0.0) return -2;
	w0[0] = r->d[0] / d0;
	w0[1] = r->d[1] / d0;
	w0[2] = r->d[2] / d0;

	d1 = lensy_mag3(a);			// the ruling spacing
	d2 = lensy_inner3(a, n1);
	a1[0] = a[0] - d2 * n1[0];	// ensure a1 is perpendicular to n
	a1[1] = a[1] - d2 * n1[1];
	a1[2] = a[2] - d2 * n1[2];

	d2 = lensy_mag3(a1);
	if (d2 == 0.0) return -2;
	a1[0] /= d2;
	a1[1] /= d2;
	a1[2] /= d2;

	lensy_cross3(a1, n1, t1);

	d4 = lensy_inner3(w0, n1);	// (-) for transmit, (+) for reflect
	d5 = lensy_inner3(w0, a1);
	d6 = lensy_inner3(w0, t1);

	if (d4 == 0.0) return -2;
	if (d6 == 1.0) return -2;
	d10 = 1.0 / sqrt(1.0 - d6 * d6);
	wli1 = wli * d10;
	wlt1 = wlt * d10;
	d7 = atan2(d5, -d4);
	d8 = (sin(d7) / wli1 + m / d1) * wlt1;
	if (fabs(d8) >= 1.0) return -2;
	d9 = asin(d8);

	w1[0] = d6 * t1[0] + (cos(d9) / d10) * n1[0] + (sin(d9) / d10) * a1[0];
	w1[1] = d6 * t1[1] + (cos(d9) / d10) * n1[1] + (sin(d9) / d10) * a1[1];
	w1[2] = d6 * t1[2] + (cos(d9) / d10) * n1[2] + (sin(d9) / d10) * a1[2];

	r->d[0] = d0 * w1[0];
	r->d[1] = d0 * w1[1];
	r->d[2] = d0 * w1[2];

	return (0);
}



/*----------------------------------------------------- lensy_redirect_impact
 * The ray reaches the intersect point (focal plane).
 *
 * 'r' is the ray to be impacted.
 * 'q' is the 3-D intersect point.
 * 'n' is the unit normal vector to the surface.
 */
void lensy_redirect_impact(struct lensy_ray_struct *r,
						double q[3], double n[3]) {
	r->p[0]	= q[0];
	r->p[1]	= q[1];
	r->p[2]	= q[2];
}



//---------- coefficients for calculating the index of refraction
double CaF2[6] = { 2.0388472e0, -3.2320997e-3, 6.1568960e-3,
                   5.6612714e-5, -4.0951444e-9, 2.2406560e-8 };

double tsu2[6] = { 2.5310795e0, -1.0750804e-2, 1.4091541e-2,
                   2.4479041e-4, -4.3396907e-6, 4.2269287e-7 };

double tsu4[6] = { 2.5310397e0, -1.0751078e-2, 1.4089396e-2,
                   2.4455705e-4, -4.3189009e-6, 4.2184152e-7 };

double tsu5[6] = { 2.2182723e0, -5.2937745e-3, 8.4751835e-3,
                   9.0035648e-5, -2.1638749e-7, 8.8532657e-8 };

double tsu6[6] = { 2.3863743e0, -9.2750923e-3, 1.2963764e-2,
                   2.6012532e-4, -7.1806739e-6, 6.4902518e-7 };

double tsu7[6] = { 2.5309288e0, -1.0751176e-2, 1.4087125e-2,
                   2.4433615e-4, -4.2994607e-6, 4.2104219e-7 };

double fsilica[6] = { 2.1045254e0, 9.5251763e-3, 8.5795589e-3,
                   1.2770234e-4, -2.2841020e-6, 1.2397250e-7 };


/*------------------------------------------------- lensy_index_of_refraction
 * Return the index of refraction for the wavelength wl in meters.
 */
double lensy_index_of_refraction(double wl, double a[6])
{
	if ((wl < 0.3e-6) || (wl > 2.0e-6)) {
		fprintf(stderr, "%s: wavelength outside limits\n", __func__);
		exit(-1);
	}

	wl *= 1.0e6;
	return (sqrt(a[0] + a[1] * pow(wl, 2) + a[2] / pow(wl, 2) +
		a[3] / pow(wl, 4) + a[4] / pow(wl, 6) + a[5] / pow(wl, 8)));
}


//----- values for calculating index of refraction using sellmeier formula
struct lensy_sellmeier_struct N_BAF10 = {
   .b1 = 1.58514950e+00,
   .b2 = 1.43559385e-01,
   .b3 = 1.08521269e+00,
   .c1 = 9.26681282e-03,
   .c2 = 4.24489805e-02,
   .c3 = 1.05613573e+02
};

struct lensy_sellmeier_struct N_SF6 = {
   .b1 = 1.77931763e+00,
   .b2 = 3.38149866e-01,
   .b3 = 2.08734474e+00,
   .c1 = 1.33714182e-02,
   .c2 = 6.17533621e-02,
   .c3 = 1.74017590e+02
};

struct lensy_sellmeier_struct N_BK7 = {
   .b1 = 1.03961212E+00,
   .b2 = 2.31792344E-01,
   .b3 = 1.01046945E+00,
   .c1 = 6.00069867E-03,
   .c2 = 2.00179144E-02,
   .c3 = 1.03560653E+02
};

struct lensy_sellmeier_struct SF2 = {
   .b1 = 1.40301821E+00,
   .b2 = 2.31767504E-01,
   .b3 = 9.39056586E-01,
   .c1 = 1.05795466E-02,
   .c2 = 4.93226978E-02,
   .c3 = 1.12405955E+02
};


/*------------------------------------------------------ lensy_index_sellmeier
 * Return the index of refraction using the Sellmeier formula.
 * Wavelength wl is in meters.
 */
double lensy_index_sellmeier(double wl, struct lensy_sellmeier_struct *a)
{
	double d0, d1;

	wl *= 1.0e6;
	d0 = wl * wl;
	d1 = (a->b1 * d0) / (d0 - a->c1) + (a->b2 * d0) / (d0 - a->c2) +
		 (a->b3 * d0) / (d0 - a->c3);
	return (sqrt(d1 + 1.0));
}


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
				double cone_dia, double cone_step)
{
	int32_t j, k, m, n, rc;
	double d0, d1, d2, d3, d4, d5, d10;
	struct lensy_ray_struct *pray;
	double w0[3], w1[3], w2[3], w3[3], u0[3], u1[3], u2[3];

	rc = 0;

	d10 = lensy_mag3(pr->d);
	if (d10 == 0.0) {
		fprintf(stderr, "%s: lensy_mag3(d) == 0.0 (ray direction is null)\n", __func__);
		return rc;
	}

	//--- Add the center ray to the list first
	pray = (struct lensy_ray_struct *) malloc(sizeof(struct lensy_ray_struct));
	if (pray != NULL) {
		pray->p[0] = pr->p[0];
		pray->p[1] = pr->p[1];
		pray->p[2] = pr->p[2];
		pray->d[0] = pr->d[0];
		pray->d[1] = pr->d[1];
		pray->d[2] = pr->d[2];
		pray->wavelength = pr->wavelength;
		snprintf(pray->pathkey, sizeof(pray->pathkey), "%e%e%e%e",
			pray->p[0], pray->p[1], pray->p[2], pray->wavelength);
		list_add(&(pray->raylist), pl);
		rc++;
	}

	//------------- Now add rays inside a cone to the raylist
	d0 = DEG2RAD * cone_step;
	n = floor((DEG2RAD * cone_dia / 2) / d0);
	for (j = 1; j <= n; j++) {
		d1 = j * d0;
		m = floor(sin(d1) * 2 * PI / d0);
		for (k = 0; k < m; k++) {
			d2 = k * (2 * PI / m);
			d3 = PI / 2 - d1;

			//------ w0, w1, w2 are in 'unprimed' coordinates
			w0[0] = 0.0;
			w0[1] = 0.0;
			w0[2] = 1.0;

			w1[0] = cos(d3) * cos(d2);
			w1[1] = cos(d3) * sin(d2);
			w1[2] = sin(d3);

			w2[0] = w1[0] - w0[0];
			w2[1] = w1[1] - w0[1];
			w2[2] = w1[2] - w0[2];

			d4 = atan2(pr->d[1], pr->d[0]);
			d5 = asin(pr->d[2]);

			//------- unit vectors in the 'primed' coordinates
			u0[0] = cos(PI / 2 - d4);
			u0[1] = cos(d4) * cos(PI / 2 - d5);
			u0[2] = cos(d4) * sin(PI / 2 - d5);

			u1[0] = sin(-(PI / 2 - d4));
			u1[1] = cos(-(PI / 2 - d4)) * cos(PI / 2 - d5);
			u1[2] = cos(-(PI / 2 - d4)) * sin(PI / 2 - d5);

			u2[0] = 0.0;
			u2[1] = sin(-(PI / 2 - d5));
			u2[2] = cos(-(PI / 2 - d5));

			w3[0] = pr->d[0] + d10 * lensy_inner3(w2, u0);
			w3[1] = pr->d[1] + d10 * lensy_inner3(w2, u1);
			w3[2] = pr->d[2] + d10 * lensy_inner3(w2, u2);

			//--- Add the ray to the list
			pray = (struct lensy_ray_struct *) malloc(sizeof(struct lensy_ray_struct));
			if (pray != NULL) {
				pray->p[0] = pr->p[0];
				pray->p[1] = pr->p[1];
				pray->p[2] = pr->p[2];

				pray->d[0] = w3[0];
				pray->d[1] = w3[1];
				pray->d[2] = w3[2];

				pray->wavelength = pr->wavelength;

				pray->red = pr->red;
				pray->green = pr->green;
				pray->blue = pr->blue;

				snprintf(pray->pathkey, sizeof(pray->pathkey), "%e%e%e%e",
					pray->p[0], pray->p[1], pray->p[2], pray->wavelength);
				list_add(&(pray->raylist), pl);
				rc++;
			}
		}
	}

	return rc;
}



/*----------------------------------------------------- lensy_beam
 * Create a circular beam of parallel rays in a list 'pl', centered around the
 * ray 'pr'.
 *
 * Input parameters:
 *
 *	beam_dia	- The beam diameter in meters.
 *	beam_step	- The spacing between rays in the beam, in meters.
 *
 */
int32_t lensy_beam(struct list_head *pl, struct lensy_ray_struct *pr,
				double beam_dia, double beam_step)
{
	int32_t rc;
	double d0, d1, d2, d10;
	struct lensy_ray_struct *pray;
	double w0[3], u0[3], u1[3];

	rc = 0;

	d10 = lensy_mag3(pr->d);
	if (d10 == 0.0) {
		fprintf(stderr, "%s: lensy_mag3(d) == 0.0 (ray direction is null)\n", __func__);
		return rc;
	}

	w0[0] = pr->d[0] / d10;
	w0[1] = pr->d[1] / d10;
	w0[2] = pr->d[2] / d10;

	u0[0] = w0[1] / sqrt(w0[0] * w0[0] + w0[1] * w0[1]);
	u0[1] = sqrt(1 - u0[0] * u0[0]);
	u0[2] = 0.0;

	lensy_cross3(w0, u0, u1);

	for (d0 = -beam_dia / 2; d0 < +beam_dia / 2; d0 += beam_step) {
		for (d1 = -beam_dia / 2; d1 < +beam_dia / 2; d1 += beam_step) {
			d2 = sqrt(d0 * d0 + d1 * d1);
			if (d2 > beam_dia / 2) continue;

			//--- Add a ray to the list
			pray = (struct lensy_ray_struct *) malloc(sizeof(struct lensy_ray_struct));
			if (pray != NULL) {
				pray->p[0] = pr->p[0] + d0 * u0[0] + d1 * u1[0];
				pray->p[1] = pr->p[1] + d0 * u0[1] + d1 * u1[1];
				pray->p[2] = pr->p[2] + d0 * u0[2] + d1 * u1[2];

				pray->d[0] = pr->d[0];
				pray->d[1] = pr->d[1];
				pray->d[2] = pr->d[2];

				pray->wavelength = pr->wavelength;

				pray->red = pr->red;
				pray->green = pr->green;
				pray->blue = pr->blue;

				snprintf(pray->pathkey, sizeof(pray->pathkey), "%e%e%e%e",
					pray->d[0], pray->d[1], pray->d[2], pray->wavelength);
				list_add(&(pray->raylist), pl);
				rc++;
			}
		}
	}

	return rc;
}


/*------------------------------------------ lensy_init_ccd
 * This function fills in the values for the associated plane structure,
 * and allocates space for an image buffer.
 */
void lensy_init_ccd(struct lensy_ccd_struct *ccd) {
	double d0, w0[3];

	ccd->p.v[0] = ccd->v[0];
	ccd->p.v[1] = ccd->v[1];
	ccd->p.v[2] = ccd->v[2];

	lensy_cross3(ccd->vx, ccd->vy, w0);

	d0 = lensy_mag3(w0);
	if (d0 == 0.0) {
		fprintf(stderr, "%s: invalid ccd parameters vx, vy\n", __func__);
		exit(-1);
	}

	ccd->p.n[0] = w0[0] / d0;
	ccd->p.n[1] = w0[1] / d0;
	ccd->p.n[2] = w0[2] / d0;

	ccd->b_size = sizeof(*ccd->b) * ccd->x_nmax * ccd->y_nmax;
	ccd->b = (uint16_t *) malloc(ccd->b_size);
	if (ccd->b == NULL) {
		fprintf(stderr, "%s: malloc ccd buffer failed\n", __func__);
		exit(-1);
	}

	ccd->p.aperture = 2 * (ccd->x_nmax * lensy_mag3(ccd->vx) +
				 ccd->y_nmax * lensy_mag3(ccd->vy));
}
