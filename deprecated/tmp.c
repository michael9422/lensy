/*
 * Compile with:
 *
 *	# gcc tmp.c -lSDL2
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

//--------- SDL-2 header files
#include <SDL2/SDL.h>
/*
#include <SDL/SDL_hints.h>
#include <SDL/SDL_error.h>
#include <SDL/SDL_log.h>
#include <SDL/SDL_assert.h>
#include <SDL/SDL_version.h>
#include <SDL/SDL_video.h>
#include <SDL/SDL_render.h>
#include <SDL/SDL_pixels.h>
#include <SDL/SDL_rect.h>
#include <SDL/SDL_surface.h>
#include <SDL/SDL_syswm.h>
#include <SDL/SDL_clipboard.h>
#include <SDL/SDL_events.h>
#include <SDL/SDL_keyboard.h>
#include <SDL/SDL_keycode.h>
#include <SDL/SDL_scancode.h>
#include <SDL/SDL_mouse.h>
#include <SDL/SDL_joystick.h>
#include <SDL/SDL_gamecontroller.h>
*/
/*
#include <SDL/SDL_video.h>
#include <SDL/begin_code.h>
#include <SDL/close_code.h>
*/

SDL_Window *optic_win;
SDL_Window *focus_win;
SDL_Renderer *optic_ren;


//---------------------------------------------------- main
int main(int argc, char *argv[]) {
	int n;

	n = SDL_Init(SDL_INIT_EVERYTHING);
	if (n != 0) {
		fprintf(stderr, "initialize SDL failed: %s\n", SDL_GetError());
		return -1;
	}
	atexit(SDL_Quit);

	optic_win = SDL_CreateWindow("Optics",
					SDL_WINDOWPOS_UNDEFINED,
					SDL_WINDOWPOS_UNDEFINED,
					640, 480, 0);

	optic_ren = SDL_CreateRenderer(optic_win, -1, SDL_RENDERER_ACCELERATED);
	if ((optic_win == NULL) || (optic_ren == NULL)) {
		fprintf(stderr, "create SDL window failed\n");
		exit(-1);
	}
/*
	focus_win = SDL_CreateWindow("Focal Plane",
					SDL_WINDOWPOS_UNDEFINED,
					SDL_WINDOWPOS_UNDEFINED,
					640, 480, 0);
*/
	SDL_SetHint(SDL_HINT_RENDER_SCALE_QUALITY, "linear");
	SDL_RenderSetLogicalSize(optic_ren, 640, 480);
	SDL_SetRenderDrawColor(optic_ren, 0, 0, 0, 255);
	SDL_RenderClear(optic_ren);
	SDL_RenderPresent(optic_ren);

	SDL_SetRenderDrawColor(optic_ren, 255, 0, 0, 255);
	SDL_RenderDrawPoint(optic_ren, 20, 240);
	SDL_SetRenderDrawColor(optic_ren, 0, 255, 0, 255);
	SDL_RenderDrawLine(optic_ren, -300, -300, 11300, 11300);
	SDL_RenderPresent(optic_ren);

	sleep(3);

	exit(0);
}
