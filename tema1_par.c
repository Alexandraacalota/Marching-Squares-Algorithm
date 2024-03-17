// Author: APD team, except where source was noted

#include "helpers.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>

#define CONTOUR_CONFIG_COUNT    16
#define FILENAME_MAX_SIZE       50
#define STEP                    8
#define SIGMA                   200
#define RESCALE_X               2048
#define RESCALE_Y               2048

#define CLAMP(v, min, max) if(v < min) { v = min; } else if(v > max) { v = max; }
#define min(a, b) ((a) < (b) ? (a) : (b))

// Definying a type with all the necessary information for each thread's
// function
typedef struct {
	ppm_image *image;
	ppm_image *new_image;  // rescaled image
	ppm_image **contour_map;
	unsigned char **grid;
	int step_x, step_y;
	int id;  // id of the thread
	int num_threads;  // total number of threads
	pthread_barrier_t *barrier;
} func_arg;

// Creates a map between the binary configuration (e.g. 0110_2) and the
// corresponding pixels that need to be set on the output image. An array is
// used for this map since the keys are binary numbers in 0-15. Contour images
// are located in the './contours' directory.
ppm_image **init_contour_map() {
	ppm_image **map = (ppm_image **)malloc(CONTOUR_CONFIG_COUNT *
										   sizeof(ppm_image *));
	if (!map) {
		fprintf(stderr, "Unable to allocate memory\n");
		exit(1);
	}

	for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
		char filename[FILENAME_MAX_SIZE];
		sprintf(filename, "./contours/%d.ppm", i);
		map[i] = read_ppm(filename);
	}

	return map;
}

// Corresponds to step 1 of the marching squares algorithm, which focuses on
// sampling the image.
// Builds a p x q grid of points with values which can be either 0 or 1,
// depending on how the pixel values compare to the `sigma` reference value.
// The points are taken at equal distances in the original image, based on the
// `step_x` and `step_y` arguments.
unsigned char **sample_grid(ppm_image *image, int step_x, int step_y,
							unsigned char sigma, unsigned char **grid,
							int id, int num_threads)
{
	int p = image->x / step_x;
	int q = image->y / step_y;

	// calculate the lines of the grid that each thread works on
	// based on the id of the thread
	int start = id * (double)p / num_threads;
	int end = min((id + 1) * (double)p / num_threads, p);

	// each thread works on its own interval [start, end) of the grid
	for (int i = start; i < end; i++) {
		for (int j = 0; j < q; j++) {
			ppm_pixel curr_pixel = image->data[i * step_x * image->y +
											   j * step_y];

			unsigned char curr_color = (curr_pixel.red + curr_pixel.green +
										curr_pixel.blue) / 3;

			if (curr_color > sigma) {
				grid[i][j] = 0;
			} else {
				grid[i][j] = 1;
			}
		}
	}
	grid[p][q] = 0;

	// last sample points have no neighbors below / to the right, so we use
	// pixels on the last row / column of the input image for them
	
	
	for (int i = start; i < end; i++) {
		ppm_pixel curr_pixel = image->data[i * step_x * image->y +
										   image->x - 1];

		unsigned char curr_color = (curr_pixel.red + curr_pixel.green +
									curr_pixel.blue) /
								   3;

		if (curr_color > sigma) {
			grid[i][q] = 0;
		} else {
			grid[i][q] = 1;
		}
	}

	// recalculate end for q values instead of p
	end = min((id + 1) * (double)q / num_threads, q);
	for (int j = start; j < end; j++) {
		ppm_pixel curr_pixel = image->data[(image->x - 1) * image->y +
										   j * step_y];

		unsigned char curr_color = (curr_pixel.red + curr_pixel.green +
									curr_pixel.blue) / 3;

		if (curr_color > sigma) {
			grid[p][j] = 0;
		} else {
			grid[p][j] = 1;
		}
	}

	return grid;
}

// Updates a particular section of an image with the corresponding contour pixels.
// Used to create the complete contour image.
void update_image(ppm_image *image, ppm_image *contour, int x, int y) {
	for (int i = 0; i < contour->x; i++) {
		for (int j = 0; j < contour->y; j++) {
			int contour_pixel_index = contour->x * i + j;
			int image_pixel_index = (x + i) * image->y + y + j;

			image->data[image_pixel_index].red = contour->
												 data[contour_pixel_index].red;
			image->data[image_pixel_index].green = contour->
												   data[contour_pixel_index].green;
			image->data[image_pixel_index].blue = contour->
												  data[contour_pixel_index].blue;
		}
	}
}

// Corresponds to step 2 of the marching squares algorithm, which focuses on
// identifying the type of contour which corresponds to each subgrid. It
// determines the binary value of each sample fragment of the original image
// and replaces the pixels in the original image with the pixels of the
// corresponding contour image accordingly.
void march(ppm_image *image, unsigned char **grid, ppm_image **contour_map,
		   int step_x, int step_y, int id, int num_threads)
{
	int p = image->x / step_x;
	int q = image->y / step_y;

	int start = id * (double)p / num_threads;
	int end = min((id + 1) * (double)p / num_threads, p);

	// the lines on which each thread operates
	for (int i = start; i < end; i++) {
		for (int j = 0; j < q; j++) {
			unsigned char k = 8 * grid[i][j] + 4 * grid[i][j + 1] + 2 *
							  grid[i + 1][j + 1] + 1 * grid[i + 1][j];
			update_image(image, contour_map[k], i * step_x, j * step_y);
		}
	}
}

// Calls `free` method on the utilized resources.
void free_resources(ppm_image *image, ppm_image **contour_map,
					unsigned char **grid, int step_x)
{
	for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
		free(contour_map[i]->data);
		free(contour_map[i]);
	}
	free(contour_map);

	for (int i = 0; i <= image->x / step_x; i++) {
		free(grid[i]);
	}
	free(grid);

	free(image->data);
	free(image);
}

// Called to allocate memory for a new_image in case the initial one is
// too large and it needs rescaling
ppm_image *allocate_new_image() {
	// alloc memory for image
	ppm_image *new_image = (ppm_image *)malloc(sizeof(ppm_image));
	if (!new_image) {
		fprintf(stderr, "Unable to allocate memory\n");
		exit(1);
	}
	new_image->x = RESCALE_X;
	new_image->y = RESCALE_Y;

	new_image->data = (ppm_pixel *)malloc(new_image->x * new_image->y *
										  sizeof(ppm_pixel));
	if (!new_image) {
		fprintf(stderr, "Unable to allocate memory\n");
		exit(1);
	}
	return new_image;
}

// Function called to allocate memory for the grid before the parallel
// part of the program begins
unsigned char **allocate_grid(ppm_image *image, int step_x, int step_y) {
	int p = min(image->x / step_x, RESCALE_X);
	int q = min(image->y / step_y, RESCALE_Y);

	unsigned char **grid = (unsigned char **)malloc((p + 1) *
													sizeof(unsigned char *));
	if (!grid) {
		fprintf(stderr, "Unable to allocate memory\n");
		exit(1);
	}

	for (int i = 0; i <= p; i++) {
		grid[i] = (unsigned char *)malloc((q + 1) * sizeof(unsigned char));
		if (!grid[i]) {
			fprintf(stderr, "Unable to allocate memory\n");
			exit(1);
		}
	}
	return grid;
}

// Function called if the initial image's dimensions are beyond RESCALE_X
// and RESCALE_Y
ppm_image *rescale_image(ppm_image *image, ppm_image *new_image, int thread_id,
						 int num_threads)
{
	uint8_t sample[3];

	int start = thread_id * (double)new_image->x / num_threads;
	int end = min((thread_id + 1) * (double)new_image->x / num_threads,
				  new_image->x);

	// use bicubic interpolation for scaling
	// each thread works on its own lines of the image in [start, end)
	for (int i = start; i < end; i++) {
		for (int j = 0; j < new_image->y; j++) {
			float u = (float)i / (float)(new_image->x - 1);
			float v = (float)j / (float)(new_image->y - 1);
			sample_bicubic(image, u, v, sample);

			new_image->data[i * new_image->y + j].red = sample[0];
			new_image->data[i * new_image->y + j].green = sample[1];
			new_image->data[i * new_image->y + j].blue = sample[2];
		}
	}

	return new_image;
}

// Allocates memory for the grid and, if necessary, for the new_image that
// has the rescaled dimensions. Creates the threads and sets the argument
// for the function of each thread. Also initializez a barrier and distroys
// it at the end of the parallel part of the program.
void apply(ppm_image *image, ppm_image **contour_map, int num_threads,
		   int step_x, int step_y, void *(*f)(void *), func_arg *argument,
		   int counter)
{
	pthread_t threads[num_threads];
	int r;
	int id;
	void *status;
	pthread_barrier_t barrier;

	pthread_barrier_init(&barrier, NULL, num_threads);

	unsigned char **grid = allocate_grid(image, step_x, step_y);

	ppm_image *new_image;

	if (counter == 2) {
		new_image = allocate_new_image();
	}
	
	for (id = 0; id < num_threads; id++) {
		argument[id].id = id;
		argument[id].step_x = step_x;
		argument[id].step_y = step_y;
		argument[id].num_threads = num_threads;
		argument[id].image = image;
		argument[id].barrier = &barrier;
		if (counter == 2) {
			argument[id].new_image = new_image;
		}
		argument[id].grid = grid;
		argument[id].contour_map = contour_map;
		r = pthread_create(&threads[id], NULL, f, &argument[id]);

		if (r) {
			printf("Error while creating the thread %d\n", id);
			exit(-1);
		}
	}

	for (id = 0; id < num_threads; id++) {
		r = pthread_join(threads[id], &status);

		if (r) {
			printf("Error while waiting the thread %d\n", id);
			exit(-1);
		}
	}

	pthread_barrier_destroy(&barrier);
}

// Threads function for small images
void *f1(void *arg) {
	func_arg argument = *(func_arg *)arg;

	// Sample the grid
	unsigned char **grid = sample_grid(argument.image, argument.step_x,
									   argument.step_y, SIGMA, argument.grid,
									   argument.id, argument.num_threads);

	pthread_barrier_wait(argument.barrier);

	// March the squares
	march(argument.image, grid, argument.contour_map, argument.step_x,
		  argument.step_y, argument.id, argument.num_threads);

	pthread_exit(NULL);
}

// Threads function for larger sized images that need rescaling
void *f2(void *arg) {
	func_arg argument = *(func_arg *)arg;

	// Rescale the image
	argument.new_image = rescale_image(argument.image, argument.new_image,
									   argument.id, argument.num_threads);

	pthread_barrier_wait(argument.barrier);

	// Sample the grid
	unsigned char **grid = sample_grid(argument.new_image, argument.step_x,
									   argument.step_y, SIGMA, argument.grid,
									   argument.id, argument.num_threads);

	pthread_barrier_wait(argument.barrier);

	// March the squares
	march(argument.new_image, grid, argument.contour_map, argument.step_x,
		  argument.step_y, argument.id, argument.num_threads);
	pthread_exit(NULL);
}

int main(int argc, char *argv[]) {
	if (argc < 4) {
		fprintf(stderr, "Usage: ./tema1 <in_file> <out_file> <P>\n");
		return 1;
	}

	int num_threads = atoi(argv[3]);
	ppm_image *image = read_ppm(argv[1]);
	int step_x = STEP;
	int step_y = STEP;

	// Initialize contour map
	ppm_image **contour_map = init_contour_map();

	func_arg argument[num_threads];

	if (image->x <= RESCALE_X && image->y <= RESCALE_Y) {
		// apply functions without rescaling the initial image
		apply(image, contour_map, num_threads, step_x, step_y, f1, argument, 1);
		// Write output
		write_ppm(argument[0].image, argv[2]);
	} else {
		// apply functions including rescaling the initial image
		apply(image, contour_map, num_threads, step_x, step_y, f2, argument, 2);
		// Write output
		write_ppm(argument[0].new_image, argv[2]);
		free(argument[0].new_image->data);
		free(argument[0].new_image);
	}

	free_resources(argument[0].image, contour_map, argument[0].grid, step_x);

	return 0;
}
