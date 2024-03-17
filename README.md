Copyright 2023 Alexandra-Maria Calota
# Marching Squares Algorithm

## Overview

The project focuses on parallelizing the contour curves drawing process using the Marching Squares Algorithm. Initially implemented sequentially, the algorithm is optimized through parallelization using pthreads to accelerate key steps in the program.

## Implementation Details

- **Memory Management:** Memory allocations are moved outside the thread function to facilitate shared memory access. Memory allocated for the initial image is freed once at the end of the program, in the sequential section.
- **Thread Creation:** Threads are created using `pthread_create()` in the `apply` function, with parameters including the array of threads, function pointers, and thread-specific arguments.
- **Function Refactoring:** Two new functions, `allocate_new_image()` and `allocate_grid()`, are introduced for memory allocation purposes.
- **Mathematical Formulas:** Uniform distribution of work among threads is ensured using mathematical formulas to calculate the start and end indices for each thread's processing.
- **Thread Synchronization:** `pthread_join()` is employed to synchronize the main thread with the termination of all created threads. Additionally, a barrier is implemented using `pthread_barrier_init()` to synchronize threads between specific function calls.
- **Barrier Usage:** The barrier is used to ensure synchronization between the `rescale_image()`, `sample_grid()`, and `march()` functions, allowing threads to progress only after completing the preceding step.
- **Barrier Destruction:** At the end of the `apply()` function, the barrier is destroyed using `pthread_barrier_destroy()`.

## Conclusion

The parallel contour curves drawing project demonstrates the effective use of pthreads to enhance the performance of the Marching Squares Algorithm. Through meticulous thread management and synchronization, the algorithm's execution is optimized, resulting in significant speedup and improved scalability.