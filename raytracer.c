#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Constants
#define WIDTH 800
#define HEIGHT 600
#define MAX_DEPTH 5
#define EPSILON 1e-4

// Vector structure and operations
typedef struct {
    double x, y, z;
} Vector;

// Function prototypes for vector operations
Vector vec_add(Vector a, Vector b);
Vector vec_sub(Vector a, Vector b);
Vector vec_mul(Vector a, double scalar);
Vector vec_div(Vector a, double scalar);
double vec_dot(Vector a, Vector b);
Vector vec_normalize(Vector a);
Vector vec_reflect(Vector v, Vector n);
Vector vec_hadamard(Vector a, Vector b);

// Vector operations implementation
Vector vec_add(Vector a, Vector b) {
    return (Vector){a.x + b.x, a.y + b.y, a.z + b.z};
}

Vector vec_sub(Vector a, Vector b) {
    return (Vector){a.x - b.x, a.y - b.y, a.z - b.z};
}

Vector vec_mul(Vector a, double scalar) {
    return (Vector){a.x * scalar, a.y * scalar, a.z * scalar};
}

Vector vec_div(Vector a, double scalar) {
    return (Vector){a.x / scalar, a.y / scalar, a.z / scalar};
}

double vec_dot(Vector a, Vector b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vector vec_normalize(Vector a) {
    double len = sqrt(vec_dot(a, a));
    return vec_div(a, len);
}

Vector vec_reflect(Vector v, Vector n) {
    return vec_sub(v, vec_mul(n, 2 * vec_dot(v, n)));
}

// Hadamard product (component-wise multiplication)
Vector vec_hadamard(Vector a, Vector b) {
    return (Vector){a.x * b.x, a.y * b.y, a.z * b.z};
}

// Sphere structure
typedef struct {
    Vector center;
    double radius;
    Vector color;
    double reflectivity;
} Sphere;

// Light structure
typedef struct {
    Vector position;
    Vector color;
} Light;

// Ray structure
typedef struct {
    Vector origin;
    Vector direction;
} Ray;

// Scene objects
Sphere spheres[] = {
    {{0, 0, -20}, 4, {1, 0, 0}, 0.5},  // Red sphere
    {{5, -1, -15}, 2, {0, 1, 0}, 0.5}, // Green sphere
    {{-5, 0, -25}, 3, {0, 0, 1}, 0.5}, // Blue sphere
};

Light light = {{-10, 10, 10}, {1, 1, 1}};

// Function prototypes for ray tracing
int intersect_sphere(Ray ray, Sphere sphere, double *t);
Vector trace_ray(Ray ray, int depth);
void save_ppm(const char *filename, unsigned char *image);

// Intersect a ray with a sphere
int intersect_sphere(Ray ray, Sphere sphere, double *t) {
    Vector oc = vec_sub(ray.origin, sphere.center);
    double a = vec_dot(ray.direction, ray.direction);
    double b = 2.0 * vec_dot(oc, ray.direction);
    double c = vec_dot(oc, oc) - sphere.radius * sphere.radius;
    double discriminant = b * b - 4 * a * c;

    if (discriminant < 0) return 0;
    double sqrt_disc = sqrt(discriminant);

    double t0 = (-b - sqrt_disc) / (2.0 * a);
    double t1 = (-b + sqrt_disc) / (2.0 * a);
    *t = (t0 < t1 && t0 > EPSILON) ? t0 : t1;

    return *t > EPSILON;
}

// Trace a ray and compute its color
Vector trace_ray(Ray ray, int depth) {
    double t_min = INFINITY;
    Sphere *hit_sphere = NULL;

    // Find the nearest intersection
    for (int i = 0; i < sizeof(spheres) / sizeof(Sphere); i++) {
        double t;
        if (intersect_sphere(ray, spheres[i], &t) && t < t_min) {
            t_min = t;
            hit_sphere = &spheres[i];
        }
    }

    if (!hit_sphere) return (Vector){0, 0, 0}; // Background color

    Vector hit_point = vec_add(ray.origin, vec_mul(ray.direction, t_min));
    Vector normal = vec_normalize(vec_sub(hit_point, hit_sphere->center));
    Vector light_dir = vec_normalize(vec_sub(light.position, hit_point));

    // Ambient lighting
    Vector ambient = vec_mul(light.color, 0.1);

    // Diffuse lighting
    double diffuse_intensity = fmax(0, vec_dot(normal, light_dir));
    Vector diffuse = vec_mul(vec_hadamard(hit_sphere->color, light.color), diffuse_intensity);

    // Specular lighting
    Vector view_dir = vec_mul(ray.direction, -1);
    Vector reflect_dir = vec_reflect(vec_mul(light_dir, -1), normal);
    double specular_intensity = pow(fmax(0, vec_dot(view_dir, reflect_dir)), 32);
    Vector specular = vec_mul(light.color, specular_intensity);

    // Combine lighting components
    Vector color = vec_add(ambient, vec_add(diffuse, specular));

    // Reflection
    if (depth < MAX_DEPTH && hit_sphere->reflectivity > 0) {
        Vector reflect_dir = vec_reflect(ray.direction, normal);
        Ray reflect_ray = {vec_add(hit_point, vec_mul(normal, EPSILON)), reflect_dir};
        Vector reflected_color = trace_ray(reflect_ray, depth + 1);
        color = vec_add(color, vec_mul(reflected_color, hit_sphere->reflectivity));
    }

    return color;
}

// Save the image as a PPM file
void save_ppm(const char *filename, unsigned char *image) {
    FILE *file = fopen(filename, "wb");
    if (!file) {
        perror("Failed to save image");
        exit(1);
    }

    fprintf(file, "P6\n%d %d\n255\n", WIDTH, HEIGHT);
    fwrite(image, 1, WIDTH * HEIGHT * 3, file);
    fclose(file);
}

// Main function
int main() {
    unsigned char *image = malloc(WIDTH * HEIGHT * 3);
    if (!image) {
        perror("Failed to allocate memory for image");
        return 1;
    }

    Vector camera = {0, 0, 0};
    double viewport_size = 1.0;
    double aspect_ratio = (double)WIDTH / HEIGHT;
    double focal_length = 1.0;

    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            double px = (2 * ((x + 0.5) / WIDTH) - 1) * aspect_ratio * viewport_size;
            double py = (1 - 2 * ((y + 0.5) / HEIGHT)) * viewport_size;
            Vector direction = vec_normalize((Vector){px, py, -focal_length});
            Ray ray = {camera, direction};

            Vector color = trace_ray(ray, 0);
            int index = (y * WIDTH + x) * 3;
            image[index] = (unsigned char)(fmin(color.x, 1) * 255);
            image[index + 1] = (unsigned char)(fmin(color.y, 1) * 255);
            image[index + 2] = (unsigned char)(fmin(color.z, 1) * 255);
        }
    }

    save_ppm("output.ppm", image);
    free(image);

    printf("Ray tracing complete! Image saved as output.ppm\n");
    return 0;
}

