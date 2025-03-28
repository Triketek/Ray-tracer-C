#include <stdio.h>
#include <math.h>

// Structures
typedef struct {
    double x, y, z;
} Vector;

typedef struct {
    Vector origin, direction;
} Ray;

typedef struct {
    Vector center, color;
    double radius;
} Sphere;

typedef struct {
    Vector position, normal, color;
} Plane;

typedef struct {
    Vector position, color;
} Light;

// Constants
#define MAX_DEPTH 3
#define EPSILON 1e-4

// Function declarations
Vector vec_add(Vector a, Vector b);
Vector vec_sub(Vector a, Vector b);
Vector vec_mul(Vector a, double scalar);
Vector vec_hadamard(Vector a, Vector b);
Vector vec_normalize(Vector a);
double vec_dot(Vector a, Vector b);
double vec_length(Vector a);
int intersect_sphere(Ray ray, Sphere sphere, double *t);
int intersect_plane(Ray ray, Plane plane, double *t);
int is_in_shadow(Vector point, Vector light_dir);
Vector trace_ray(Ray ray, int depth);

// Global variables
Sphere spheres[] = {
    {{0, -0.5, 3}, {1, 0, 0}, 0.5},
    {{-1, -0.5, 4}, {0, 1, 0}, 0.5},
    {{1, -0.5, 4}, {0, 0, 1}, 0.5}
};

Plane plane = {{0, -1, 0}, {0, 1, 0}, {0.3, 0.3, 0.3}}; // Floor at y = -1, facing upwards

Light light = {{2, 5, 0}, {1, 1, 1}};

// Vector operations
Vector vec_add(Vector a, Vector b) {
    return (Vector){a.x + b.x, a.y + b.y, a.z + b.z};
}

Vector vec_sub(Vector a, Vector b) {
    return (Vector){a.x - b.x, a.y - b.y, a.z - b.z};
}

Vector vec_mul(Vector a, double scalar) {
    return (Vector){a.x * scalar, a.y * scalar, a.z * scalar};
}

Vector vec_hadamard(Vector a, Vector b) {
    return (Vector){a.x * b.x, a.y * b.y, a.z * b.z};
}

Vector vec_normalize(Vector a) {
    double length = sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
    return vec_mul(a, 1.0 / length);
}

double vec_dot(Vector a, Vector b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

double vec_length(Vector a) {
    return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

// Sphere intersection
int intersect_sphere(Ray ray, Sphere sphere, double *t) {
    Vector oc = vec_sub(ray.origin, sphere.center);
    double a = vec_dot(ray.direction, ray.direction);
    double b = 2.0 * vec_dot(oc, ray.direction);
    double c = vec_dot(oc, oc) - sphere.radius * sphere.radius;
    double discriminant = b * b - 4 * a * c;

    if (discriminant < 0) return 0;
    double sqrt_d = sqrt(discriminant);

    double t1 = (-b - sqrt_d) / (2 * a);
    double t2 = (-b + sqrt_d) / (2 * a);

    if (t1 > EPSILON) *t = t1;
    else if (t2 > EPSILON) *t = t2;
    else return 0;

    return 1;
}

// Plane intersection
int intersect_plane(Ray ray, Plane plane, double *t) {
    double denom = vec_dot(ray.direction, plane.normal);
    if (fabs(denom) < EPSILON) return 0;

    *t = vec_dot(vec_sub(plane.position, ray.origin), plane.normal) / denom;
    return *t > EPSILON;
}

// Shadow detection
int is_in_shadow(Vector point, Vector light_dir) {
    Ray shadow_ray = {point, light_dir};
    double t;

    // Check sphere shadows
    for (int i = 0; i < sizeof(spheres) / sizeof(Sphere); i++) {
        if (intersect_sphere(shadow_ray, spheres[i], &t) && t > EPSILON) {
            return 1;
        }
    }

    // Check plane shadows
    if (intersect_plane(shadow_ray, plane, &t) && t > EPSILON) {
        return 1;
    }

    return 0;
}

// Ray tracing
Vector trace_ray(Ray ray, int depth) {
    Sphere *hit_sphere = NULL;
    Plane *hit_plane = NULL;
    double t_min = 1e6, t;
    Vector color = {0, 0, 0};

    // Check sphere intersections
    for (int i = 0; i < sizeof(spheres) / sizeof(Sphere); i++) {
        if (intersect_sphere(ray, spheres[i], &t) && t < t_min) {
            t_min = t;
            hit_sphere = &spheres[i];
        }
    }

    // Check plane intersection
    if (intersect_plane(ray, plane, &t) && t < t_min) {
        t_min = t;
        hit_plane = &plane;
        hit_sphere = NULL;
    }

    if (hit_sphere) {
        Vector hit_point = vec_add(ray.origin, vec_mul(ray.direction, t_min));
        Vector normal = vec_normalize(vec_sub(hit_point, hit_sphere->center));

        Vector light_dir = vec_normalize(vec_sub(light.position, hit_point));
        Vector ambient = vec_mul(hit_sphere->color, 0.1);

        Vector adjusted_point = vec_add(hit_point, vec_mul(normal, EPSILON));  // Push point slightly along normal
        int in_shadow = is_in_shadow(adjusted_point, light_dir);
        if (in_shadow) return ambient;

        double diffuse_intensity = fmax(0, vec_dot(normal, light_dir));
        Vector diffuse = vec_mul(vec_hadamard(hit_sphere->color, light.color), diffuse_intensity);

        color = vec_add(ambient, diffuse);
    } else if (hit_plane) {
        Vector hit_point = vec_add(ray.origin, vec_mul(ray.direction, t_min));
        Vector normal = plane.normal;

        Vector light_dir = vec_normalize(vec_sub(light.position, hit_point));
        Vector ambient = vec_mul(plane.color, 0.1);

        Vector adjusted_point = vec_add(hit_point, vec_mul(plane.normal, EPSILON));
        int in_shadow = is_in_shadow(adjusted_point, light_dir);
        if (in_shadow) return ambient;

        double diffuse_intensity = fmax(0, vec_dot(normal, light_dir));
        Vector diffuse = vec_mul(vec_hadamard(plane.color, light.color), diffuse_intensity);

        color = vec_add(ambient, diffuse);
    }

    return color;
}

// Render scene
void render_scene(int width, int height) {
    FILE *file = fopen("output.ppm", "w");
    fprintf(file, "P3\n%d %d\n255\n", width, height);

    double fov = M_PI / 3.0;

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            double px = (2 * (x + 0.5) / (double)width - 1) * tan(fov / 2) * width / (double)height;
            double py = (1 - 2 * (y + 0.5) / (double)height) * tan(fov / 2);

            Ray ray = {{0, 0, 0}, vec_normalize((Vector){px, py, 1})};
            Vector color = trace_ray(ray, 0);

            fprintf(file, "%d %d %d ", (int)(fmin(color.x, 1) * 255), (int)(fmin(color.y, 1) * 255), (int)(fmin(color.z, 1) * 255));
        }
    }

    fclose(file);
}

int main() {
    int width = 800, height = 600;
    render_scene(width, height);
    return 0;
}

