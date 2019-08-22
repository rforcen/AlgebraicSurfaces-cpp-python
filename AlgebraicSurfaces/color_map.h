#include <stdlib.h>

typedef struct {
    float r,g,b;
} Color;

Color color_map(float v, float vmin, float vmax, int type);
