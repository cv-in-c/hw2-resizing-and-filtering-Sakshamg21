#include <math.h>
#include "image.h"
image make_image(int w, int h, int c)
{
    image im;
    im.w = w;
    im.h = h;
    im.c = c;
    im.data = calloc(w * h * c, sizeof(float));
    return im;
}
float get_pixel(image im, int x, int y, int c)
{
    // Handle out-of-bounds indices by clamping to valid range
    x = fmin(fmax(x, 0), im.w - 1);
    y = fmin(fmax(y, 0), im.h - 1);
    return im.data[x + y * im.w + c * im.w * im.h];
}
void set_pixel(image im, int x, int y, int c, float value)
{
    // Handle out-of-bounds indices by ignoring the operation
    if (x < 0 || x >= im.w || y < 0 || y >= im.h)
        return;
    im.data[x + y * im.w + c * im.w * im.h] = value;
}
float nn_interpolate(image im, float x, float y, int c)
{
    // TODO Fill in
    int nearest_x = (int)round(x);
    int nearest_y = (int)round(y);

    // Get the pixel value at the rounded coordinates
    return get_pixel(im, nearest_x, nearest_y, c);
}

image nn_resize(image im, int w, int h)
{
    image resized = make_image(new_w, new_h, im.c);

    float w_scale = (float)im.w / new_w;
    float h_scale = (float)im.h / new_h;

    for (int c = 0; c < im.c; c++) {
        for (int y = 0; y < new_h; y++) {
            for (int x = 0; x < new_w; x++) {
                int nearest_x = (int)round(x * w_scale);
                int nearest_y = (int)round(y * h_scale);
                float value = get_pixel(im, nearest_x, nearest_y, c);
                set_pixel(resized, x, y, c, value);
            }
        }
    }

    return resized;
}

float bilinear_interpolate(image im, float x, float y, int c)
{
    // TODO
    int x1 = (int)floor(x);
    int x2 = x1 + 1;
    int y1 = (int)floor(y);
    int y2 = y1 + 1;

    // Get the pixel values at the four corners
    float q11 = get_pixel(im, x1, y1, c);
    float q12 = get_pixel(im, x1, y2, c);
    float q21 = get_pixel(im, x2, y1, c);
    float q22 = get_pixel(im, x2, y2, c);

    // Calculate the fractional part of the coordinates
    float dx = x - x1;
    float dy = y - y1;

    // Perform bilinear interpolation
    float value = (1 - dx) * (1 - dy) * q11 +
                  dx * (1 - dy) * q21 +
                  (1 - dx) * dy * q12 +
                  dx * dy * q22;
    return value;
}

image bilinear_resize(image im, int w, int h)
{
    // TODO
    image resized = make_image(new_w, new_h, im.c);

    float w_scale = (float)im.w / new_w;
    float h_scale = (float)im.h / new_h;

    for (int c = 0; c < im.c; c++) {
        for (int y = 0; y < new_h; y++) {
            for (int x = 0; x < new_w; x++) {
                float src_x = (x + 0.5) * w_scale - 0.5;
                float src_y = (y + 0.5) * h_scale - 0.5;

                int x1 = (int)floor(src_x);
                int x2 = x1 + 1;
                int y1 = (int)floor(src_y);
                int y2 = y1 + 1;

                float q11 = get_pixel(im, x1, y1, c);
                float q12 = get_pixel(im, x1, y2, c);
                float q21 = get_pixel(im, x2, y1, c);
                float q22 = get_pixel(im, x2, y2, c);

                float dx = src_x - x1;
                float dy = src_y - y1;

                float value = (1 - dx) * (1 - dy) * q11 +
                              dx * (1 - dy) * q21 +
                              (1 - dx) * dy * q12 +
                              dx * dy * q22;

                set_pixel(resized, x, y, c, value);
            }
        }
    }

    return resized;
}

