#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

void l1_normalize(image im)
{
    // TODO
    float sum = 0.0;
 
    // Compute the sum of all pixel values
    for (int i = 0; i < im.width * im.height * im.channels; i++) {
        sum += im.data[i];
    }

    // Normalize each pixel by dividing by the sum
    for (int i = 0; i < im.width * im.height * im.channels; i++) {
        im.data[i] /= sum;
    }
}

image make_image(int width, int height, int channels)
{
    image im;
    im.width = width;
    im.height = height;
    im.channels = channels;
    im.data = (float*)malloc(width * height * channels * sizeof(float));
    return im;
}

image make_box_filter(int w)
{
    int size = w * w;
    image filter = make_image(w, w, 1);

    // Compute the value for each element in the filter
    float value = 1.0 / size;
    for (int i = 0; i < size; i++) {
        filter.data[i] = value;
    }

    return filter;
}

image convolve_image(image im, image filter, int preserve)
{
    // TODO
    int outWidth = (preserve) ? im.width : im.width - filter.width + 1;
    int outHeight = (preserve) ? im.height : im.height - filter.height + 1;
    int outChannels = im.channels;

    // Create the output image
    image output = make_image(outWidth, outHeight, outChannels);

    // Perform convolution for each pixel
    for (int c = 0; c < outChannels; c++) {
        for (int h = 0; h < outHeight; h++) {
            for (int w = 0; w < outWidth; w++) {
                // Compute the convolved value for the current pixel
                float convValue = 0.0;
                for (int fc = 0; fc < filter.channels; fc++) {
                    for (int fh = 0; fh < filter.height; fh++) {
                        for (int fw = 0; fw < filter.width; fw++) {
   int imX = w + fw;
                            int imY = h + fh;
                            int imC = c + fc;
                            int filterX = fw;
                            int filterY = fh;
                            int filterC = fc;

                            // Handle edge cases for out-of-bounds indices
                            if (imX < 0 || imX >= im.width || imY < 0 || imY >= im.height || imC < 0 || imC >= im.channels) {
                                continue;
                            }

                            // Compute the convolution value
                            convValue += im.data[(imY * im.width + imX) * im.channels + imC] *
                                         filter.data[(filterY * filter.width + filterX) * filter.channels + filterC];
                        }
                    }
                }
                 output.data[(h * outWidth + w) * outChannels + c] = convValue;
            }
        }
    }

    return output;
}

image make_highpass_filter()
{
     int filterSize = 3;
    image filter = make_image(filterSize, filterSize, 1);

    // Define the high-pass filter values
    float filterData[9] = {-1, -1, -1, -1, 8, -1, -1, -1, -1};

    // Set the filter values
    for (int i = 0; i < filterSize * filterSize; i++) {
        filter.data[i] = filterData[i];
    }

    return filter;
}

image make_sharpen_filter()
{
    int filterSize = 3;
    image filter = make_image(filterSize, filterSize, 1);

    // Define the sharpen filter values
    float filterData[9] = {0, -1, 0, -1, 5, -1, 0, -1, 0};

    // Set the filter values
    for (int i = 0; i < filterSize * filterSize; i++) {
        filter.data[i] = filterData[i];
    }

    return filter;
}

image make_emboss_filter()
{
    int filterSize = 3;
    image filter = make_image(filterSize, filterSize, 1);

    // Define the emboss filter values
    float filterData[9] = {-2, -1, 0, -1, 1, 1, 0, 1, 2};

    // Set the filter values
    for (int i = 0; i < filterSize * filterSize; i++) {
        filter.data[i] = filterData[i];
    }

    return filter;
}

// Question 2.2.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: Preserving the image size: In some cases, you may want the output image to have the same dimensions as the input image after convolution.
//This is useful when the spatial information needs to be preserved, such as in tasks like edge detection or object recognition. In these cases, the 
//preserve flag should be set to 1.

//Reducing the image size: In other cases, the filter may be larger than the input image, and the intention is to perform convolution with the filter sliding
//over the entire image. This process will result in an output image that is smaller in size compared to the input image. This is often used in tasks like 
//downsampling, blurring, or feature extraction. In these cases, the preserve flag should be set to 0.

// Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: Box Filter: The box filter is a simple averaging filter that does not require any post-processing. It calculates the average value of the neighboring pixels within the filter window.

//Highpass Filter: The highpass filter enhances high-frequency components in an image. It may produce negative values or values outside the usual image range. In such cases, post-processing can be applied to normalize the output or clip values outside the valid range.

//Sharpen Filter: The sharpen filter enhances edges and details in an image. It involves subtracting a blurred version of the image from the original. Post-processing may involve clipping values to ensure they remain within the valid range or normalizing the output.

//Emboss Filter: The emboss filter creates a 3D embossed effect in an image. It often produces output values outside the usual image range. Post-processing can involve normalizing the output or mapping the values back to the valid range.

//Gaussian Filter: The Gaussian filter blurs an image and reduces high-frequency noise. It does not usually require post-processing, but if the filter parameters are not properly chosen, it may lead to overblurring or loss of image details.

//Gradient Filters (GX and GY): The gradient filters compute the gradient magnitude and angle of an image. They typically do not require post-processing, but if the gradient values need to be used for specific purposes (e.g., edge detection), additional processing may be needed.
float gaussian(float x, float y, float sigma)
{
    float twoSigmaSq = 2 * sigma * sigma;
    float coefficient = 1.0 / (M_PI * twoSigmaSq);
    float exponent = -(x * x + y * y) / twoSigmaSq;
    return coefficient * exp(exponent);
}
image make_gaussian_filter(float sigma)
{
    int filterSize = (int)(6 * sigma + 1); // Determine filter size based on sigma
    int center = filterSize / 2; // Center position of the filter

    image filter = make_image(filterSize, filterSize, 1);

    // Compute the Gaussian filter values
    for (int i = 0; i < filterSize; i++) {
        for (int j = 0; j < filterSize; j++) {
            float x = i - center;
            float y = j - center;
            filter.data[i * filterSize + j] = gaussian(x, y, sigma);
        }
    }

    return filter;
}

image add_image(image a, image b)
{
    if (a.width != b.width || a.height != b.height || a.channels != b.channels) {
        printf("Error: Images have incompatible dimensions.\n");
        exit(1);
    }

    // Create the output image with the same dimensions
    image result = make_image(a.width, a.height, a.channels);

    // Perform element-wise addition of the images
    for (int i = 0; i < a.width * a.height * a.channels; i++) {
        result.data[i] = a.data[i] + b.data[i];
    }

    return result;
}

image sub_image(image a, image b)
{
    if (a.width != b.width || a.height != b.height || a.channels != b.channels) {
        printf("Error: Images have incompatible dimensions.\n");
        exit(1);
    }

    // Create the output image with the same dimensions
    image result = make_image(a.width, a.height, a.channels);

    // Perform element-wise subtraction of the images
    for (int i = 0; i < a.width * a.height * a.channels; i++) {
        result.data[i] = a.data[i] - b.data[i];
    }

    return result;
}

image make_gx_filter()
{
    int filterSize = 3;
    image filter = make_image(filterSize, filterSize, 1);

    // Define the GX filter values for edge detection in the x-direction
    float filterData[9] = {-1, 0, 1, -2, 0, 2, -1, 0, 1};

    // Set the filter values
    for (int i = 0; i < filterSize * filterSize; i++) {
        filter.data[i] = filterData[i];
    }

    return filter;
}

image make_gy_filter()
{
    int filterSize = 3;
    image filter = make_image(filterSize, filterSize, 1);

    // Define the GY filter values for edge detection in the y-direction
    float filterData[9] = {-1, -2, -1, 0, 0, 0, 1, 2, 1};

    // Set the filter values
    for (int i = 0; i < filterSize * filterSize; i++) {
        filter.data[i] = filterData[i];
    }

    return filter;
}

void feature_normalize(image im)
{
    float minVal = im.data[0];
    float maxVal = im.data[0];

    // Find the minimum and maximum values in the image
    for (int i = 0; i < im.width * im.height * im.channels; i++) {
        if (im.data[i] < minVal) {
            minVal = im.data[i];
        }
        if (im.data[i] > maxVal) {
            maxVal = im.data[i];
        }
    }

    // Normalize the image values
    float range = maxVal - minVal;
    for (int i = 0; i < im.width * im.height * im.channels; i++) {
        im.data[i] = (im.data[i] - minVal) / range;
    }
}

image *sobel_image(image im)
{
    // Create the sobel filters
    image gxFilter = make_image(3, 3, 1);
    image gyFilter = make_image(3, 3, 1);

    // Define the sobel filters for edge detection
    float gxFilterData[9] = {-1, 0, 1, -2, 0, 2, -1, 0, 1};
    float gyFilterData[9] = {-1, -2, -1, 0, 0, 0, 1, 2, 1};

    // Set the filter values
    for (int i = 0; i < 9; i++) {
        gxFilter.data[i] = gxFilterData[i];
        gyFilter.data[i] = gyFilterData[i];
    }

    // Create the output images for gradient magnitude and angle
    image *output = calloc(2, sizeof(image));
    output[0] = make_image(im.width, im.height, 1); // Gradient magnitude image
    output[1] = make_image(im.width, im.height, 1); // Gradient angle image
    // Perform convolution with the sobel filters
    for (int y = 1; y < im.height - 1; y++) {
        for (int x = 1; x < im.width - 1; x++) {
            float gx = 0.0, gy = 0.0;

            // Convolve with gxFilter
            for (int fy = -1; fy <= 1; fy++) {
                for (int fx = -1; fx <= 1; fx++) {
                    int pixelIdx = (y + fy) * im.width + (x + fx);
                    gx += im.data[pixelIdx] * gxFilter.data[(fy + 1) * 3 + (fx + 1)];
                }
            }

            // Convolve with gyFilter
            for (int fy = -1; fy <= 1; fy++) {
                for (int fx = -1; fx <= 1; fx++) {
                    int pixelIdx = (y + fy) * im.width + (x + fx);
                    gy += im.data[pixelIdx] * gyFilter.data[(fy + 1) * 3 + (fx + 1)];
                }
            }
            float magnitude = sqrt(gx * gx + gy * gy);
            float angle = atan2(gy, gx);

            // Store the results in the output images
            output[0].data[y * im.width + x] = magnitude;
            output[1].data[y * im.width + x] = angle;
        }
    }

    // Free the memory used by the sobel filters
    free(gxFilter.data);
    free(gyFilter.data);

    return output;
}

image colorize_sobel(image im)
{
    // Create the output image with RGB channels
    image result = make_image(im.width, im.height, 3);

    // Iterate over the image pixels
    for (int i = 0; i < im.width * im.height; i++) {
        float magnitude = im.data[i];

        // Compute the hue value based on the angle
        float angle = im.data[i + im.width * im.height];
        float hue = (angle + M_PI) / (2 * M_PI);

        // Set the RGB values based on the hue and magnitude
        result.data[3 * i] = hue;
        result.data[3 * i + 1] = magnitude;
        result.data[3 * i + 2] = magnitude;
    }
 return result
}
