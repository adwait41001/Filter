#include "helpers.h"
#include <math.h>
#include <stdio.h>

// Convert image to grayscale
void grayscale(int height, int width, RGBTRIPLE image[height][width])
{

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            RGBTRIPLE pixel; //declares an RGB triple that constitutes one pixel.
            pixel = image[i][j]; //assigns a pixel to the variable pixel by locating it using its height and width.
            int avg_value = 0;
            avg_value = round((pixel.rgbtRed + pixel.rgbtBlue + pixel.rgbtGreen)/(3.0)); //Calculates the average of the Red,Green and Blue values.
            image[i][j].rgbtRed = avg_value;   //Sets the value of RED to the average
            image[i][j].rgbtBlue = avg_value;  //Sets the value of BLUE to the average
            image[i][j].rgbtGreen = avg_value; //Sets the value of GREEN to the average
        }
    }

    return;
}

// Reflect image horizontally
void reflect(int height, int width, RGBTRIPLE image[height][width])
{
    for (int i = 0; i < height ; i++)
    {
        for (int j = 0; j < width/2; j++)
        {
            RGBTRIPLE temp;
            temp = image[i][j]; // selects a pixel located in the ith row and the jth column and stores it in temporary memory for future use for swapping.
            image[i][j] = image[i][width - j - 1]; //swaps the ith element of the jth column with the width - i - 1 th element of the jth column
            image[i][width - j - 1] = temp;
        }
    }
    return;
}

// Blur Image

int *find_neighbours_sum(int i, int j, int height, int width,
                         RGBTRIPLE image[height][width])  // find sum of neighbour pixel values for all the channels respectively
{
    int sum = 0;
    int count = 0;
    static int return_values[4];
    enum ret_values {Red_values, Green_values, Blue_values, pixel_count};

    return_values[Red_values] = 0;
    return_values[Green_values] = 0;
    return_values[Blue_values] = 0;
    return_values[pixel_count] = 0;

    for (int k = -1; k <= 1; k++)
    {
        if (i + k >= 0 && i + k < height)
        {
            for (int l = -1; l <= 1; l++)
            {
                if (j + l >= 0 && j + l < width)
                {
                    count++;
                    return_values[Red_values] += image[i + k][j + l].rgbtRed;
                    return_values[Green_values] += image[i + k][j + l].rgbtGreen;
                    return_values[Blue_values] += image[i + k][j + l].rgbtBlue;
                }
                else
                {
                    continue;
                }
            }
        }
        else
        {
            continue;
        }
    }

    return_values[pixel_count] = count;

    return return_values;
}

void blur(int height, int width, RGBTRIPLE image[height][width])
{
    RGBTRIPLE copy[height][width];
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            copy[i][j] = image[i][j]; //copies the image to another image
        }
    }

    int *neighbours; //declares a pointer to store the address of the sum of neighbours
    enum ret_values {Red_values, Green_values, Blue_values, pixel_count};

    for (int i = 0; i < height; i++)
    {

        for (int j = 0; j < width; j++)
        {
            neighbours = find_neighbours_sum(i, j, height, width, image); // stores the address of the sum of neighbours in neighbours.

            copy[i][j].rgbtRed = round(neighbours[Red_values] / (float)neighbours[pixel_count]);
            copy[i][j].rgbtGreen = round(neighbours[Green_values] / (float)neighbours[pixel_count]);
            copy[i][j].rgbtBlue = round(neighbours[Blue_values] / (float)neighbours[pixel_count]);
        }
    }

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            image[i][j] = copy[i][j];
        }
    }

    return;
}

// Edges

void sobel_operator(char direction, int sobel[3][3]) // Select x or y direction sobel operator
{
    int sobel_x[3][3] = {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};
    int sobel_y[3][3] = {{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}};
    if (direction == 'x') //Sobel operator for x direction
    {
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                sobel[i][j] = sobel_x[i][j];
            }
        }
        return;
    }
    else // sobel operator for y direction
    {
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                sobel[i][j] = sobel_y[i][j];
            }
        }
        return;
    }
}

void find_sobel_product(int sobel[3][3], int neighbours[3][3][3],
                        int sobel_product[3][3][3]) // Find product of sobel operator with neighbours
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                sobel_product[i][j][k] = neighbours[i][j][k] * sobel[j][k];
            }
        }
    }
    return;
}

void edges(int height, int width, RGBTRIPLE image[height][width])
{

    int sobel_x[3][3];  // For sobel operator in x direction
    int sobel_y[3][3];  // For sobel operator in y direction

    // For all 3d storage arrays used below, [0][0][0] corresponds to 0th pixel's red channel, [1][0][0] corresponds to 0th pixel's green channel and [2][0][0] corresponds to 0th pixel's blue channel


    int neighbours[3][3][3];  // For neighbour pixel value storage
    int sobel_product_x[3][3][3];
    int sobel_product_y[3][3][3];
    int final_gx_red = 0;  // Final red channel value for a pixel in x direction
    int final_gx_blue = 0;
    int final_gx_green = 0;
    int final_gy_red = 0;  // Final red channel value for a pixel in y direction
    int final_gy_green = 0;
    int final_gy_blue = 0;
    int final_red = 0;  // Final red channel value computed from x and y directions
    int final_green = 0;
    int final_blue = 0;

    RGBTRIPLE copy[height][width];

    for (int i = 0; i < height; i++) // Copy to a temporary variable image
    {
        for (int j = 0; j < width; j++)
        {
            copy[i][j] = image[i][j];
        }
    }

    for (int i = 0; i < height; i++) // Loop over the image pixel by pixel
    {
        for (int j = 0; j < width; j++)
        {
            final_gx_red = 0;
            final_gx_blue = 0;
            final_gx_green = 0;
            final_gy_red = 0;
            final_gy_green = 0;
            final_gy_blue = 0;
            final_red = 0;
            final_green = 0;
            final_blue = 0;
            sobel_operator('x', sobel_x);  // Get sobel operator for x direction
            sobel_operator('y', sobel_y);  // Get sobel operator for y direction

            for (int k = -1; k <= 1; k++) // A double loop to find neighbouring pixels of a pixel
            {
                for (int l = -1; l <= 1; l++)
                {
                    if (i + k < 0 || i + k >= height || j + l < 0 || j + l >= width)  // Violating boundary conditions
                    {
                        neighbours[0][k + 1][l + 1] = 0;  // Handle edge cases for red channel
                        neighbours[1][k + 1][l + 1] = 0;
                        neighbours[2][k + 1][l + 1] = 0;
                    }
                    else
                    {
                        neighbours[0][k + 1][l + 1] = image[i + k][j + l].rgbtRed;  // Store red channel value of neighbour pixel
                        neighbours[1][k + 1][l + 1] = image[i + k][j + l].rgbtGreen;
                        neighbours[2][k + 1][l + 1] = image[i + k][j + l].rgbtBlue;
                    }
                }
            }

            find_sobel_product(sobel_x, neighbours, sobel_product_x); //Find product of neighbours and sobel operator
            find_sobel_product(sobel_y, neighbours, sobel_product_y);

            for (int s = 0; s < 3; s++) // Compute sum final Gx and Gy value for a channel
            {
                for (int t = 0; t < 3; t++)
                {
                    final_gx_red += sobel_product_x[0][s][t];
                    final_gx_green += sobel_product_x[1][s][t];
                    final_gx_blue += sobel_product_x[2][s][t];
                    final_gy_red += sobel_product_y[0][s][t];
                    final_gy_green += sobel_product_y[1][s][t];
                    final_gy_blue += sobel_product_y[2][s][t];
                }
            }

            final_red = round(sqrt(pow(final_gx_red, 2) + pow(final_gy_red,
                                   2))); // Compute final channel value as per formula sqrt(Gx^2 + Gy^2)
            final_green = round(sqrt(pow(final_gx_green, 2) + pow(final_gy_green, 2)));
            final_blue = round(sqrt(pow(final_gx_blue, 2) + pow(final_gy_blue, 2)));

            if (final_red > 255) //Cap each channel value at 255
            {
                final_red = 255;
            }
            if (final_green > 255)
            {
                final_green = 255;
            }
            if (final_blue > 255)
            {
                final_blue = 255;
            }


            copy[i][j].rgbtRed = final_red;
            copy[i][j].rgbtGreen = final_green;
            copy[i][j].rgbtBlue = final_blue;
        }
    }

    for (int i = 0; i < height; i++) //Copy new image to original
    {
        for (int j = 0; j < width; j++)
        {
            image[i][j] = copy[i][j];
        }
    }

    return;
}