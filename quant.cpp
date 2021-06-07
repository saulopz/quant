#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>

typedef unsigned char byte;
typedef unsigned int u_int;
typedef unsigned short int word;
typedef unsigned long u_long;

// 5 bits to red, 5 bits to green and 5 bits to blue.
// This macro convert RGB to a long int to create a
// otimized palette not considering very close colors.
#define RGB(r, g, b) (word)(((b) & ~7) << 7) | (((g) & ~7) << 2) | ((r) >> 3)

// this macro converts from the optimized palette to RGB
#define RED(x) (byte)(((x)&31) << 3)
#define GREEN(x) (byte)((((x) >> 5) & 255) << 3)
#define BLUE(x) (byte)((((x) >> 10) & 255) << 3)

// Returns only the most significant 5 bits and cancels the least 3
#define SIG5(x) (byte)((x)&248)

// Colors vector ******************************************************

enum
{
    RED = 0,   // Red color vector
    GREEN = 1, // Green color vector
    BLUE = 2,  // Blue color vector
};

// Frequency histogram ***********************************************

struct Hist
{
    byte red, green, blue;
    u_long freq;
    float erro;
} * histogram;

u_long hist_size;

// Image information **********************************************

int Width;                // Image width
int Height;               // Image height
int Type;                 // Image type (0=RGB, 1=MAP)
byte *Red, *Green, *Blue; // Image RGB
byte *Map;                // Image Map
long *Colors;             // Image color table

int err;
double quantization_error;

// Timer information ************************************************

typedef long int __time_t;
typedef __time_t time_t;

struct tm *the_time;
time_t t;

typedef struct t_dec
{
    int sec;
    int min;
    int hour;
} tmp;

static tmp t_ini, t_fin, t_tot;

// Frequency histogram with float colors instead of byte ****************

struct Hist_Float
{
    float red, green, blue;
    u_long freq;
    float error;
} * hist_float;

u_long new_hist_size;

// Image functions ************************************************************************

SDL_Surface *load_image(char *filename, int *Width, int *Height);
int load_image_rgb(SDL_Surface *surface, byte *Red, byte *Green, byte *Blue);
int save_image_rgb(SDL_Surface *surface, byte *Red, byte *Green, byte *Blue, char *filename);

// Functions for general use *************************************************************

void error_message(int err);
char *str_to_upper(char *s);
int check_parameters(char *Param);
void save_time(tmp &T);
void calculate_time(void);
void information(void);
void create_histgram(long size, byte *red, byte *green, byte *blue);
double euclidian_distance(Hist_Float ci, Hist_Float cj);
double euclidian_distance(Hist ci, Hist cj);

// Information for the Double Cluster Quantization Algorithm ****************

// Error matrix information - E[i][j]
float *coluna;
float **E;

// Functions for the Double Cluster Quantization Algorithm *********************

void double_cluster(int width, int height, byte *red, byte *green, byte *blue, u_long newpalsize);
void select_color_pairs(u_long &i_min, u_long &j_min);
void calculate_error_matrix(u_long tam_matriz);
void join_color_pairs(u_long ci, u_long cj, Hist_Float &cor);
void change_colors(u_long ci, u_long cj, Hist_Float cor);
void apply_quantization(u_long tam_imagem, byte *red, byte *green, byte *blue);

// Vector used in the Recursive Subdivision Algorithm *****************

struct VET
{
    byte red, green, blue;
    u_long bigger, smaller;
    u_long freq;
} * vector;

// Functions for the Recursive Subdivision Quantization Algorithm

void change_histogram_position(byte *color, u_int pos1, u_int pos2);
void sort(byte *color, u_long esq, u_long dir);
void recursive_subdivision(int width, int height, byte *red,
                           byte *green, byte *blue,
                           u_long newpalsize);

// Information for the octree *****************************************

#define prof_max = 8; // Maximum depth to the tree

// Information of a octree node
typedef struct _Node
{
    struct _Node *father, *child[8];
    float r, g, b;   // colors red, green and blue
    byte er, eg, eb; // smaller colors of subcube
    byte dr, dg, db; // bigger colors of subcube
    byte level;      // level where is node in the octree
    u_long n_colors; // quantity of colores
} Node;

// RGB cube information
typedef struct _RGBCube
{
    Node *root;            // Root node
    u_long idx_colors;     // Size colors of RGB cube
    u_long maximum_colors; // Maximun size of colors
    Node **list;           // Leafs list of RGB cube
} RGBCube;

static RGBCube cubo;

// Functions for the Octree Quantization Algorithm **************************

void octree(int width, int height, byte *red, byte *green, byte *blue,
            u_long newpalsize);
void inicialize_octree(u_long qt_cores);
void reduces_octree(void);
void insert_node(Node *no_pai, byte index, byte red, byte green, byte blue);
void track_octree(byte red, byte green, byte blue);
void subcube_area(Node *pont);

// Main function ***************************************************************************

int main(int argc, char *argv[])
{
    u_long quantcolors = 0;
    if (argc < 5)
    {
        information();
        return 0;
    }
    SDL_Surface *surface = load_image(argv[1], &Width, &Height);
    if (surface)
    {
        printf("Image information %s\n", argv[1]);
        printf("Width: %d, Height: %d\n", Width, Height);
        Red = (byte *)malloc(Width * Height);
        Green = (byte *)malloc(Width * Height);
        Blue = (byte *)malloc(Width * Height);
        if (!Red | !Green | !Blue)
        {
            printf("Alloc memory error.");
            return 0;
        }
        err = load_image_rgb(surface, Red, Green, Blue);
        if (!err)
        {
            if (argc > 3)
            {
                quantcolors = atol(argv[3]);
                save_time(t_ini);
                switch (atoi(argv[4]))
                {
                case 0:
                    recursive_subdivision(Width, Height, Red, Green, Blue, quantcolors);
                    break;
                case 1:
                    octree(Width, Height, Red, Green, Blue, quantcolors);
                    break;
                case 2:
                    double_cluster(Width, Height, Red, Green, Blue, quantcolors);
                    break;
                default:
                    information();
                };
                save_time(t_fin);
                calculate_time();
                err = save_image_rgb(surface, Red, Green, Blue, argv[2]);
                if (err)
                    return 1;
                free(Red);
                free(Green);
                free(Blue);
            }
        }
        else
            return 2;
    }
    else
        return 3;
    return 0;
}

// Image functions **********************************************************************

SDL_Surface *load_image(char *filename, int *Width, int *Height)
{
    SDL_Surface *surface = IMG_Load(filename);
    if (surface)
    {
        *Width = surface->w;
        *Height = surface->h;
        // Only 32 bits image
        if (surface->format->BitsPerPixel != 32)
        {
            fprintf(stderr, "Image has %d-bits. Use only 32-bit images.\n",
                    surface->format->BitsPerPixel);
            SDL_FreeSurface(surface);
            return nullptr;
        }
    }
    return surface;
}

int load_image_rgb(SDL_Surface *surface, byte *Red, byte *Green, byte *Blue)
{
    if (!surface)
    {
        return 1;
    }
    SDL_Color *pixels = (SDL_Color *)surface->pixels;
    for (int i = 0; i < surface->w * surface->h; i++)
    {
        Red[i] = pixels[i].r;
        Green[i] = pixels[i].g;
        Blue[i] = pixels[i].b;
    }
    return 0;
}

int save_image_rgb(SDL_Surface *surface, byte *Red, byte *Green, byte *Blue, char *filename)
{
    if (!surface)
    {
        return 1;
    }
    SDL_Color *pixels = (SDL_Color *)surface->pixels;
    for (int i = 0; i < surface->w * surface->h; i++)
    {
        pixels[i].r = Red[i];
        pixels[i].g = Green[i];
        pixels[i].b = Blue[i];
    }
    if (SDL_SaveBMP(surface, filename) != 0)
    {
        printf("SDL_SaveBMP failed: %s\n", SDL_GetError());
    }
    return 0;
}

// Functions for general use ***********************************************************

// Transform a string to uppercase

char *str_to_upper(char *s)
{
    int i, len = strlen(s);

    for (i = 0; i < len; i++)
        s[i] = toupper(s[i]);

    return s;
}

void save_time(tmp &T)
{
    t = time(NULL);
    the_time = localtime(&t);
    T.sec = the_time->tm_sec;
    T.min = the_time->tm_min;
    T.hour = the_time->tm_hour;
}

// Calculate execution time of an algorithm

void calculate_time(void)
{

    if (t_fin.sec < t_ini.sec)
    {
        t_fin.min--;
        t_fin.sec += 60;
    }
    t_tot.sec = t_fin.sec - t_ini.sec;

    if (t_fin.min < t_ini.min)
    {
        t_fin.hour--;
        t_fin.min += 60;
    }
    t_tot.min = t_fin.min - t_ini.min;
    t_tot.hour = t_fin.hour - t_ini.hour;

    printf("Elapsed time: %d:%d:%d\n", t_tot.hour, t_tot.min, t_tot.sec);
}

/***************************************************************** 
 * If not enter with correct parameters, shows on the screen
 * information about program utilization.
******************************************************************/

void information(void)
{
    printf("IMAGE QUANTIZATION                                   \n");
    printf("                                                     \n");
    printf("This program is part of a Bachelor's Degree Work in  \n");
    printf("Computer Science from UNOESC in 1998.                \n");
    printf("                                                     \n");
    printf("ACADEMIC: Saulo Popov Zambiasi.                      \n");
    printf("  E-MAIL: saulopz@gmail.com                          \n");
    printf(" ADVISOR: Marcos Vinicius Rayol Sobreiro             \n");
    printf("                                                     \n");
    printf("This program aims to transform an image from N colors\n");
    printf("to M colors, where N > M, calculating the elapsed    \n");
    printf("time. Three color image quantization algorithms are  \n");
    printf("used.                                                \n");
    printf("                                                     \n");
    printf("Utilization:                                         \n");
    printf("                                                     \n");
    printf("quant <img1.ext> <img2.ext> <colors> <algorithm>     \n");
    printf("                                                     \n");
    printf("   img1.ext  - original image with N colors.         \n");
    printf("   img2.ext  - result image with M colors (N > M).   \n");
    printf("   colors    - quantity of colors in img2.ext.       \n");
    printf("   algorithm - quantization algorithm:               \n");
    printf("      0 - Recursive Subdivision.                     \n");
    printf("      1 - Octree.                                    \n");
    printf("      2 - Double Clusters.                           \n");
    printf("                                                     \n");
    printf(" Example: quant image.png result.png 256 0           \n");
}

/***************************************************************** 
 * Calculate the number of colors existing on image and generates
 * the frequency histogram.
 * Parameters:
 *    long size	 - image size (width * height
 *    byte red   - vector for red color
 *    byte green - vector for green color
 *    byte blue  - vector for blue color
 * ***************************************************************/

void create_histgram(long size, byte *red, byte *green, byte *blue)
{
    long pimage, phist, qt;
    int pcent, old;
    bool exist;
    word *cor;
    u_long *Freq;

    cor = (word *)malloc(1000000);
    Freq = (u_long *)malloc(1000000);

    /******************************************************
     * First quantizes the image to 15 bits:
     *   red   - 5 bits
     *   green - 5 bits
     *   blue  - 5 bits
     * and calculate frequency histogram size
     * ****************************************************/

    pcent = old = 0;
    printf("Quantizes to 15 bits. Computing colors...\n");
    qt = 0;
    for (pimage = 0; pimage < size; pimage++)
    {
        phist = 0;
        exist = false;
        while ((!exist) && (phist < qt))
        {
            exist = ((SIG5(red[pimage]) == RED(cor[phist])) &&
                     (SIG5(green[pimage]) == GREEN(cor[phist])) &&
                     (SIG5(blue[pimage]) == BLUE(cor[phist])));
            phist++;
        }
        if (!exist)
        {
            cor[qt] = RGB(red[pimage], green[pimage], blue[pimage]);
            Freq[qt] = 1;
            qt++;
        }
        else
            Freq[phist - 1]++;

        pcent = ((100 * pimage) / size) / 2;
        if (pcent > old)
        {
            old = pcent;
            printf(".");
        }
    }
    printf("%ld\n", qt);
    hist_size = qt;

    // Creating the Frequency Histogram
    pcent = old = 0;
    printf("Gerar Histograma");

    histogram = (struct Hist *)malloc((hist_size) * sizeof(struct Hist));
    if (!histogram)
    {
        printf("...insuficient memory...\n");
        return;
    }
    old = 0;
    for (phist = 0; phist < qt; phist++)
    {
        histogram[phist].red = RED(cor[phist]);
        histogram[phist].green = GREEN(cor[phist]);
        histogram[phist].blue = BLUE(cor[phist]);
        histogram[phist].freq = Freq[phist];
        histogram[phist].erro = 0;

        pcent = (50 * phist) / hist_size;
        if (pcent > old)
        {
            old = pcent;
            printf(".");
        }
    }
    printf("Ok\n");
    free(cor);
    free(Freq);
}
/*********************************************************************
 * Calculate the distance of two colors usinng Euclidean metric square
 * 
 *    d(ci, cj) = ||ci - cj||^ 2.
 * 
 * *******************************************************************/

double euclidian_distance(Hist_Float ci, Hist_Float cj)
{
    double dist;

    dist = pow(ci.red - cj.red, 2) +
           pow(ci.green - cj.green, 2) +
           pow(ci.blue - cj.blue, 2);
    return ((double)dist);
}

double euclidian_distance(Hist ci, Hist cj)
{
    double dist;

    dist = pow(ci.red - cj.red, 2) +
           pow(ci.green - cj.green, 2) +
           pow(ci.blue - cj.blue, 2);
    return ((double)dist);
}

/******************************************************************
 * Double Cluster Quantization Function
 * Parameters:
 *   int width         - image width
 *   int height        - image height
 *   byte red          - vector for red color
 *   byte green        - vector for green color
 *   byte blue         - vector for blue color
 *   u_long newpalsize - image size of colors
 * ****************************************************************/

void double_cluster(int width, int height, byte *red, byte *green, byte *blue, u_long newpalsize)
{
    u_long tam_imagem, i;
    tam_imagem = width * height;
    printf("Quantization for Double Clusters.\n");

    // Compute quantity of colors on image and generate the histogram
    create_histgram(width * height, red, green, blue);
    new_hist_size = hist_size;

    // Create the histogram with float values (for calcuations purpose).
    hist_float = (struct Hist_Float *)malloc((new_hist_size) * sizeof(struct Hist_Float));
    for (i = 0; i < new_hist_size; i++)
    {
        hist_float[i].red = histogram[i].red;
        hist_float[i].green = histogram[i].green;
        hist_float[i].blue = histogram[i].blue;
        hist_float[i].freq = histogram[i].freq;
    }

    // If new quantity of colors of new image is bigger or equals
    // the quantity of original image, then is not necessary operation
    if (new_hist_size <= newpalsize)
    {
        printf("Nao e necessario quantizar\n");
        exit(1);
    }

    u_long ci, cj;
    Hist_Float cor;
    int pcent, old = 0;

    // Compute the error matrix[new_hist_size][new_hist_size]
    calculate_error_matrix(new_hist_size);

    printf("Quantizing...");

    while (new_hist_size > newpalsize)
    {
        select_color_pairs(ci, cj);
        join_color_pairs(ci, cj, cor);
        change_colors(ci, cj, cor);

        pcent = (int)(((new_hist_size - newpalsize) * 50) / hist_size - newpalsize);
        if (pcent < old)
        {
            printf(".");
            old = pcent;
        } /* if */
    }
    printf("Ok\n");
    apply_quantization(tam_imagem, red, green, blue);

    free(histogram);
}

// Select the pair of colors (ci, cj) with minimum error
// of quantizing Ei,j in matrix E

void select_color_pairs(u_long &i_min, u_long &j_min)
{
    u_long i, j, im, jm;
    float Em;

    im = 0;
    jm = 1;
    Em = E[im][jm];

    for (i = 0; i < new_hist_size; i++)
    {
        for (j = 0; j < new_hist_size; j++)
        {
            if ((i < j) && (Em > E[i][j]))
            {
                im = i;
                jm = j;
                Em = E[i][j];
            }
        }
    }
    i_min = im;
    j_min = jm;
}

void join_color_pairs(u_long ci, u_long cj, Hist_Float &cor)
{
    float a, h1, h2;
    a = (float)hist_float[ci].freq + hist_float[cj].freq;
    h1 = hist_float[ci].freq / a;
    h2 = hist_float[cj].freq / a;

    cor.red = ((h1 * hist_float[ci].red) + (h2 * hist_float[cj].red));
    cor.green = ((h1 * hist_float[ci].green) + (h2 * hist_float[cj].green));
    cor.blue = ((h1 * hist_float[ci].blue) + (h2 * hist_float[cj].blue));
    cor.freq = hist_float[ci].freq + hist_float[cj].freq;
}

void change_colors(u_long ci, u_long cj, Hist_Float cor)
{
    u_long i, j;

    // I put the new color on the color ci place
    hist_float[ci].red = cor.red;
    hist_float[ci].green = cor.green;
    hist_float[ci].blue = cor.blue;
    hist_float[ci].freq = cor.freq;

    // I decrease the size of histogram
    new_hist_size--;

    // I remove the color cj of histogram
    for (i = cj; i < new_hist_size; i++)
    {
        hist_float[i] = hist_float[i + 1];
    };

    // I remove the color cj of error matrix, line and column
    for (i = cj; i < new_hist_size; i++)
    {
        for (j = 0; j < new_hist_size; j++)
        {
            E[i][j] = E[i + 1][j];
        }
    }
    /* free(E[tam_hist_new]); */
    for (i = 0; i < new_hist_size; i++)
    {
        for (j = cj; j < new_hist_size; j++)
        {
            E[i][j] = E[i][j + 1];
        }
        /* E[i][tam_hist_new] = NULL; */
    }

    // Recalculate the error in ci, line and column
    j = ci;
    for (i = 0; i < new_hist_size; i++)
    {
        if (i < j)
        {
            E[i][j] = (float)(((hist_float[i].freq + pow(hist_float[j].freq, 2)) +
                               (hist_float[j].freq + pow(hist_float[i].freq, 2))) /
                              pow(hist_float[i].freq + hist_float[j].freq, 2)) *
                      (float)euclidian_distance(hist_float[i], hist_float[j]);
        }
        else
            E[i][j] = 0;
        if (j < i)
        {
            E[j][i] = (float)(((hist_float[i].freq + pow(hist_float[j].freq, 2)) +
                               (hist_float[j].freq + pow(hist_float[i].freq, 2))) /
                              pow(hist_float[i].freq + hist_float[j].freq, 2)) *
                      (float)euclidian_distance(hist_float[i], hist_float[j]);
        }
        else
            E[j][i] = 0;
    }
}

void calculate_error_matrix(u_long tam_matriz)
{
    u_long i, j;
    int pcent, old = 0;

    // Alloc memory to a matrix with matrix_zise x matrix_size
    printf("Matrix error");

    E = (float **)malloc(sizeof(E) * tam_matriz);
    if (!E)
    {
        printf("\nNot possible alloc memory to error matrix.");
        exit(1);
    }
    for (i = 0; i < tam_matriz; i++)
    {
        E[i] = (float *)malloc(sizeof(coluna) * tam_matriz);
        if (!E[i])
        {
            printf("\nNot possible alloc memory to error matrix.");
            exit(1);
        }
    }

    // compute the errors of matrix
    for (i = 0; i < tam_matriz; i++)
    {
        for (j = 0; j < tam_matriz; j++)
        {
            if (i < j)
            {
                E[i][j] = (float)(((hist_float[i].freq + pow(hist_float[j].freq, 2)) +
                                   (hist_float[j].freq + pow(hist_float[i].freq, 2))) /
                                  pow(hist_float[i].freq + hist_float[j].freq, 2)) *
                          (float)euclidian_distance(hist_float[i], hist_float[j]);
            }
            else
                E[i][j] = 0;
        }
        pcent = (int)((i * 50) / tam_matriz);
        if (pcent > old)
        {
            printf(".");
            old = pcent;
        }
    }
    printf("Ok\n");
}

void apply_quantization(u_long tam_imagem, byte *red, byte *green, byte *blue)
{
    u_long i, j;
    Hist_Float *hist2;
    int pcent, old = 0;

    hist2 = (struct Hist_Float *)malloc((hist_size) * sizeof(struct Hist_Float));
    for (i = 0; i < hist_size; i++)
    {
        hist2[i].red = histogram[i].red;
        hist2[i].green = histogram[i].green;
        hist2[i].blue = histogram[i].blue;
        hist2[i].freq = histogram[i].freq;
    }

    u_long menor_dist;
    double dist, dist2;

    printf("Close colors");

    // For each color on histogram
    for (i = 0; i < hist_size; i++)
    {
        // Find the smaller distance
        menor_dist = 0;
        dist = euclidian_distance(hist2[i], hist_float[menor_dist]);
        for (j = 0; j < new_hist_size; j++)
        {
            dist2 = euclidian_distance(hist2[i], hist_float[j]);
            if (dist2 < dist)
            {
                dist = dist2;
                menor_dist = j;
            }
        }
        // Replaces the color by with quantized cell of smaller distance
        hist2[i].red = hist_float[menor_dist].red;
        hist2[i].green = hist_float[menor_dist].green;
        hist2[i].blue = hist_float[menor_dist].blue;
        hist2[i].freq = hist_float[menor_dist].freq;

        pcent = (int)((i * 50) / hist_size);
        if (pcent > old)
        {
            printf(".");
            old = pcent;
        }
    }
    printf("Ok\n");

    // Maps the new frequency histogram to the image
    bool achou;

    printf("Mapping colors");
    old = 0;

    Hist_Float H;

    for (i = 0; i < tam_imagem; i++)
    {
        achou = false;
        j = 0;
        while ((!achou) && (j < hist_size))
        {
            if ((SIG5(red[i]) == histogram[j].red) &&
                (SIG5(green[i]) == histogram[j].green) &&
                (SIG5(blue[i]) == histogram[j].blue))
            {
                H.red = red[i];
                H.green = green[i];
                H.blue = blue[i];
                // Calculate the distance from the color to its quantized color
                dist = (float)(euclidian_distance(H, hist2[j]));
                histogram[j].erro = (float)(dist / histogram[j].freq);

                red[i] = (byte)hist2[j].red;
                green[i] = (byte)hist2[j].green;
                blue[i] = (byte)hist2[j].blue;
                achou = true;
            }
            j++;
        }
        pcent = (int)((i * 50) / tam_imagem);
        if (pcent > old)
        {
            printf(".");
            old = pcent;
        } /* if */
    }     /* for */
    printf("Ok\n");

    printf("Quantization error");
    quantization_error = 0;
    old = 0;

    for (i = 0; i < hist_size; i++)
    {
        quantization_error += histogram[i].erro;

        pcent = (int)((i * 50) / hist_size);
        if (pcent > old)
        {
            printf(".");
            old = pcent;
        } /* if */
    }
    quantization_error = quantization_error / hist_size;
    printf("%f\n", quantization_error);

    free(hist2);
}

/*****************************************************************
 * Change two elements of histogram and of color vector in a
 * specific position
 * ***************************************************************/

void change_histogram_position(byte *color, u_int pos1, u_int pos2)
{
    byte c;
    byte red;
    byte green;
    byte blue;
    float erro;
    u_long freq;

    c = color[pos1];
    red = histogram[pos1].red;
    green = histogram[pos1].green;
    blue = histogram[pos1].blue;
    freq = histogram[pos1].freq;
    erro = histogram[pos1].erro;

    color[pos1] = color[pos2];
    histogram[pos1].red = histogram[pos2].red;
    histogram[pos1].green = histogram[pos2].green;
    histogram[pos1].blue = histogram[pos2].blue;
    histogram[pos1].freq = histogram[pos2].freq;
    histogram[pos1].erro = histogram[pos2].erro;

    color[pos2] = c;
    histogram[pos2].red = red;
    histogram[pos2].green = green;
    histogram[pos2].blue = blue;
    histogram[pos2].freq = freq;
    histogram[pos2].erro = erro;
}

/*****************************************************************
 * Sort a sequence on frequency histogram, idexing by a color.
 * The algorithm used is the selection ordenation.
 * ***************************************************************/

void sort(byte *color, u_long esq, u_long dir)
{
    u_long i, j, min;

    for (i = esq; i < dir; i++)
    {
        min = i;
        for (j = i + 1; j <= dir; j++)
        {
            if (color[j] < color[min])
                min = j;
        }
        change_histogram_position(color, min, i);
    }
}

/******************************************************************
 * Find the bigger vector from left to right.
 * returns:
 *   RED   - if bigger vector is RED
 *   GREEN - if bigger vector is GREEN
 *   BLUE  - if bigger vector is BLUE
 * ****************************************************************/

int bigger_vector(u_int esq, u_int dir)
{
    u_long size_red, size_green, size_blue, i;
    byte min_red, min_green, min_blue;
    byte max_red, max_green, max_blue;

    max_red = histogram[esq].red;
    max_green = histogram[esq].green;
    max_blue = histogram[esq].blue;

    min_red = histogram[esq].red;
    min_green = histogram[esq].green;
    min_blue = histogram[esq].blue;

    for (i = esq; i <= dir; i++)
    {
        if (histogram[i].red > max_red)
            max_red = histogram[i].red;
        if (histogram[i].green > max_green)
            max_green = histogram[i].green;
        if (histogram[i].blue > max_blue)
            max_blue = histogram[i].blue;

        if (histogram[i].red < min_red)
            min_red = histogram[i].red;
        if (histogram[i].green < min_green)
            min_green = histogram[i].green;
        if (histogram[i].blue < min_blue)
            min_blue = histogram[i].blue;
    }
    size_red = max_red - min_red;
    size_green = max_green - min_green;
    size_blue = max_blue - min_blue;
    if ((size_red > size_green) && (size_red > size_blue))
        return RED;
    else if (size_green > size_blue)
        return GREEN;
    else
        return BLUE;
}

/*******************************************************************
 * Quantization for Recursive Subdivision
 * Inputs:
 *   int width         - image width
 *   int height        - image height
 *   byte red          - vector for red color
 *   byte gree         - vector for green color
 *   byte blue         - vector for blue color
 *   u_long newpalsize - colors quantity of new image
 * ****************************************************************/

void recursive_subdivision(int width, int height, byte *red,
                           byte *green, byte *blue,
                           u_long newpalsize)
{
    u_long ipos, iultimo, ipont, msoma, i, soma, size, image_size;
    int mvec, pcent, old;

    printf("Algorithm of Recursive Subdivision\n");

    // Calculates the colors number on image and create histogram
    image_size = width * height;
    create_histgram(image_size, red, green, blue);

    // If new color numbers is bigger or equals old image, exit
    if (hist_size <= newpalsize)
    {
        printf("It is not necessary to quantize.\n");
        exit(1);
    }
    vector = (struct VET *)malloc((newpalsize) * sizeof(struct VET));
    if (!vector)
    {
        printf("Unable to allocate memory for quantization vector.");
        exit(2);
    }

    ipont = iultimo = 0;
    // Insert first element into the vector
    vector[ipont].bigger = hist_size - 1;
    vector[ipont].smaller = 0;

    printf("Mid cut");
    pcent = old = 0;

    byte *color;
    color = (byte *)malloc(hist_size);
    bool inseriu = true;

    while ((iultimo < newpalsize - 1) && (inseriu))
    {
        ipont = 0;
        ipos = iultimo;
        inseriu = false;
        while ((ipont <= ipos) && (iultimo < newpalsize - 1))
        {
            // Find bigger vector on histogram
            mvec = bigger_vector(vector[ipont].smaller, vector[ipont].bigger);

            // Sort by bigger vector
            for (i = vector[ipont].smaller; i <= vector[ipont].bigger; i++)
            {
                switch (mvec)
                {
                case RED:
                    color[i] = histogram[i].red;
                    break;
                case GREEN:
                    color[i] = histogram[i].green;
                    break;
                case BLUE:
                    color[i] = histogram[i].blue;
                    break;
                default:
                    color[i] = histogram[i].red;
                }
            }
            if (vector[ipont].smaller < vector[ipont].bigger - 1)
            {
                sort(color, vector[ipont].smaller, vector[ipont].bigger);
                // Find middle
                // First find the totam sum
                soma = 0;
                for (i = vector[ipont].smaller; i <= vector[ipont].bigger; i++)
                {
                    soma += histogram[i].freq;
                }
                soma = (u_long)soma / 2;
                // Find middle
                msoma = 0;
                for (i = vector[ipont].smaller; i <= vector[ipont].bigger; i++)
                {
                    if (msoma >= soma)
                        break;
                    msoma += histogram[i].freq;
                }
                // Replaces the current vector with the first
                // Add a second vector
                iultimo++;
                vector[iultimo].bigger = vector[ipont].bigger;
                vector[iultimo].smaller = i;
                vector[ipont].bigger = i;
                inseriu = true;
            }
            ipont++;
            while ((vector[ipont].smaller == vector[ipont].bigger) && (ipont <= ipos))
                ipont++;
            pcent = (iultimo * 50) / newpalsize;
            if (pcent > old)
            {
                printf(".");
                old = pcent;
            }
        }
    }
    //	free(color);
    printf("Ok\nColor average");
    pcent = old = 0;

    float cr, cg, cb;

    // Find the color averages in each cell
    for (ipont = 0; ipont < newpalsize; ipont++)
    {
        cr = cg = cb = 0;
        soma = 0;
        for (i = vector[ipont].smaller; i <= vector[ipont].bigger; i++)
        {
            // Add cell colors
            cr += histogram[i].red * histogram[i].freq;
            cg += histogram[i].green * histogram[i].freq;
            cb += histogram[i].blue * histogram[i].freq;
            soma += histogram[i].freq;
        }
        // Divide the sum of the cells by the amount of colors
        // in it, resulting in the average.
        if (soma < 1)
        {
            soma = 1;
        }
        vector[ipont].red = (byte)(cr / soma);
        vector[ipont].green = (byte)(cg / soma);
        vector[ipont].blue = (byte)(cb / soma);
        vector[ipont].freq = soma;

        pcent = (ipont * 50) / newpalsize;
        if (pcent > old)
        {
            printf(".");
            old = pcent;
        }
    }
    printf("Ok\n");

    /*  THIS PART IS WITH PROBLEMS - I ISOLATED
        UPDATE 2021.06.07:
            And I don't even remember why.
            So let's leave it at that.

    printf("Mapear Bordas"); pcent = old = 0;

	bool primeira;
	int  dist1, dist2;

        // Mapear as cores nas bordas para a celula mais proxima
        //	for (i=0; i<tam_hist; i++){
        // Encontra a primeira referencia
        //		primeira = true;
		for(ipos=0; ipos<newpalsize; ipos++){
			if ((i == vetor[ipos].menor) || (i == vetor[ipos].maior)){
				if (primeira){
					primeira = false;
					ipont = ipos;
					dist1 = 
						abs(histograma[i].red   - vetor[ipos].red)   +
						abs(histograma[i].green - vetor[ipos].green) +
						abs(histograma[i].blue  - vetor[ipos].blue);
				} else {
					dist2 =
						abs(histograma[i].red   - vetor[ipos].red)   +
						abs(histograma[i].green - vetor[ipos].green) +
						abs(histograma[i].blue  - vetor[ipos].blue);
					if (dist2 < dist1){
						if (vetor[ipont].menor < vetor[ipont].maior){
							if (i == vetor[ipont].menor){
								vetor[ipont].menor++;
							} else vetor[ipont].maior--;
							dist1 = dist2;
							ipont = ipos;
						} else if (vetor[ipos].menor < vetor[ipos].maior){
							if (i == vetor[ipos].menor){
								vetor[ipos].menor++;
							} else vetor[ipos].maior--;
						} 
					} else{
						if (vetor[ipos].menor < vetor[ipos].maior){
							if (i == vetor[ipos].menor){
								vetor[ipos].menor++;
							} else vetor[ipos].maior--;
						} else{ 
							if (i == vetor[ipont].menor){
								vetor[ipont].menor++;
							} else vetor[ipont].maior--;
							dist1 = dist2;
							ipont = ipos;
						}					
					}
				}
			}
		}
		pcent = (int) ((i*50)/tam_hist);
		if (pcent > old){
			printf(".");
			old = pcent;
		}
	}
	printf("Ok\n");
    */

    // Create new frequency histogram
    struct Hist *hist2;
    hist2 = (struct Hist *)malloc((hist_size) * sizeof(struct Hist));
    if (!hist2)
    {
        printf("Unable to allocate memory...");
        exit(2);
    }

    // Map Colors to New Frequency Histogram
    for (ipos = 0; ipos < newpalsize; ipos++)
    {
        for (i = vector[ipos].smaller; i <= vector[ipos].bigger; i++)
        {
            hist2[i].red = (byte)vector[ipos].red;
            hist2[i].green = (byte)vector[ipos].green;
            hist2[i].blue = (byte)vector[ipos].blue;
        }
    }

    // Mapping the colors to the new frequency histogram
    printf("Mapping the colors");
    pcent = old = 0;

    bool found;
    size = width * height;
    for (i = 0; i < size; i++)
    {
        found = false;
        ipos = 0;
        while ((!found) && (ipos < hist_size))
        {
            if ((SIG5(red[i]) == histogram[ipos].red) &&
                (SIG5(green[i]) == histogram[ipos].green) &&
                (SIG5(blue[i]) == histogram[ipos].blue))
            {
                red[i] = hist2[ipos].red;
                green[i] = hist2[ipos].green;
                blue[i] = hist2[ipos].blue;
                found = true;
            }
            ipos++;
        }
        pcent = (int)((i * 50) / size);
        if (pcent > old)
        {
            printf(".");
            old = pcent;
        }
    }
    // Calculate the quantization error
    printf("Ok\nCalculating error");
    pcent = old = 0;

    quantization_error = 0;

    float dist;
    Hist H;
    old = 0;
    for (ipos = 0; ipos < (image_size); ipos++)
    {
        found = false;
        ipont = 0;
        while ((!found) && (ipont < hist_size))
        {
            if ((SIG5(red[ipos]) == histogram[ipont].red) &&
                (SIG5(green[ipos]) == histogram[ipont].green) &&
                (SIG5(blue[ipos]) == histogram[ipont].blue))
            {
                H.red = red[ipos];
                H.green = green[ipos];
                H.blue = blue[ipos];
                dist = (float)(euclidian_distance(H, hist2[ipont]));
                // histogram[ipont].erro += (float) (dist / histograma[ipont].freq);
                quantization_error += dist;
                found = true;
            }
            ipont++;
        }
        pcent = (int)((ipos * 50) / image_size);
        if (pcent > old)
        {
            printf(".");
            old = pcent;
        }
    }

    /*	quantization_error = 0;
	for (i=0; i<hist_size; i++){
		quantization_error += histograma[i].erro;
	}
    */
    quantization_error /= image_size;
    printf("%f\n", quantization_error);

    free(histogram);
}

/******************************************************************
 * Quantization algorithm by Octree
 * Inputs:
 *   int width         - image width
 *   int height        - image height
 *   byte red          - vector for red color
 *   byte gree         - vector for green color
 *   byte blue         - vector for blue color
 *   u_long newpalsize - colors quantity of new image
 * ****************************************************************/

void octree(int width, int height, byte *red, byte *green, byte *blue,
            u_long newpalsize)
{
    u_long image_size, index;
    int pcent, old;

    image_size = width * height;

    printf("Quantization by Octree\n");

    // Calculates the colors number on image and create histogram
    create_histgram(image_size, red, green, blue);

    // If colors size in new image is bigger or equals old, exit
    if (hist_size <= newpalsize)
    {
        printf("It is not necessary to quantize\n");
        exit(1);
    }

    /* Initializing Octree */
    inicialize_octree(newpalsize);

    // Determination of representative colors
    printf("Inserting");
    pcent = old = 0;

    for (index = 0; index < hist_size; index++)
    {
        // If the octree is already complete, make a reduction on the octree.
        track_octree(histogram[index].red, histogram[index].green, histogram[index].blue);
        if (cubo.idx_colors >= cubo.maximum_colors)
            reduces_octree();

        pcent = (int)((index * 50) / hist_size);
        if (pcent > old)
        {
            printf(".");
            old = pcent;
        }
    }
    printf("Ok\n");

    // Mapping from original colors to representative colors
    Node *pont;
    byte i;

    printf("Mapping");
    pcent = old = 0;

    bool found;
    byte sr, sg, sb;
    Hist_Float h1, h2;

    quantization_error = 0;

    // For each color of image
    for (index = 0; index < image_size; index++)
    {
        // Start search on root node
        pont = cubo.root;
        // Look on octree for your occurrence
        i = 0;
        while (i < 8)
        {
            found = false;
            if (pont->child[i] != NULL)
            {
                sr = SIG5(red[index]);
                sg = SIG5(green[index]);
                sb = SIG5(blue[index]);
                if ((sr >= pont->child[i]->er) && (sr <= pont->child[i]->dr) &&
                    (sg >= pont->child[i]->eg) && (sg <= pont->child[i]->dg) &&
                    (sb >= pont->child[i]->eb) && (sb <= pont->child[i]->db))
                {
                    pont = pont->child[i];
                    i = 0;
                    found = true;
                }
            }
            if (!found)
                i++;
        }

        // Manipulation to supply error
        // If not found, go to first to looking for it
        i = 0;
        while (i < 8)
        {
            if (pont->child[i] != NULL)
            {
                pont = pont->child[i];
                i = 0;
            }
            else
                i++;
        }

        h1.red = red[index];
        h1.green = green[index];
        h1.blue = blue[index];
        h2.red = pont->r;
        h2.green = pont->g;
        h2.blue = pont->b;
        quantization_error += euclidian_distance(h1, h2);

        // Found its cell of quantity and replaces original color by quantized
        red[index] = (byte)pont->r;
        green[index] = (byte)pont->g;
        blue[index] = (byte)pont->b;

        pcent = (int)((index * 50) / image_size);
        if (pcent > old)
        {
            printf(".");
            old = pcent;
        }
    }
    printf("Ok\n");
    quantization_error /= image_size;
    printf("Quantization error : %f\n", quantization_error);

    free(histogram);
}

// Initialization and dynamic allocation for octree variables
void inicialize_octree(u_long qt_cores)
{
    byte i;
    u_long j;

    cubo.root = (Node *)malloc(sizeof(Node));
    cubo.root->er = cubo.root->eg = cubo.root->eb = 0;
    cubo.root->dr = cubo.root->dg = cubo.root->db = 255;
    cubo.root->level = 0;
    cubo.root->father = NULL;
    cubo.root->r = 0;
    cubo.root->g = 0;
    cubo.root->b = 0;
    for (i = 0; i < 8; i++)
        cubo.root->child[i] = NULL;
    cubo.idx_colors = 0;
    cubo.maximum_colors = qt_cores;
    cubo.list = (Node **)malloc(sizeof(Node *) * qt_cores + 2);
    for (j = 0; j <= qt_cores + 2; j++)
        cubo.list[j] = NULL;
}

// Defines area of subcube using datas of RGB colors by father subcube
void subcube_area(Node *pont)
{
    byte meio;

    // Red color area of subcube
    meio = (byte)(((pont->father->dr - pont->father->er) / 2) + pont->father->er);
    if (pont->r > meio)
    {
        pont->er = meio + 1;
        pont->dr = pont->father->dr;
    }
    else
    {
        pont->er = pont->father->er;
        pont->dr = meio;
    }

    // Green color area of subcube
    meio = (byte)(((pont->father->dg - pont->father->eg) / 2) + pont->father->eg);
    if (pont->g > meio)
    {
        pont->eg = meio + 1;
        pont->dg = pont->father->dg;
    }
    else
    {
        pont->eg = pont->father->eg;
        pont->dg = meio;
    }

    // Blue color area of subcube
    meio = (byte)(((pont->father->db - pont->father->eb) / 2) + pont->father->eb);
    if (pont->b > meio)
    {
        pont->eb = meio + 1;
        pont->db = pont->father->db;
    }
    else
    {
        pont->eb = pont->father->eb;
        pont->db = meio;
    }
}

// Allocate, insert and initialize a new node on index position with RGB colors
void insert_node(Node *pont, byte indice, byte red, byte green, byte blue)
{
    byte i;

    pont->child[indice] = (Node *)malloc(sizeof(Node));
    pont->child[indice]->father = pont;
    pont = pont->child[indice];
    pont->r = red;
    pont->g = green;
    pont->b = blue;
    pont->level = pont->father->level + 1;
    pont->n_colors = 1;

    // Adjust subucube area
    subcube_area(pont);
    // Put NULL to all children of new node
    for (i = 0; i < 8; i++)
    {
        pont->child[i] = NULL;
    }
    // insert the new node on leaf list of tree
    cubo.list[cubo.idx_colors] = pont;
    cubo.idx_colors++;
};

// Search a position to insert a RGB color and insert.
void track_octree(byte red, byte green, byte blue)
{
    Node *pont, *pont_aux, new_node;
    bool found, son_exists, entred;
    byte i;

    pont = cubo.root;
    son_exists = false;
    entred = false;
    i = 0;

    // Go down to greatest depth
    for (i = 0; i < 8; i++)
    {
        if (entred)
        {
            i = 0;
            entred = false;
        }
        if (pont->child[i] != NULL)
        {
            son_exists = true;
            if ((red >= pont->child[i]->er) && (red <= pont->child[i]->dr) &&
                (green >= pont->child[i]->eg) && (green <= pont->child[i]->dg) &&
                (blue >= pont->child[i]->eb) && (blue <= pont->child[i]->db))
            {
                pont = pont->child[i];
                i = 0;
                son_exists = false;
                entred = true;
            }
        }
    }
    found = ((!son_exists) && (pont != cubo.root));

    // If not found subcube of entred color...
    // It is the first case described on monography
    if (!found)
    {
        i = 0;
        // Find a void child of node to insert a new node
        while (pont->child[i] != NULL)
        {
            i++;
        }
        insert_node(pont, i, red, green, blue);
    }
    else
    {
        // If the color exists on octree, just increment color counter
        // It is second case described on monography
        if ((pont->r == red) && (pont->g == green) && (pont->b == blue))
        {
            pont->n_colors++;
        }
        // Else, if exists the cube, but not exists the color, create
        // a cube to existing color and insert other subcube to the color
        // to insert then.
        // It is third case described on monography
        else
        {
            // While colors are not on distinct subcubes
            do
            {
                // Create a new node and insert the new node
                pont_aux = (Node *)malloc(sizeof(Node));
                // Set NULL to all chindren of new node
                for (i = 0; i < 8; i++)
                {
                    pont_aux->child[i] = NULL;
                }
                // Assumes the area of existing cube
                pont_aux->er = pont->er;
                pont_aux->dr = pont->dr;
                pont_aux->eg = pont->eg;
                pont_aux->dg = pont->dg;
                pont_aux->eb = pont->eb;
                pont_aux->db = pont->db;
                // Takes father as the existin node fater and turns father of existing node
                pont_aux->father = pont->father;
                pont->father = pont_aux;
                pont_aux->child[0] = pont;
                i = 0;
                while (pont_aux->father->child[i] != pont)
                {
                    i++;
                }
                pont_aux->father->child[i] = pont_aux;

                pont_aux->level = pont->level;
                pont->level++;
                // adjust the new area of subcube of chid node
                subcube_area(pont);

                // Insert a brother
                new_node.father = pont_aux;
                new_node.r = red;
                new_node.g = green;
                new_node.b = blue;
                subcube_area(&new_node);
            } while ((pont->er == new_node.er) && (pont->dr == new_node.dr) &&
                     (pont->eg == new_node.eg) && (pont->dg == new_node.dg) &&
                     (pont->eb == new_node.eb) && (pont->db == new_node.db));
            insert_node(pont_aux, 1, red, green, blue);
        }
    }
}

// Reduces the octree in case the insertion of colors exceeds
// the maximum number of colors.
void reduces_octree(void)
{
    u_long i, i_aux, smaller_n_colors, smaller_aux, bigger_prof;
    byte j, n_children;
    Node *node_reducible, *node_aux;
    bool exists;

    // Initialize with first element of list with more than one child
    exists = false;
    i = 0;
    while ((i < cubo.idx_colors) && (!exists))
    {
        // If node has brothers, its father can be reduced
        n_children = 0;
        for (j = 0; j < 8; j++)
        {
            if (cubo.list[i]->father->child[j] != NULL)
                n_children++;
        }
        if (n_children > 1)
        {
            node_reducible = cubo.list[i]->father;
            bigger_prof = cubo.list[i]->level;
            smaller_n_colors = 0;
            for (j = 0; j < 8; j++)
            {
                if (node_reducible->child[j] != NULL)
                {
                    smaller_n_colors += node_reducible->child[j]->n_colors;
                }
            }
            exists = true;
        }
        i++;
    }

    if (!exists)
    {
        exit(1);
    }
    // Find node with greater depth and smaller colors number until now
    int aux = i;
    for (i = aux; i < cubo.idx_colors; i++)
    {
        // If node level is the greater depth until now
        if (cubo.list[i]->level > bigger_prof)
        {
            // If node has brothers, than your father can be reduced
            n_children = 0;
            for (j = 0; j < 8; j++)
            {
                if (cubo.list[i]->father->child[j] != NULL)
                    n_children++;
            }
            if (n_children > 1)
            {
                node_reducible = cubo.list[i]->father;
                bigger_prof = cubo.list[i]->level;
                smaller_n_colors = 0;
                // Sum the brothers frequence
                for (j = 0; j < 8; j++)
                {
                    if (node_reducible->child[j] != NULL)
                    {
                        smaller_n_colors = node_reducible->child[j]->n_colors;
                    }
                }
            }
        }
        else
        {
            // Else if it is equals and has less colors until now
            if (cubo.list[i]->level == bigger_prof)
            {
                // If it is not a brother
                if (node_reducible != cubo.list[i]->father)
                {
                    smaller_aux = 0;
                    // Sum the frequence of all brothers
                    n_children = 0;
                    for (j = 0; j < 8; j++)
                    {
                        if (cubo.list[i]->father->child[j] != NULL)
                        {
                            smaller_aux += cubo.list[i]->father->child[j]->n_colors;
                            n_children++;
                        }
                    }
                    // If has less colors and has brothers, its father can be reduced
                    if ((smaller_n_colors > smaller_aux) && (n_children > 1))
                    {
                        node_reducible = cubo.list[i]->father;
                        bigger_prof = cubo.list[i]->level;
                        smaller_n_colors = smaller_aux;
                    }
                }
            }
        }
    }

    // Reduces the founded node
    float cR, cG, cB;
    cR = cG = cB = 0;
    node_reducible->n_colors = 0;

    for (j = 0; j < 8; j++)
    {
        if (node_reducible->child[j] != NULL)
        {
            cR += node_reducible->child[j]->r * node_reducible->child[j]->n_colors;
            cG += node_reducible->child[j]->g * node_reducible->child[j]->n_colors;
            cB += node_reducible->child[j]->b * node_reducible->child[j]->n_colors;
            node_reducible->n_colors += node_reducible->child[j]->n_colors;
        }
    }
    node_reducible->r = (float)(cR / node_reducible->n_colors);
    node_reducible->g = (float)(cG / node_reducible->n_colors);
    node_reducible->b = (float)(cB / node_reducible->n_colors);

    // Remove children and delete from the list
    bool found;
    byte reduced_size = 0;

    for (j = 0; j < 8; j++)
    {
        if (node_reducible->child[j] != NULL)
        {
            reduced_size++;
            found = false;
            i = 0;
            // Remove the node child from the list
            while ((!found) && (i < cubo.idx_colors))
            {
                found = (cubo.list[i] == node_reducible->child[j]);
                if (found)
                {
                    cubo.list[i] = NULL;
                }
                i++; // Next element of the list
            }
            free(node_reducible->child[j]);
            node_reducible->child[j] = NULL;
        }
    }

    // removes empty spaces left in the list by removing the children of the reducible node
    i = 0;
    while (i < cubo.idx_colors)
    {
        if (cubo.list[i] == NULL)
        {
            // Cleaning void spaces
            i_aux = i;
            while (i_aux < cubo.idx_colors - 1)
            {
                cubo.list[i_aux] = cubo.list[i_aux + 1];
                cubo.list[i_aux + 1] = NULL;
                i_aux++;
            }
            cubo.idx_colors--;
        }
        else
        {
            i++; // Next element on the list
        }
    }

    // Insert a new node reduced on the leaf list
    cubo.list[cubo.idx_colors] = node_reducible;
    cubo.idx_colors++;

    // Reduces the levels until the node that has brothes and is son of the root
    n_children = 0;
    while ((n_children < 2) && (node_reducible->father != cubo.root))
    {
        for (j = 0; j < 8; j++)
        {
            if (node_reducible->father->child[j] != NULL)
                n_children++;
        }
        if (n_children < 2)
        {
            i = 0;
            while (node_reducible->father->father->child[i] != node_reducible->father)
                i++;
            node_reducible->father->father->child[i] = node_reducible;
            node_aux = node_reducible->father;
            node_reducible->father = node_reducible->father->father;
            node_reducible->level--;
            subcube_area(node_reducible);
            free(node_aux);
            n_children = 0;
        }
    }
}
