#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES

#include <Windows.h>
#include <iostream>
#include <GL/glew.h>
#include <GL/GL.h>
#include <GL/freeglut.h>

#define GLFW_INCLUDE_GLU
#define GLFW_DLL
#include <GLFW/glfw3.h>
#include <vector>

#define GLM_SWIZZLE
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>


#include <cstdio>
#include <cmath>
#include <algorithm>

using namespace glm;

// -------------------------------------------------
// Global Variables
// -------------------------------------------------
int Width = 512;
int Height = 512;
std::vector<float> OutputImage;
// -------------------------------------------------



struct Vec3 {
    float x, y, z;
};

struct Vec2 {
    float x, y;
};

struct Color {
    unsigned char r, g, b;
};

int     gNumVertices = 0;    // Number of 3D vertices.
int     gNumTriangles = 0;    // Number of triangles.
int* gIndexBuffer = NULL; // Vertex indices for the triangles.
Vec3* gVertexBuffer = NULL; // 3D Vertex postion (x, y, z)
Vec3* gNormalBuffer = NULL;

const int SCREEN_WIDTH = 512;
const int SCREEN_HEIGHT = 512;
Color framebuffer[SCREEN_HEIGHT][SCREEN_WIDTH]; // 출력 버퍼
float depthbuffer[SCREEN_HEIGHT][SCREEN_WIDTH]; // Z-buffer

const Vec3 ka = { 0.0f, 1.0f, 0.0f };
const Vec3 kd = { 0.0f, 0.5f, 0.0f };
const Vec3 ks = { 0.5f, 0.5f, 0.5f };
const float p = 32.0f;
const float Ia = 0.2f;
const Vec3 lightPos = { -4.0f, 4.0f, -3.0f };
const Vec3 lightIntensity = { 1.0f, 1.0f, 1.0f };
const float gamma = 2.2f;


Vec3 operator+(const Vec3& a, const Vec3& b) {
    return { a.x + b.x, a.y + b.y, a.z + b.z };
}

Vec3 operator-(const Vec3& a, const Vec3& b) {
    return { a.x - b.x, a.y - b.y, a.z - b.z };
}

Vec3 operator*(const Vec3& a, const Vec3& b) {
    return { a.x * b.x, a.y * b.y, a.z * b.z };
}

Vec3 operator*(const Vec3& a, float scalar) {
    return { a.x * scalar, a.y * scalar, a.z * scalar };
}

Vec3 operator*(float scalar, const Vec3& a) {
    return a * scalar;
}

float dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 normalize(Vec3 v) {
    float len = std::sqrt(dot(v, v));
    return { v.x / len, v.y / len, v.z / len };
}


void create_scene()
{
    int width = 32;
    int height = 16;

    float theta, phi;
    int t;

    gNumVertices = (height - 2) * width + 2;
    gNumTriangles = (height - 2) * (width - 1) * 2;

    // TODO: Allocate an array for gNumVertices vertices.

    gIndexBuffer = new int[3 * gNumTriangles];
    gVertexBuffer = new Vec3[gNumVertices];

    t = 0;
    for (int j = 1; j < height - 1; ++j)
    {
        for (int i = 0; i < width; ++i)
        {
            theta = (float)j / (height - 1) * M_PI;
            phi = (float)i / (width - 1) * M_PI * 2;

            float   x = sinf(theta) * cosf(phi);
            float   y = cosf(theta);
            float   z = -sinf(theta) * sinf(phi);

            // TODO: Set vertex t in the vertex array to {x, y, z}.

            gVertexBuffer[t] = { x, y, z };
            t++;
        }
    }

    // TODO: Set vertex t in the vertex array to {0, 1, 0}.

    gVertexBuffer[t] = { 0, 1, 0 };
    t++;

    // TODO: Set vertex t in the vertex array to {0, -1, 0}.

    gVertexBuffer[t] = { 0, -1, 0 };
    t++;

    t = 0;
    for (int j = 0; j < height - 3; ++j)
    {
        for (int i = 0; i < width - 1; ++i)
        {
            gIndexBuffer[t++] = j * width + i;
            gIndexBuffer[t++] = (j + 1) * width + (i + 1);
            gIndexBuffer[t++] = j * width + (i + 1);
            gIndexBuffer[t++] = j * width + i;
            gIndexBuffer[t++] = (j + 1) * width + i;
            gIndexBuffer[t++] = (j + 1) * width + (i + 1);
        }
    }
    for (int i = 0; i < width - 1; ++i)
    {
        gIndexBuffer[t++] = (height - 2) * width;
        gIndexBuffer[t++] = i;
        gIndexBuffer[t++] = i + 1;
        gIndexBuffer[t++] = (height - 2) * width + 1;
        gIndexBuffer[t++] = (height - 3) * width + (i + 1);
        gIndexBuffer[t++] = (height - 3) * width + i;
    }

    gNormalBuffer = new Vec3[gNumVertices];

    for (int i = 0; i < gNumVertices; ++i) {
        gNormalBuffer[i] = normalize(gVertexBuffer[i]);
    }

    // The index buffer has now been generated. Here's how to use to determine
    // the vertices of a triangle. Suppose you want to determine the vertices
    // of triangle i, with 0 <= i < gNumTriangles. Define:
    //
    // k0 = gIndexBuffer[3*i + 0]
    // k1 = gIndexBuffer[3*i + 1]
    // k2 = gIndexBuffer[3*i + 2]
    //
    // Now, the vertices of triangle i are at positions k0, k1, and k2 (in that
    // order) in the vertex array (which you should allocate yourself at line
    // 27).
    //
    // Note that this assumes 0-based indexing of arrays (as used in C/C++,
    // Java, etc.) If your language uses 1-based indexing, you will have to
    // add 1 to k0, k1, and k2.
}

Vec3 modeling_transform(Vec3 v)
{
    v = v * 2;
    v.z -= 7;

    return v;
}

Vec3 camera_transform(Vec3 v) {
    const Vec3 e = { 0.0f, 0.0f, 0.0f };
    const Vec3 u = { 1.0f, 0.0f, 0.0f };
    const Vec3 v_dir = { 0.0f, 1.0f, 0.0f };
    const Vec3 w = { 0.0f, 0.0f, 1.0f };

    Vec3 eye_to_point = { v.x - e.x, v.y - e.y, v.z - e.z };

    float x = dot(eye_to_point, u);
    float y = dot(eye_to_point, v_dir);
    float z = -dot(eye_to_point, w);

    return { x, y, z };
}

Vec3 perspective_projection(Vec3 v)
{
    float l = -0.1f, r = 0.1f;
    float b = -0.1f, t = 0.1f;
    float n = -0.1f, f = -1000.0f;

    float A = (2 * n) / (r - l);
    float B = (2 * n) / (t - b);
    float C = (r + l) / (r - l);
    float D = (t + b) / (t - b);
    float E = (f + n) / (n - f);
    float F = (2 * f * n) / (n - f);

    float x_proj = A * v.x + C * v.z;
    float y_proj = B * v.y + D * v.z;
    float z_proj = E * v.z + F;
    float w_proj = -v.z;

    return {
        x_proj / w_proj,
        y_proj / w_proj,
        z_proj / w_proj
    };
}

Vec2 viewport_transform(Vec3 v)
{
    float x_prime = (v.x + 1.0f) * 0.5f * SCREEN_WIDTH;
    float y_prime = (1.0f - (1.0f - v.y) * 0.5f) * SCREEN_HEIGHT;

    return { x_prime, y_prime };
}

Vec2 world_to_screen(Vec3 v)
{
    v = modeling_transform(v);
    v = camera_transform(v);
    v = perspective_projection(v);
    return viewport_transform(v);
}

bool inside_triangle(Vec2 p, Vec2 a, Vec2 b, Vec2 c, float& alpha, float& beta, float& gamma)
{
    float denomBeta = (a.y - c.y) * b.x + (c.x - a.x) * b.y + a.x * c.y - c.x * a.y;
    float denomGamma = (a.y - b.y) * c.x + (b.x - a.x) * c.y + a.x * b.y - b.x * a.y;

    if (denomBeta == 0 || denomGamma == 0)
        return false;

    beta = ((a.y - c.y) * p.x + (c.x - a.x) * p.y + a.x * c.y - c.x * a.y) / denomBeta;
    gamma = ((a.y - b.y) * p.x + (b.x - a.x) * p.y + a.x * b.y - b.x * a.y) / denomGamma;
    alpha = 1.0f - beta - gamma;

    return (alpha >= 0.0f && beta >= 0.0f && gamma >= 0.0f &&
        alpha <= 1.0f && beta <= 1.0f && gamma <= 1.0f);
}

Vec3 phong_shading(Vec3 pos, Vec3 normal)
{
    Vec3 L = normalize(lightPos - pos);
    Vec3 V = normalize(pos);
    Vec3 H = normalize(L + V);

    float NL = std::max(0.0f, dot(normal, L));
    float NH = std::max(0.0f, dot(normal, H));

    Vec3 ambient = ka * Ia;
    Vec3 diffuse = kd * lightIntensity * NL;
    Vec3 specular = ks * lightIntensity * std::pow(NH, p);
    Vec3 linearRGB = ambient + diffuse + specular;

    Vec3 corrected = {
        std::pow(std::min(linearRGB.x, 1.0f), 1.0f / gamma),
        std::pow(std::min(linearRGB.y, 1.0f), 1.0f / gamma),
        std::pow(std::min(linearRGB.z, 1.0f), 1.0f / gamma)
    };

    return corrected;
}

void draw_triangle(Vec3 v0, Vec3 v1, Vec3 v2, Vec3 n0, Vec3 n1, Vec3 n2)
{
    Vec2 a = world_to_screen(v0);
    Vec2 b = world_to_screen(v1);
    Vec2 c = world_to_screen(v2);

    int minX = std::max(0, (int)floorf(std::min({ a.x, b.x, c.x })));
    int maxX = std::min(SCREEN_WIDTH - 1, (int)ceilf(std::max({ a.x, b.x, c.x })));
    int minY = std::max(0, (int)floorf(std::min({ a.y, b.y, c.y })));
    int maxY = std::min(SCREEN_HEIGHT - 1, (int)ceilf(std::max({ a.y, b.y, c.y })));

    for (int y = minY; y <= maxY; ++y) {
        for (int x = minX; x <= maxX; ++x) {
            float alpha, beta, gamma;

            if (inside_triangle({ (float)x, (float)y }, a, b, c, alpha, beta, gamma)) {
                float z_ndc = alpha * v0.z + beta * v1.z + gamma * v2.z;
                float z_screen = z_ndc * 0.5f + 0.5f;

                if (z_screen < depthbuffer[y][x]) {
                    Vec3 pos = v0 * alpha + v1 * beta + v2 * gamma;
                    Vec3 normal = normalize(n0 * alpha + n1 * beta + n2 * gamma);
                    Vec3 color = phong_shading(pos, normal);

                    framebuffer[y][x] = {
                        (unsigned char)(std::min(color.x * 255.0f, 255.0f)),
                        (unsigned char)(std::min(color.y * 255.0f, 255.0f)),
                        (unsigned char)(std::min(color.z * 255.0f, 255.0f))
                    };

                    depthbuffer[y][x] = z_screen;
                }
            }
        }
    }
}

void render_scene()
{
    for (int y = 0; y < SCREEN_HEIGHT; ++y)
    {
        for (int x = 0; x < SCREEN_WIDTH; ++x)
        {
            framebuffer[y][x] = { 0, 0, 0 };
            depthbuffer[y][x] = 1e9;
        }
    }

    for (int i = 0; i < gNumTriangles; ++i)
    {
        int k0 = gIndexBuffer[3 * i + 0];
        int k1 = gIndexBuffer[3 * i + 1];
        int k2 = gIndexBuffer[3 * i + 2];

        draw_triangle(
            gVertexBuffer[k0],
            gVertexBuffer[k1],
            gVertexBuffer[k2],
            gNormalBuffer[k0],
            gNormalBuffer[k1],
            gNormalBuffer[k2]
        );
    }
}

void delete_memory()
{
    delete gIndexBuffer;
    delete gVertexBuffer;
    delete gNormalBuffer;
}




void render()
{
    OutputImage.clear();
    for (int j = 0; j < Height; ++j)
    {
        for (int i = 0; i < Width; ++i)
        {
            Color color = framebuffer[j][i];

            OutputImage.push_back(color.r / 255.0f); // R
            OutputImage.push_back(color.g / 255.0f); // G
            OutputImage.push_back(color.b / 255.0f); // B
        }
    }
}


void resize_callback(GLFWwindow*, int nw, int nh) 
{
	//This is called in response to the window resizing.
	//The new width and height are passed in so we make 
	//any necessary changes:
	Width = nw;
	Height = nh;
	//Tell the viewport to use all of our screen estate
	glViewport(0, 0, nw, nh);

	//This is not necessary, we're just working in 2d so
	//why not let our spaces reflect it?
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	glOrtho(0.0, static_cast<double>(Width)
		, 0.0, static_cast<double>(Height)
		, 1.0, -1.0);

	//Reserve memory for our render so that we don't do 
	//excessive allocations and render the image
	OutputImage.reserve(Width * Height * 3);
	render();
}


int main(int argc, char* argv[])
{
	// -------------------------------------------------
	// Initialize Window
	// -------------------------------------------------

    create_scene();
    render_scene();

	GLFWwindow* window;

	/* Initialize the library */
	if (!glfwInit())
		return -1;

	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(Width, Height, "OpenGL Viewer", NULL, NULL);
	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	/* Make the window's context current */
	glfwMakeContextCurrent(window);

	//We have an opengl context now. Everything from here on out 
	//is just managing our window or opengl directly.

	//Tell the opengl state machine we don't want it to make 
	//any assumptions about how pixels are aligned in memory 
	//during transfers between host and device (like glDrawPixels(...) )
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glPixelStorei(GL_PACK_ALIGNMENT, 1);

	//We call our resize function once to set everything up initially
	//after registering it as a callback with glfw
	glfwSetFramebufferSizeCallback(window, resize_callback);
	resize_callback(NULL, Width, Height);

	/* Loop until the user closes the window */
	while (!glfwWindowShouldClose(window))
	{
		//Clear the screen
		glClear(GL_COLOR_BUFFER_BIT);

		// -------------------------------------------------------------
		//Rendering begins!
		glDrawPixels(Width, Height, GL_RGB, GL_FLOAT, &OutputImage[0]);
		//and ends.
		// -------------------------------------------------------------

		/* Swap front and back buffers */
		glfwSwapBuffers(window);

		/* Poll for and process events */
		glfwPollEvents();

		//Close when the user hits 'q' or escape
		if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS
			|| glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
		{
			glfwSetWindowShouldClose(window, GL_TRUE);
		}
	}

	glfwDestroyWindow(window);
	glfwTerminate();


    delete_memory();


	return 0;
}
