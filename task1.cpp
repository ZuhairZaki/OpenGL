#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<vector>

#include <windows.h>
#include <GL\freeglut.h>

#define pi (2*acos(0.0))

double cube_size;
double cylinder_len;
double sphere_radius;

struct point {
    double x, y, z;
};

/**
* glVertex2f(float x, float y).
* The point (0.0, 0.0) represents the middle of the window (not the top left corner).
* The "2f" suffix means 2 values of float type (x and y).
*/
void initGL()
{
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClearDepth(1.0f);
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);
}

void drawAxes()
{
    glColor3f(1.0f, 1.0f, 1.0f);
    glBegin(GL_LINES); {
        glVertex3f(100.0f, 0.0f, 0.0f);
        glVertex3f(-100.0f, 0.0f, 0.0f);

        glVertex3f(0.0f, -100.0f, 0.0f);
        glVertex3f(0.0f, 100.0f, 0.0f);

        glVertex3f(0.0f, 0.0f, 100.0f);
        glVertex3f(0.0f, 0.0f, -100.0f);
    }glEnd();
}

void drawSphereOneEighth(double radius, int slices, int stacks)
{   
    std::vector<std::vector<point>> points(stacks+1);

    
    for (int i = 0; i <= stacks; i++)
    {
        double z = radius * sin(((double)i / (double)stacks) * (pi / 2));
        double xy_radius = radius * cos(((double)i / (double)stacks) * (pi / 2));
        for (int j = 0; j <= slices; j++)
        {   
            point p;
            p.x = xy_radius * cos(((double)j / (double)slices) * (pi/2));
            p.y = xy_radius * sin(((double)j / (double)slices) * (pi/2));
            p.z = z;
            points[i].push_back(p);
        }
    }
    
    glColor3f(1.0f, 0.0f, 0.0f);
    for (int i = 0; i < stacks; i++)
    {
        for (int j = 0; j < slices; j++)
        {
            glBegin(GL_QUADS); {
                glVertex3f(points[i][j].x, points[i][j].y, points[i][j].z);
                glVertex3f(points[i][j + 1].x, points[i][j + 1].y, points[i][j + 1].z);
                glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, points[i + 1][j + 1].z);
                glVertex3f(points[i + 1][j].x, points[i + 1][j].y, points[i + 1][j].z);
            }glEnd();
        }
    }
}

void display(void) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    gluLookAt(50,50,50,0,0,0,0,1,0);

    glPushMatrix();
    drawAxes();
    glPopMatrix();

    glPushMatrix();
    drawSphereOneEighth(10, 100, 100);
    glPopMatrix();

    glutSwapBuffers();
}

void reshape(GLsizei width, GLsizei height) {
    if (height == 0) height = 1;                
    GLfloat aspect = (GLfloat)width / (GLfloat)height;

    glViewport(0, 0, width, height);

    glMatrixMode(GL_PROJECTION); 
    glLoadIdentity();

    gluPerspective(80.0f, aspect, 0.1f, 1000.0f);
}

int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(640, 640);                    // window size
    glutInitWindowPosition(200, 200);                // distance from the top-left screen
    glutCreateWindow("Camera Movement");    // message displayed on top bar window

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    initGL();

    glutMainLoop();
    return 0;
}