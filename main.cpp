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

double camera_speed;
double yaw;
double yaw_step;

struct point {
    double x, y, z;
};

class Vector
{
public:
    double x, y, z;

    Vector() {}
    Vector(double, double, double);
    Vector operator + (Vector const&);
    Vector operator - (Vector const&);
    Vector cross(Vector const&);
    Vector normalize();
};

Vector::Vector(double comp_x, double comp_y, double comp_z)
{
    x = comp_x;
    y = comp_y;
    z = comp_z;
}

Vector Vector::operator + (Vector const& vec)
{
    Vector res;
    res.x = x + vec.x;
    res.y = y + vec.y;
    res.z = z + vec.z;
    return res;
}

Vector Vector::operator - (Vector const& vec)
{
    Vector res;
    res.x = x - vec.x;
    res.y = y - vec.y;
    res.z = z - vec.z;
    return res;
}

Vector Vector::cross(Vector const& vec)
{
    Vector res;
    res.x = y * vec.z - z * vec.y;
    res.y = z * vec.x - x * vec.z;
    res.z = x * vec.y - y * vec.x;
    return res;
}

Vector Vector::normalize()
{
    double mag = sqrt(x * x + y * y + z * z);

    Vector res;
    res.x = x / mag;
    res.y = y / mag;
    res.z = z / mag;
    return res;
}

Vector operator*(double scale, Vector const& vec)
{
    Vector res;
    res.x = scale * vec.x;
    res.y = scale * vec.y;
    res.z = scale * vec.z;
    return res;
}

Vector camera_pos, direc, camera_up;

void initGL()
{   
    cube_size = 50.0;
    sphere_radius = 10.0;
    cylinder_len = cube_size - 2*sphere_radius;

    camera_speed = 2.0;
    yaw = 0;
    yaw_step = 5;

    camera_pos.x = 0;
    camera_pos.y = 0;
    camera_pos.z = 100.0;
    
    direc.x = 0;
    direc.y = 0;
    direc.z = -100.0;

    camera_up.x = 0;
    camera_up.y = 1.0;
    camera_up.z = 0;

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

void drawCylinderOneFourth(double radius,double len,int slices)
{
    std::vector<point> points;

    for (int i = 0; i <= slices; i++)
    {
        point p;
        p.x = radius * cos(((double)i / (double)slices) * (pi / 2));
        p.y = radius * sin(((double)i / (double)slices) * (pi / 2));
        p.z = 0;
        points.push_back(p);
    }

    glColor3f(0.0f, 1.0f, 0.0f);
    for (int i = 0; i < slices; i++)
    {
        glBegin(GL_QUADS);
        {
            glVertex3f(points[i].x, points[i].y, points[i].z);
            glVertex3f(points[i].x, points[i].y, points[i].z-len);
            glVertex3f(points[i+1].x, points[i+1].y, points[i+1].z-len);
            glVertex3f(points[i+1].x, points[i+1].y, points[i+1].z);
        }
        glEnd();
    }
}

void drawSquare(double len)
{   
    glColor3f(1.0f, 1.0f, 1.0f);
    glBegin(GL_QUADS);
    {
        glVertex3f(-len / 2, -len / 2, 0);
        glVertex3f(len / 2, -len / 2, 0);
        glVertex3f(len/2,len/2,0);
        glVertex3f(-len / 2, len / 2, 0);
    }
    glEnd();
}

void drawCubeFaces()
{
    glPushMatrix();
    glTranslatef(0, 0, cube_size/2);
    drawSquare(cylinder_len);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(cube_size/2, 0, 0);
    glRotatef(90.0, 0, 1.0, 0);
    drawSquare(cylinder_len);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0, 0, -cube_size/2);
    drawSquare(cylinder_len);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-cube_size/2, 0, 0);
    glRotatef(90.0, 0, 1.0, 0);
    drawSquare(cylinder_len);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0, cube_size/2, 0);
    glRotatef(90.0, 1.0, 0, 0);
    drawSquare(cylinder_len);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0, -cube_size/2, 0);
    glRotatef(90.0, 1.0, 0, 0);
    drawSquare(cylinder_len);
    glPopMatrix();
}

void drawCornerSpheres()
{
    glPushMatrix();
    glTranslatef(cylinder_len / 2, cylinder_len / 2, cylinder_len / 2);
    drawSphereOneEighth(sphere_radius, 100, 100);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(cylinder_len / 2, -cylinder_len / 2, cylinder_len / 2);
    glScalef(1.0f, -1.0f, 1.0f);
    drawSphereOneEighth(sphere_radius, 100, 100);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-cylinder_len / 2, cylinder_len / 2, cylinder_len / 2);
    glScalef(-1.0f, 1.0f, 1.0f);
    drawSphereOneEighth(sphere_radius, 100, 100);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-cylinder_len / 2, -cylinder_len / 2, cylinder_len / 2);
    glScalef(-1.0f, -1.0f, 1.0f);
    drawSphereOneEighth(sphere_radius, 100, 100);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(cylinder_len / 2, cylinder_len / 2, -cylinder_len / 2);
    glScalef(1.0f, 1.0f, -1.0f);
    drawSphereOneEighth(sphere_radius, 100, 100);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(cylinder_len / 2, -cylinder_len / 2, -cylinder_len / 2);
    glScalef(1.0f, -1.0f, -1.0f);
    drawSphereOneEighth(sphere_radius, 100, 100);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-cylinder_len / 2, cylinder_len / 2, -cylinder_len / 2);
    glScalef(-1.0f, 1.0f, -1.0f);
    drawSphereOneEighth(sphere_radius, 100, 100);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-cylinder_len / 2, -cylinder_len / 2, -cylinder_len / 2);
    glScalef(-1.0f, -1.0f, -1.0f);
    drawSphereOneEighth(sphere_radius, 100, 100);
    glPopMatrix();
}

void drawCylinderSides()
{   
    glPushMatrix();
    glTranslatef(cylinder_len / 2, cylinder_len / 2, cylinder_len / 2);
    drawCylinderOneFourth(sphere_radius, cylinder_len, 100);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(cylinder_len / 2, -cylinder_len / 2, cylinder_len / 2);
    glRotatef(90.0, 1.0, 0, 0);
    drawCylinderOneFourth(sphere_radius, cylinder_len, 100);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(cylinder_len / 2, -cylinder_len / 2, -cylinder_len / 2);
    glRotatef(-180.0, 1.0, 0, 0);
    drawCylinderOneFourth(sphere_radius, cylinder_len, 100);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(cylinder_len / 2, cylinder_len / 2, -cylinder_len / 2);
    glRotatef(-90.0, 1.0, 0, 0);
    drawCylinderOneFourth(sphere_radius, cylinder_len, 100);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-cylinder_len / 2, cylinder_len / 2, cylinder_len / 2);
    glRotatef(-90.0, 0, 1.0, 0);
    drawCylinderOneFourth(sphere_radius, cylinder_len, 100);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-cylinder_len / 2, -cylinder_len / 2, cylinder_len / 2);
    glRotatef(-90.0, 0, 1.0, 0);
    glRotatef(90.0, 1.0, 0, 0);
    drawCylinderOneFourth(sphere_radius, cylinder_len, 100);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(cylinder_len / 2, -cylinder_len / 2, cylinder_len / 2);
    glRotatef(-90.0, 0, 1.0, 0);
    glRotatef(-180.0, 1.0, 0, 0);
    drawCylinderOneFourth(sphere_radius, cylinder_len, 100);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-cylinder_len / 2, cylinder_len / 2, -cylinder_len / 2);
    glRotatef(-180.0, 0, 1.0, 0);
    drawCylinderOneFourth(sphere_radius, cylinder_len, 100);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-cylinder_len / 2, -cylinder_len / 2, -cylinder_len / 2);
    glRotatef(-180.0, 0, 1.0, 0);
    glRotatef(90.0, 1.0, 0, 0);
    drawCylinderOneFourth(sphere_radius, cylinder_len, 100);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-cylinder_len / 2, -cylinder_len / 2, cylinder_len / 2);
    glRotatef(180.0, 0, 0, 1.0);
    drawCylinderOneFourth(sphere_radius, cylinder_len, 100);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(cylinder_len / 2, cylinder_len / 2, -cylinder_len / 2);
    glRotatef(90.0, 0, 1.0, 0);
    drawCylinderOneFourth(sphere_radius, cylinder_len, 100);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-cylinder_len / 2, -cylinder_len / 2, -cylinder_len / 2);
    glRotatef(90.0, 0, 1.0, 0);
    glRotatef(-180.0, 1.0, 0, 0);
    drawCylinderOneFourth(sphere_radius, cylinder_len, 100);
    glPopMatrix();
}

void display(void) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    gluLookAt(camera_pos.x,camera_pos.y,camera_pos.z, 
        camera_pos.x+direc.x,camera_pos.y+direc.y,camera_pos.z+direc.z,
        camera_up.x,camera_up.y,camera_up.z);

    drawAxes();
    drawCubeFaces();
    drawCornerSpheres();
    drawCylinderSides();

    glutSwapBuffers();
}

void keyboardListener(unsigned char key, int x, int y) {

    double yaw_radian;

    switch (key) {

    case '1':
        break;
    case '2':
        break;
    case '3':
        break;
    case '4':
        break;
    case '5':
        break;
    case '6':
        break;
    default:
        break;
    }
}

void specialKeyListener(int key, int x, int y) {
    Vector right;

    switch (key) {
    case GLUT_KEY_DOWN:		
        camera_pos = camera_pos - camera_speed * direc.normalize();
        break;
    case GLUT_KEY_UP:		
        camera_pos = camera_pos + camera_speed * direc.normalize();
        break;

    case GLUT_KEY_RIGHT:
        right = direc.cross(camera_up).normalize();
        camera_pos = camera_pos + right;
        break;
    case GLUT_KEY_LEFT:
        right = direc.cross(camera_up).normalize();
        camera_pos = camera_pos - right;
        break;

    case GLUT_KEY_PAGE_UP:
        camera_pos = camera_pos + camera_up.normalize();
        break;
    case GLUT_KEY_PAGE_DOWN:
        camera_pos = camera_pos - camera_up.normalize();
        break;

    case GLUT_KEY_HOME:
        if (sphere_radius < (cube_size / 2)) {
            sphere_radius += 1;
            cylinder_len = cube_size - 2 * sphere_radius;
        }
        break;
    case GLUT_KEY_END:
        if (sphere_radius > 0) {
            sphere_radius -= 1;
            cylinder_len = cube_size - 2 * sphere_radius;
        }
        break;

    default:
        break;
    }
}

void mouseListener(int button, int state, int x, int y) {
    switch (button) {
    case GLUT_LEFT_BUTTON:
        break;

    case GLUT_RIGHT_BUTTON:
        break;

    case GLUT_MIDDLE_BUTTON:
        break;

    default:
        break;
    }
}


void reshape(GLsizei width, GLsizei height) {
    if (height == 0) height = 1;                
    GLfloat aspect = (GLfloat)width / (GLfloat)height;

    glViewport(0, 0, width, height);

    glMatrixMode(GL_PROJECTION); 
    glLoadIdentity();

    gluPerspective(80.0f, aspect, 0.1f, 1000.0f);
}

void animate() {
    glutPostRedisplay();
}


int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(640, 640);                    
    glutInitWindowPosition(200, 200);                
    glutCreateWindow("Camera Movement"); 

    initGL();

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutIdleFunc(animate);

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutMouseFunc(mouseListener);

    glutMainLoop();
    return 0;
}