#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include <windows.h>
#include <GL\freeglut.h>

#define pi (2*acos(0.0))

struct point
{
    double x, y, z;
};

class Vector
{
public:
    double x, y, z;

    Vector() {}
    Vector(double, double, double);
    double Magnitude();
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

double Vector::Magnitude()
{
    return sqrt(x * x + y * y + z * z);
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

Vector operator*(Vector const& vec, double scale)
{
    Vector res;
    res.x = vec.x * scale;
    res.y = vec.y * scale;
    res.z = vec.z * scale;
    return res;
}

double wheel_radius, wheel_len;
double turn, rot, rot_step, turn_step;
Vector pos, direc;

Vector camera_pos, camera_up;
double camera_angle, camera_speed;

void initGL()
{   
    wheel_radius = 20.0;
    wheel_len = 10.0;

    turn = 0;
    rot = 0;
    rot_step = 2.0;
    turn_step = 2.0;

    camera_angle = pi / 180;
    camera_speed = 2.0;

    pos = Vector(0, wheel_radius, 0);
    direc = Vector(1, 0, 0);

    camera_pos = Vector(-100,50,100);
    camera_up = Vector(0, 1, 0);


    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClearDepth(1.0f);
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);
}

void drawGrid()
{
    glColor3f(0.6, 0.6, 0.6);
    glBegin(GL_LINES); {
        for (int i = 0; i <= 10; i++) {
            glVertex3f(i * 10, 0, 120);
            glVertex3f(i * 10, 0, -120);

            glVertex3f(-i * 10, 0, 120);
            glVertex3f(-i * 10, 0, -120);

            glVertex3f(120, 0, i*10);
            glVertex3f(-120, 0, i*10);

            glVertex3f(120, 0, -i * 10);
            glVertex3f(-120, 0, -i * 10);
        }
    }glEnd();
}

void drawRectangle(double w,double l)
{   
    glBegin(GL_QUADS);
    {
        glVertex3f(w, 0, l);
        glVertex3f(w, 0, -l);
        glVertex3f(-w, 0, -l);
        glVertex3f(-w, 0, l);
    }
    glEnd();
}

void drawWheel(double radius,double segments, double seg_len)
{   
    glPushMatrix();
    glTranslatef(pos.x, 0, pos.z);
    glRotatef(turn, 0, 1.0, 0);

    double theta = pi / segments;
    double width = radius * tan(theta);

    for (int i = 0; i < segments; i++)
    {   
        double alpha = 2*i*theta - rot*pi/180.0;
        double x = radius * sin(alpha);
        double y = radius * (1 - cos(alpha));

        double shade = 2 * abs(int(segments / 2) - i);
        glColor3f(shade / segments, shade / segments, shade / segments);

        glPushMatrix();
        glTranslatef(x, y,0);
        glRotatef(alpha * 180.0 / pi, 0, 0, 1.0f);
        drawRectangle(width, seg_len/2);
        glPopMatrix();
    }

    glColor3f(0.6, 0.6, 0.6);
    
    glPushMatrix();
    glTranslatef(0,radius, 0);

    glRotatef(-90.0f, 0, 1.0f, 0);
    glRotatef(-rot, 1.0f, 0, 0);
    drawRectangle(seg_len / 4, radius);

    glPushMatrix();
    glRotatef(-90.0f, 1.0f, 0, 0);
    drawRectangle(seg_len / 4, radius);
    glPopMatrix();

    glPopMatrix();

    glPopMatrix();
}

void display(void) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    gluLookAt(camera_pos.x, camera_pos.y,camera_pos.z, 
        0,0, 0, camera_up.x, camera_up.y, camera_up.z);

    drawGrid();
    drawWheel(wheel_radius, 120, wheel_len);

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

void animate() {
    glutPostRedisplay();
}

void keyboardListener(unsigned char key, int x, int y) {
    Vector up(0, 1, 0), right;
    double s;

    switch (key) {

    case 'w':
        if (rot > 360.0)
            rot -= 360.0;

        rot += rot_step;
        s = wheel_radius * rot_step * pi / 180.0;

        pos = pos + s * direc;
        break;
    case 's':
        if (rot < -360.0)
            rot += 360.0;

        rot -= rot_step;
        s = wheel_radius * rot_step * pi / 180.0;

        pos = pos - s * direc;
        break;
    case 'a':
        if (turn > 360.0)
            turn -= 360.0;

        turn += turn_step;
        s = turn_step * pi / 180.0;

        right = direc.cross(up);
        direc = direc * cos(s) - right * sin(s);
        break;
    case 'd':
        if (turn < -360.0)
            turn += 360.0;

        turn -= turn_step;
        s = turn_step * pi / 180.0;

        right = direc.cross(up);
        direc = direc * cos(s) + right * sin(s);
        break;

    default:
        break; 
    }
}

void specialKeyListener(int key, int x, int y) {
    switch (key) {
    case GLUT_KEY_DOWN:
        camera_pos = camera_pos - camera_speed * camera_up;
        break;
    case GLUT_KEY_UP:
        camera_pos = camera_pos + camera_speed * camera_up;
        break;

    case GLUT_KEY_RIGHT:
        camera_pos.x = camera_pos.x * cos(camera_angle) + camera_pos.z * sin(camera_angle);
        camera_pos.z = -camera_pos.x * sin(camera_angle) + camera_pos.z * cos(camera_angle);
        break;
    case GLUT_KEY_LEFT:
        camera_pos.x = camera_pos.x * cos(camera_angle) - camera_pos.z * sin(camera_angle);
        camera_pos.z = camera_pos.x * sin(camera_angle) + camera_pos.z * cos(camera_angle);
        break;

    default:
        break;
    }
}

int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(640, 480);                    
    glutInitWindowPosition(200, 200);               
    glutCreateWindow("Wheel"); 

    initGL();

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutIdleFunc(animate);

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);

    glutMainLoop();
    return 0;
}