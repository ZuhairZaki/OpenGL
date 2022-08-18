#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>

#include <GL\freeglut.h>

#include "classes.hpp"
#include "bitmap_image.hpp"

#define pi (2*acos(0.0))

double camera_speed;
double yaw_step, pitch_step, roll_step;

Vector camera_pos, direc, camera_up;

std::vector<Object*> objects;
std::vector<PointLight> pointLights;
std::vector<SpotLight> spotLights;

void initGL()
{
    camera_speed = 2.0;
    yaw_step = pi / 36;
    pitch_step = pi / 36;
    roll_step = pi / 36;

    camera_pos.x = 100.0;
    camera_pos.y = 0;
    camera_pos.z = 50.0;

    direc.x = -100.0;
    direc.y = 0;
    direc.z = 0;

    camera_up.x = 0;
    camera_up.y = 0;
    camera_up.z = 1.0;
    
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClearDepth(1.0f);
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);
}

void display(void) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    gluLookAt(camera_pos.x, camera_pos.y, camera_pos.z,
        camera_pos.x + direc.x, camera_pos.y + direc.y, camera_pos.z + direc.z,
        camera_up.x, camera_up.y, camera_up.z);

    Floor f(Vector(0, 0, 0), 1000, 20);
    objects.push_back(&f);
    objects[0]->draw();

    glutSwapBuffers();
}

void keyboardListener(unsigned char key, int x, int y) {
    Vector target, right;
    target = direc.normalize();
    right = target.cross(camera_up).normalize();
    camera_up = right.cross(target);

    switch (key) {

    case '1':
        target = target * cos(yaw_step) - right * sin(yaw_step);
        direc = direc.Magnitude() * target;
        break;
    case '2':
        target = target * cos(yaw_step) + right * sin(yaw_step);
        direc = direc.Magnitude() * target;
        break;
    case '3':
        target = target * cos(pitch_step) + camera_up * sin(pitch_step);
        direc = direc.Magnitude() * target;
        break;
    case '4':
        target = target * cos(pitch_step) - camera_up * sin(pitch_step);
        direc = direc.Magnitude() * target;
        break;
    case '5':
        camera_up = camera_up * cos(roll_step) + right * sin(roll_step);
        break;
    case '6':
        camera_up = camera_up * cos(roll_step) - right * sin(roll_step);
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
    glutCreateWindow("Ray Tracing");

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