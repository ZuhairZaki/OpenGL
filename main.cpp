#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>

#include <GL\freeglut.h>
#include "classes.hpp"
#include "bitmap_image.hpp"

#define pi (2*acos(0.0))
#define T_MAX 100000

double windowHeight, windowWidth, viewAngle;
int capture_no;

double camera_speed;
double yaw_step, pitch_step, roll_step;

Vector camera_pos, direc, camera_up;

int recursion_level, img_dim, no_objects, no_point_light, no_spot_light;

std::vector<Object*> objects;
std::vector<PointLight> pointLights;
std::vector<SpotLight> spotLights;

void initGL()
{   
    windowHeight = 100;
    windowWidth = 100;
    viewAngle = 80.0;
    capture_no = 11;

    camera_speed = 4.0;
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

void loadData()
{
    std::ifstream fin("scene.txt");

    fin >> recursion_level;
    fin >> img_dim;
    fin >> no_objects;

    for (int i = 0; i < no_objects; i++)
    {
        std::string obj;
        fin >> obj;

        if (obj == "sphere")
        {
            Vector ref_point;
            double sphere_radius;
            fin >> ref_point.x >> ref_point.y >> ref_point.z;
            fin >> sphere_radius;
            Sphere *s = new Sphere(ref_point,sphere_radius);

            double sphere_color[3], sphere_coefficients[4], sphere_shine;
            for (int i = 0; i < 3; i++)
                fin >> sphere_color[i];
            for (int i = 0; i < 4; i++)
                fin >> sphere_coefficients[i];
            fin >> sphere_shine;
            s->SetColor(sphere_color);
            s->SetCoefficients(sphere_coefficients);
            s->SetShine(sphere_shine);

            objects.push_back(s);
        }
        else if (obj == "triangle")
        {
            Vector tri_points[3];
            for(int i=0;i<3;i++)
                fin >> tri_points[i].x >> tri_points[i].y >> tri_points[i].z;
            Triangle* t = new Triangle(Vector(0, 0, 0), tri_points);

            double tri_color[3], tri_coefficients[4], tri_shine;
            for (int i = 0; i < 3; i++)
                fin >> tri_color[i];
            for (int i = 0; i < 4; i++)
                fin >> tri_coefficients[i];
            fin >> tri_shine;
            t->SetColor(tri_color);
            t->SetCoefficients(tri_coefficients);
            t->SetShine(tri_shine);

            objects.push_back(t);
        }
        else if (obj == "general")
        {
            Vector ref_point;
            double coeff[10],bound_dim[3];

            for (int i = 0; i < 10; i++)
                fin >> coeff[i];
            fin >> ref_point.x >> ref_point.y >> ref_point.z;
            for (int i = 0; i < 3; i++)
                fin >> bound_dim[i];

            double gen_color[3], gen_coefficients[4], gen_shine;
            for (int i = 0; i < 3; i++)
                fin >> gen_color[i];
            for (int i = 0; i < 4; i++)
                fin >> gen_coefficients[i];
            fin >> gen_shine;

            Quadric* q = new Quadric(ref_point, bound_dim, coeff);
            q->SetColor(gen_color);
            q->SetCoefficients(gen_coefficients);
            q->SetShine(gen_shine);

            objects.push_back(q);
        }
    }

    Floor *f = new Floor(Vector(0, 0, 0), 1000, 20);
    double floor_coefficients[4] = { 0.4,0.2,0.2,0.2 };
    f->SetCoefficients(floor_coefficients);
    f->SetShine(10);
    objects.push_back(f);

    fin >> no_point_light;
    for (int i = 0; i < no_point_light; i++)
    {
        Vector light_pos;
        double light_color[3];

        fin >> light_pos.x >> light_pos.y >> light_pos.z;
        for (int i = 0; i < 3; i++)
            fin >> light_color[i];

        PointLight pl(light_pos);
        pl.SetColor(light_color);
        pointLights.push_back(pl);
    }

    fin >> no_spot_light;
    for (int i = 0; i < no_spot_light; i++)
    {
        Vector light_pos, light_dir;
        double light_color[3],cutoff;

        fin >> light_pos.x >> light_pos.y >> light_pos.z;
        for (int i = 0; i < 3; i++)
            fin >> light_color[i];
        fin >> light_dir.x >> light_dir.y >> light_dir.z;
        fin >> cutoff;

        SpotLight sl(light_pos, light_dir, light_color, cutoff);
        spotLights.push_back(sl);
    }

    fin.close();
}

void capture()
{
    bitmap_image image(img_dim, img_dim);
    for (int i = 0; i < img_dim; i++)
        for (int j = 0; j < img_dim; j++)
            image.set_pixel(i, j, 1.0, 0, 0);

    double plane_dist = (windowHeight / 2.0) / (tan(viewAngle * pi / 360));

    Vector l = direc.normalize();
    Vector u = camera_up.normalize();
    Vector r = l.cross(u).normalize();
    Vector top_left = camera_pos + l * plane_dist - r * (windowWidth / 2) + u * (windowHeight / 2);

    double du = windowWidth / img_dim;
    double dv = windowHeight / img_dim;

    top_left = top_left + r * 0.5 * du - u * 0.5 * dv;

    std::cout << "Capture Started" << std::endl;

    double t_min;
    for (int i = 0; i < img_dim; i++)
    {
        for (int j = 0; j < img_dim; j++)
        {
            Vector cur_pixel = top_left + r * i * du - u * j * dv;

            Ray r(camera_pos, (cur_pixel - camera_pos).normalize());

            int nearest_obj = -1;
            double t = T_MAX;

            double dummycolor[3] = { 0,0,0 };
            for (int k=0;k<objects.size();k++)
            {
                double t_val = objects[k]->intersect(&r, dummycolor, 0);
                if ((t_val > 0) && (t_val < t))
                {
                    t = t_val;
                    nearest_obj = k;
                }
            }

            if (nearest_obj == -1)
                continue;

            double* pixel_color = new double[3];
            for (int  k = 0; k < 3; k++)
                pixel_color[k] = 0;

            t_min = objects[nearest_obj]->intersect(&r, pixel_color, 1);
            image.set_pixel(i,j, pixel_color[0]*255 , pixel_color[1] * 255, pixel_color[2]*255);

            //std::cout << i << " " << j << std::endl;

            delete[] pixel_color;
        }
    }

    std::cout << "Capture Finished" << std::endl;

    std::string capture_file = "output_" + std::to_string(capture_no) + ".bmp";
    capture_no++;

    image.save_image(capture_file);
}

void display(void) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    gluLookAt(camera_pos.x, camera_pos.y, camera_pos.z,
        camera_pos.x + direc.x, camera_pos.y + direc.y, camera_pos.z + direc.z,
        camera_up.x, camera_up.y, camera_up.z);

    for (auto it = objects.begin(); it != objects.end(); it++)
        (*it)->draw();

    for (auto it = pointLights.begin(); it != pointLights.end(); it++)
    {
        glPointSize(5);
        glColor3f(it->color[0], it->color[1], it->color[2]);
        glBegin(GL_POINTS); {
            glVertex3f(it->pos.x, it->pos.y, it->pos.z);
        } glEnd();
    }

    for (auto it = spotLights.begin(); it != spotLights.end(); it++)
    {
        glPointSize(5);
        glColor3f(it->src.color[0], it->src.color[1], it->src.color[2]);
        glBegin(GL_POINTS); {
            glVertex3f(it->src.pos.x, it->src.pos.y, it->src.pos.z);
        } glEnd();
    }

    glutSwapBuffers();
}

void keyboardListener(unsigned char key, int x, int y) {
    Vector target, right;
    target = direc.normalize();
    right = target.cross(camera_up).normalize();
    camera_up = right.cross(target);

    switch (key) {
    
    case '0':
        capture();
        break;
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
    glutInitWindowSize(640,640);
    glutInitWindowPosition(200, 200);
    glutCreateWindow("Ray Tracing");

    initGL();
    loadData();

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutIdleFunc(animate);

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutMouseFunc(mouseListener);

    glutMainLoop();
    return 0;
}