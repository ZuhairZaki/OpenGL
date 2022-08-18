#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<vector>

#include <GL\freeglut.h>

#define pi (2*acos(0.0))
#define TOLERANCE 0.000001

struct Point {
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
    double dot(Vector const&);
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

double Vector::dot(Vector const& vec)
{
    return x * vec.x + y * vec.y + z * vec.z;
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

class Ray
{
public:
    Vector start;
    Vector dir;

    Ray() {}
    Ray(Vector, Vector);
};

Ray::Ray(Vector ray_start, Vector ray_dir)
{
    this->start = ray_start;
    this->dir = ray_dir;
}

class PointLight
{
public:
    Vector pos;
    double color[3];
    PointLight() {}
    PointLight(Vector);
    void SetLightPos(Vector);
    void SetColor(double*);
};

PointLight::PointLight(Vector light_pos)
{
    this->pos = light_pos;
}

void PointLight::SetLightPos(Vector light_pos)
{
    this->pos = light_pos;
}

void PointLight::SetColor(double* col)
{
    this->color[0] = col[0];
    this->color[1] = col[1];
    this->color[2] = col[2];
}

class SpotLight
{
public:
    PointLight src;
    Vector dir;
    double cutoff_angle;
    SpotLight() {}
    SpotLight(Vector, Vector, double*, double);
    void SetLightPos(Vector);
    void SetLightDirection(Vector);
    void SetColor(double*);
    void SetCutoffAngle(double);
};

SpotLight::SpotLight(Vector light_pos, Vector light_dir, double* col, double cutoff)
{
    src.SetLightPos(light_pos);
    src.SetColor(col);
    this->dir = light_dir;
    this->cutoff_angle = cutoff;
}

void SpotLight::SetLightPos(Vector light_pos)
{
    src.SetLightPos(light_pos);
}

void SpotLight::SetLightDirection(Vector light_dir)
{
    this->dir = light_dir;
}

void SpotLight::SetColor(double* col)
{
    src.SetColor(col);
}

void SpotLight::SetCutoffAngle(double cutoff)
{
    this->cutoff_angle = cutoff;
}

double det(double mat[][3])
{
    double res = mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
        mat[0][1] * (mat[1][0] * mat[2][2] - mat[2][0] * mat[1][2]) +
        mat[0][2] * (mat[1][0]*mat[2][1] - mat[2][0]*mat[1][1]);
    return res;
}

class Object
{   
public:
    Vector reference_point;
    double color[3];
    double coefficients[4];
    int shine;

    Object() {}
    Object(Vector);

    virtual void draw() {}
    virtual double intersect(Ray *r,double* color,int level) { return -1.0; }

    virtual double* GetColorAt(Vector p) {
        double* c = new double[3];
        c[0] = color[0];
        c[1] = color[1];
        c[2] = color[2];
        return c;
    }

    void SetColor(double*);
    void SetShine(int);
    void SetCoefficients(double*);
private:

};

Object::Object(Vector ref_point)
{   
    this->reference_point = ref_point;
}

void Object::SetColor(double* col)
{
    this->color[0] = col[0];
    this->color[1] = col[1];
    this->color[2] = col[2];
}

void Object::SetShine(int s)
{
    this->shine = s;
}

void Object::SetCoefficients(double* coef)
{
    this->coefficients[0] = coef[0];
    this->coefficients[1] = coef[1];
    this->coefficients[2] = coef[2];
    this->coefficients[3] = coef[3];
}

class Sphere : public Object
{
    double radius;
public:
    Sphere() {}
    Sphere(Vector, double);
    void draw();
    double intersect(Ray*, double*, int);
};

Sphere::Sphere(Vector ref_point, double r):Object(ref_point)
{
    this->radius = r;
}

void Sphere::draw()
{   
    glPushMatrix();
    glTranslatef(reference_point.x, reference_point.y, reference_point.z);

    int slices = 400;
    int stacks = 100;

    std::vector<std::vector<Point>> points(stacks + 1);


    for (int i = 0; i <= stacks; i++)
    {
        double z = radius * sin(((double)i / (double)stacks) * (pi / 2));
        double xy_radius = radius * cos(((double)i / (double)stacks) * (pi / 2));
        for (int j = 0; j <= slices; j++)
        {
            Point p;
            p.x = xy_radius * cos(((double)j / (double)slices) * (2*pi));
            p.y = xy_radius * sin(((double)j / (double)slices) * (2*pi));
            p.z = z;
            points[i].push_back(p);
        }
    }

    glColor3f(color[0], color[1], color[2]);
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

        for (int j = 0; j < slices; j++)
        {
            glBegin(GL_QUADS); {
                glVertex3f(points[i + 1][j].x, points[i + 1][j].y, -points[i + 1][j].z);
                glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, -points[i + 1][j + 1].z);
                glVertex3f(points[i][j + 1].x, points[i][j + 1].y, -points[i][j + 1].z);
                glVertex3f(points[i][j].x, points[i][j].y, -points[i][j].z);
            }glEnd();
        }
    }

    glPopMatrix();
}

class Triangle :public Object
{
    Vector points[3];
public:
    Triangle() {}
    Triangle(Vector,Vector*);
    void draw();
    double intersect(Ray*, double*, int);
};

Triangle::Triangle(Vector ref_point,Vector*tri_points):Object(ref_point)
{
    this->points[0] = tri_points[0];
    this->points[1] = tri_points[1];
    this->points[2] = tri_points[2];
}

void Triangle::draw()
{
    glPushMatrix();
    glTranslatef(reference_point.x, reference_point.y, reference_point.z);

    glBegin(GL_TRIANGLES); {
        glColor3f(color[0], color[1], color[2]);
        glVertex3f(points[0].x, points[0].y, points[0].z);
        glVertex3f(points[1].x, points[1].y, points[1].z);
        glVertex3f(points[2].x, points[2].y, points[2].z);
    } glEnd();

    glPopMatrix();
}

class Quadric :public Object
{   
    double eqn_coeff[10];
    double cube_dim[3];
public:
    Quadric() {}
    Quadric(Vector, double*, double*);
    double intersect(Ray*, double*, int);
};

Quadric::Quadric(Vector ref_point, double* dim, double* coeff):Object(ref_point)
{
    for (int i = 0; i < 3; i++)
        this->cube_dim[i] = dim[i];

    for (int i = 0; i < 10; i++)
        this->eqn_coeff[i] = coeff[i];
}

class Floor:public Object
{
    double floor_width;
    double tile_width;
public:
    Floor() {}
    Floor(Vector, double, double);
    void draw();
    double* GetColorAt(Vector);
    double intersect(Ray*, double*, int);
};

Floor::Floor(Vector ref_point, double fw, double tw):Object(ref_point)
{
    this->floor_width = fw;
    this->tile_width = tw;
}

void Floor::draw()
{
    glPushMatrix();
    glTranslatef(reference_point.x, reference_point.y, reference_point.z);

    int no_tiles = floor_width / (2 * tile_width);
    int row_color = 1;
    for (int i = no_tiles; i > (-no_tiles); i--)
    {   
        int col_color = row_color;
        for (int j = no_tiles; j > (-no_tiles); j--)
        {
            glBegin(GL_QUADS); {
                glColor3f((float)col_color, (float)col_color, (float)col_color);
                glVertex3f(i * tile_width, j * tile_width, 0);
                glVertex3f((i - 1) * tile_width, j * tile_width, 0);
                glVertex3f((i - 1) * tile_width, (j - 1) * tile_width, 0);
                glVertex3f(i * tile_width, (j - 1) * tile_width, 0);
            } glEnd();
            col_color = (col_color + 1) % 2;
        }
        row_color = (row_color + 1) % 2;
    }

    glPopMatrix();
}

double* Floor::GetColorAt(Vector p)
{
    int x = 0, y = 0;

    int no_tiles = floor_width / (2 * tile_width);
    for (int i = no_tiles; i > (-no_tiles); i--)
    {
        if ((p.x <= i * tile_width) && (p.x >= (i - 1) * tile_width))
        {
            x = i;
            break;
        }
    }
    for (int i = no_tiles; i > (-no_tiles); i--)
    {
        if ((p.y <= i * tile_width) && (p.y >= (i - 1) * tile_width))
        {
            y = i;
            break;
        }
    }

    int c = (1 + no_tiles - x + no_tiles - y) % 2;
    double* color = new double[3];
    color[0] = 1.0 * c;
    color[1] = 1.0 * c;
    color[2] = 1.0 * c;

    return color;
}

extern std::vector<Object*> objects;
extern std::vector<PointLight> pointLights;
extern std::vector<SpotLight> spotLights;

double Sphere::intersect(Ray* r, double* col, int level)
{
    double t;

    Vector RayOrigin = r->start - reference_point;
    double tp = -RayOrigin.dot(r->dir);

    if (tp < 0)
        t = -1.0;
    else
    {
        double d = RayOrigin.dot(RayOrigin) - tp * tp;
        if (d > radius * radius)
            t = -1.0;
        else
        {
            double td = sqrt(radius * radius - d);
            double diff = RayOrigin.dot(RayOrigin) - radius * radius;
            if (diff > TOLERANCE)
            {
                t = tp - td;
            }
            else if (diff < TOLERANCE)
            {
                t = tp + td;
            }
            else
            {
                t = 2 * tp;
            }
        }
    }

    if (level == 0)
        return t;

    Vector intersection_point = r->start + t * r->dir;

    double* intersection_color = GetColorAt(intersection_point);
    col[0] = intersection_color[0] * coefficients[0];
    col[1] = intersection_color[1] * coefficients[0];
    col[2] = intersection_color[2] * coefficients[0];

    Vector N = (intersection_point - reference_point).normalize();

    for (auto it = pointLights.begin(); it != pointLights.end(); it++)
    {
        Vector ray_direction = (intersection_point - it->pos).normalize();
        Ray point_ray(it->pos, ray_direction);

        double t_ray_min = (intersection_point.x - it->pos.x) / ray_direction.x;
        double dummyColor[3] = { 0,0,0 }, t_ray = -1.0;
        for (auto it_obj = objects.begin(); it_obj != objects.end(); it_obj++)
        {
            t_ray = (*it_obj)->intersect(&point_ray, dummyColor, 0);
            if ((t_ray > 0) && (t_ray < t_ray_min))
                break;
        }

        if ((t_ray > 0) && (t_ray < t_ray_min))
            continue;

        Vector L = (-1) * point_ray.dir;
        Vector R = 2 * (L.dot(N)) * N - L;
        Vector V = (-1) * r->dir;

        double lambertValue = L.dot(N);
        if (lambertValue < 0)
            lambertValue = 0;

        double phongValue = R.dot(V);
        if (phongValue < 0)
            phongValue = 0;

        col[0] += it->color[0] * coefficients[1] * lambertValue;
        col[1] += it->color[1] * coefficients[1] * lambertValue;
        col[2] += it->color[2] * coefficients[1] * lambertValue;

        col[0] += it->color[0] * coefficients[2] * pow(phongValue, shine);
        col[1] += it->color[1] * coefficients[2] * pow(phongValue, shine);
        col[2] += it->color[2] * coefficients[2] * pow(phongValue, shine);
    }

    for (auto it = spotLights.begin(); it != spotLights.end(); it++)
    {
        Vector ray_direction = (intersection_point - it->src.pos).normalize();
        Ray point_ray(it->src.pos, ray_direction);

        double ray_angle = acos(it->dir.dot(ray_direction)) * 180 / pi;
        if (ray_angle > it->cutoff_angle)
            continue;

        double t_ray_min = (intersection_point.x - it->src.pos.x) / ray_direction.x;
        double dummyColor[3] = { 0,0,0 }, t_ray = -1.0;
        for (auto it_obj = objects.begin(); it_obj != objects.end(); it_obj++)
        {
            t_ray = (*it_obj)->intersect(&point_ray, dummyColor, 0);
            if ((t_ray > 0) && (t_ray < t_ray_min))
                break;
        }

        if ((t_ray > 0) && (t_ray < t_ray_min))
            continue;

        Vector L = (-1) * point_ray.dir;
        Vector R = 2 * (L.dot(N)) * N - L;
        Vector V = (-1) * r->dir;

        double lambertValue = L.dot(N);
        if (lambertValue < 0)
            lambertValue = 0;

        double phongValue = R.dot(V);
        if (phongValue < 0)
            phongValue = 0;

        col[0] += it->src.color[0] * coefficients[1] * lambertValue;
        col[1] += it->src.color[1] * coefficients[1] * lambertValue;
        col[2] += it->src.color[2] * coefficients[1] * lambertValue;

        col[0] += it->src.color[0] * coefficients[2] * pow(phongValue,shine);
        col[1] += it->src.color[1] * coefficients[2] * pow(phongValue, shine);
        col[2] += it->src.color[2] * coefficients[2] * pow(phongValue, shine);
    }

    delete[] intersection_color;

    return t;
}

double Triangle::intersect(Ray* r, double* col, int level)
{
    double matA[3][3], matBeta[3][3], matGamma[3][3], matT[3][3];

    matA[0][0] = points[0].x - points[1].x;
    matA[0][1] = points[0].x - points[2].x;
    matA[0][2] = r->dir.x;
    matA[1][0] = points[0].y - points[1].y;
    matA[1][1] = points[0].y - points[2].y;
    matA[1][2] = r->dir.y;
    matA[2][0] = points[0].z - points[1].z;
    matA[2][1] = points[0].z - points[2].z;
    matA[2][2] = r->dir.z;

    matBeta[0][0] = points[0].x - r->start.x;
    matBeta[0][1] = points[0].x - points[2].x;
    matBeta[0][2] = r->dir.x;
    matBeta[1][0] = points[0].y - r->start.y;
    matBeta[1][1] = points[0].y - points[2].y;
    matBeta[1][2] = r->dir.y;
    matBeta[2][0] = points[0].z - r->start.z;
    matBeta[2][1] = points[0].z - points[2].z;
    matBeta[2][2] = r->dir.y;

    matGamma[0][0] = points[0].x - points[1].x;
    matGamma[0][1] = points[0].x - r->start.x;
    matGamma[0][2] = r->dir.x;
    matGamma[1][0] = points[0].y - points[1].y;
    matGamma[1][1] = points[0].y - r->start.y;
    matGamma[1][2] = r->dir.y;
    matGamma[2][0] = points[0].z - points[1].z;
    matGamma[2][1] = points[0].z - r->start.z;
    matGamma[2][2] = r->dir.z;

    matT[0][0] = points[0].x - points[1].x;
    matT[0][1] = points[0].x - points[2].x;
    matT[0][2] = points[0].x - r->start.x;
    matT[1][0] = points[0].y - points[1].y;
    matT[1][1] = points[0].y - points[2].y;
    matT[1][2] = points[0].y - r->start.y;
    matT[2][0] = points[0].z - points[1].z;
    matT[2][1] = points[0].z - points[2].z;
    matT[2][2] = points[0].z - r->start.z;

    double detA = det(matA);

    double beta = det(matBeta) / detA;
    double gamma = det(matGamma) / detA;
    double t = det(matT) / detA;

    if (!((beta > 0) && (gamma > 0) && (beta + gamma < 1) && (t > 0)))
        t = -1.0;

    if (level == 0)
        return t;

    Vector intersection_point = r->start + t * r->dir;

    double* intersection_color = GetColorAt(intersection_point);
    col[0] = intersection_color[0] * coefficients[0];
    col[1] = intersection_color[1] * coefficients[0];
    col[2] = intersection_color[2] * coefficients[0];

    Vector N = (points[1] - points[0]).cross(points[2] - points[0]);

    for (auto it = pointLights.begin(); it != pointLights.end(); it++)
    {
        Vector ray_direction = (intersection_point - it->pos).normalize();
        Ray point_ray(it->pos, ray_direction);

        double t_ray_min = (intersection_point.x - it->pos.x) / ray_direction.x;
        double dummyColor[3] = { 0,0,0 }, t_ray = -1.0;
        for (auto it_obj = objects.begin(); it_obj != objects.end(); it_obj++)
        {
            t_ray = (*it_obj)->intersect(&point_ray, dummyColor, 0);
            if ((t_ray > 0) && (t_ray < t_ray_min))
                break;
        }

        if ((t_ray > 0) && (t_ray < t_ray_min))
            continue;

        Vector L = (-1) * point_ray.dir;
        Vector R = 2 * (L.dot(N)) * N - L;
        Vector V = (-1) * r->dir;

        double lambertValue = L.dot(N);
        if (lambertValue < 0)
            lambertValue = 0;

        double phongValue = R.dot(V);
        if (phongValue < 0)
            phongValue = 0;

        col[0] += it->color[0] * coefficients[1] * lambertValue;
        col[1] += it->color[1] * coefficients[1] * lambertValue;
        col[2] += it->color[2] * coefficients[1] * lambertValue;

        col[0] += it->color[0] * coefficients[2] * pow(phongValue, shine);
        col[1] += it->color[1] * coefficients[2] * pow(phongValue, shine);
        col[2] += it->color[2] * coefficients[2] * pow(phongValue, shine);
    }

    for (auto it = spotLights.begin(); it != spotLights.end(); it++)
    {
        Vector ray_direction = (intersection_point - it->src.pos).normalize();
        Ray point_ray(it->src.pos, ray_direction);

        double ray_angle = acos(it->dir.dot(ray_direction)) * 180 / pi;
        if (ray_angle > it->cutoff_angle)
            continue;

        double t_ray_min = (intersection_point.x - it->src.pos.x) / ray_direction.x;
        double dummyColor[3] = { 0,0,0 }, t_ray = -1.0;
        for (auto it_obj = objects.begin(); it_obj != objects.end(); it_obj++)
        {
            t_ray = (*it_obj)->intersect(&point_ray, dummyColor, 0);
            if ((t_ray > 0) && (t_ray < t_ray_min))
                break;
        }

        if ((t_ray > 0) && (t_ray < t_ray_min))
            continue;

        Vector L = (-1) * point_ray.dir;
        Vector R = 2 * (L.dot(N)) * N - L;
        Vector V = (-1) * r->dir;

        double lambertValue = L.dot(N);
        if (lambertValue < 0)
            lambertValue = 0;

        double phongValue = R.dot(V);
        if (phongValue < 0)
            phongValue = 0;

        col[0] += it->src.color[0] * coefficients[1] * lambertValue;
        col[1] += it->src.color[1] * coefficients[1] * lambertValue;
        col[2] += it->src.color[2] * coefficients[1] * lambertValue;

        col[0] += it->src.color[0] * coefficients[2] * pow(phongValue, shine);
        col[1] += it->src.color[1] * coefficients[2] * pow(phongValue, shine);
        col[2] += it->src.color[2] * coefficients[2] * pow(phongValue, shine);
    }

    delete[] intersection_color;

    return t;
}

double Quadric::intersect(Ray* r, double* col, int level)
{
    double t;

    double coeff_a = eqn_coeff[0] * r->dir.x * r->dir.x + eqn_coeff[1] * r->dir.y * r->dir.y + eqn_coeff[2] * r->dir.z * r->dir.z +
        eqn_coeff[3] * r->dir.x * r->dir.y + eqn_coeff[4] * r->dir.x * r->dir.z + eqn_coeff[5] * r->dir.y * r->dir.z;
    double coeff_b = 2 * (eqn_coeff[0] * r->start.x * r->dir.x + eqn_coeff[1] * r->start.y * r->dir.y + eqn_coeff[2] * r->start.z * r->dir.z) +
        eqn_coeff[3] * (r->start.x * r->dir.y + r->start.y * r->dir.x) + eqn_coeff[4] * (r->start.x * r->dir.z + r->start.z * r->dir.x) + eqn_coeff[5] * (r->start.z * r->dir.y + r->start.y * r->dir.z) +
        eqn_coeff[6] * r->dir.x + eqn_coeff[7] * r->dir.y + eqn_coeff[8] * r->dir.z;
    double coeff_c = eqn_coeff[0] * r->start.x * r->start.x + eqn_coeff[1] * r->start.y * r->start.y + eqn_coeff[2] * r->start.z * r->start.z +
        eqn_coeff[3] * r->start.x * r->start.y + eqn_coeff[4] * r->start.x * r->start.z + eqn_coeff[5] * r->start.z * r->start.y +
        eqn_coeff[6] * r->start.x + eqn_coeff[7] * r->start.y + eqn_coeff[8] * r->start.z + eqn_coeff[9];

    double determinant = coeff_b * coeff_b - 4 * coeff_a * coeff_c;
    if (determinant < 0)
        t = -1.0;
    else {
        double t1 = ( -coeff_b + sqrt(determinant)) / 2 * coeff_a;
        double t2 = ( -coeff_b - sqrt(determinant)) / 2 * coeff_a;

        if (t1 < 0 && t2 < 0)
            t = -1.0;
        else
        {
            if (t1 > 0)
            {
                Vector ins = r->start + t1 * r->dir;
                if (cube_dim[0] != 0)
                    if (ins.x > (cube_dim[0] / 2) || ins.x < (-cube_dim[0] / 2))
                        t1 = -1.0;
                if (cube_dim[1] != 0)
                    if (ins.y > (cube_dim[1] / 2) || ins.y < (-cube_dim[1] / 2))
                        t1 = -1.0;
                if (cube_dim[2] != 0)
                    if (ins.z > (cube_dim[2] / 2) || ins.z < (-cube_dim[2] / 2))
                        t1 = -1.0;
            }
            if (t2 > 0)
            {
                Vector ins = r->start + t2 * r->dir;
                if (cube_dim[0] != 0)
                    if (ins.x > (cube_dim[0] / 2) || ins.x < (-cube_dim[0] / 2))
                        t2 = -1.0;
                if (cube_dim[1] != 0)
                    if (ins.y > (cube_dim[1] / 2) || ins.y < (-cube_dim[1] / 2))
                        t2 = -1.0;
                if (cube_dim[2] != 0)
                    if (ins.z > (cube_dim[2] / 2) || ins.z < (-cube_dim[2] / 2))
                        t2 = -1.0;
            }

            if (t1 < t2)
                t = t1;
            else
                t = t2;
        }
    }

    Vector intersection_point = r->start + t * r->dir;

    double* intersection_color = GetColorAt(intersection_point);
    col[0] = intersection_color[0] * coefficients[0];
    col[1] = intersection_color[1] * coefficients[0];
    col[2] = intersection_color[2] * coefficients[0];

    double Nx = 2 * eqn_coeff[0] * intersection_point.x + eqn_coeff[3] * intersection_point.y + eqn_coeff[4] * intersection_point.z + eqn_coeff[6];
    double Ny = 2 * eqn_coeff[1] * intersection_point.y + eqn_coeff[3] * intersection_point.x + eqn_coeff[5] * intersection_point.z + eqn_coeff[7];
    double Nz = 2 * eqn_coeff[2] * intersection_point.z + eqn_coeff[4] * intersection_point.x + eqn_coeff[5] * intersection_point.y + eqn_coeff[8];
    Vector N(Nx, Ny, Nz);

    for (auto it = pointLights.begin(); it != pointLights.end(); it++)
    {
        Vector ray_direction = (intersection_point - it->pos).normalize();
        Ray point_ray(it->pos, ray_direction);

        double t_ray_min = (intersection_point.x - it->pos.x) / ray_direction.x;
        double dummyColor[3] = { 0,0,0 }, t_ray = -1.0;
        for (auto it_obj = objects.begin(); it_obj != objects.end(); it_obj++)
        {
            t_ray = (*it_obj)->intersect(&point_ray, dummyColor, 0);
            if ((t_ray > 0) && (t_ray < t_ray_min))
                break;
        }

        if ((t_ray > 0) && (t_ray < t_ray_min))
            continue;

        Vector L = (-1) * point_ray.dir;
        Vector R = 2 * (L.dot(N)) * N - L;
        Vector V = (-1) * r->dir;

        double lambertValue = L.dot(N);
        if (lambertValue < 0)
            lambertValue = 0;

        double phongValue = R.dot(V);
        if (phongValue < 0)
            phongValue = 0;

        col[0] += it->color[0] * coefficients[1] * lambertValue;
        col[1] += it->color[1] * coefficients[1] * lambertValue;
        col[2] += it->color[2] * coefficients[1] * lambertValue;

        col[0] += it->color[0] * coefficients[2] * pow(phongValue, shine);
        col[1] += it->color[1] * coefficients[2] * pow(phongValue, shine);
        col[2] += it->color[2] * coefficients[2] * pow(phongValue, shine);
    }

    for (auto it = spotLights.begin(); it != spotLights.end(); it++)
    {
        Vector ray_direction = (intersection_point - it->src.pos).normalize();
        Ray point_ray(it->src.pos, ray_direction);

        double ray_angle = acos(it->dir.dot(ray_direction)) * 180 / pi;
        if (ray_angle > it->cutoff_angle)
            continue;

        double t_ray_min = (intersection_point.x - it->src.pos.x) / ray_direction.x;
        double dummyColor[3] = { 0,0,0 }, t_ray = -1.0;
        for (auto it_obj = objects.begin(); it_obj != objects.end(); it_obj++)
        {
            t_ray = (*it_obj)->intersect(&point_ray, dummyColor, 0);
            if ((t_ray > 0) && (t_ray < t_ray_min))
                break;
        }

        if ((t_ray > 0) && (t_ray < t_ray_min))
            continue;

        Vector L = (-1) * point_ray.dir;
        Vector R = 2 * (L.dot(N)) * N - L;
        Vector V = (-1) * r->dir;

        double lambertValue = L.dot(N);
        if (lambertValue < 0)
            lambertValue = 0;

        double phongValue = R.dot(V);
        if (phongValue < 0)
            phongValue = 0;

        col[0] += it->src.color[0] * coefficients[1] * lambertValue;
        col[1] += it->src.color[1] * coefficients[1] * lambertValue;
        col[2] += it->src.color[2] * coefficients[1] * lambertValue;

        col[0] += it->src.color[0] * coefficients[2] * pow(phongValue, shine);
        col[1] += it->src.color[1] * coefficients[2] * pow(phongValue, shine);
        col[2] += it->src.color[2] * coefficients[2] * pow(phongValue, shine);
    }

    delete[] intersection_color;

    return t;
}


double Floor::intersect(Ray* r, double* col, int level)
{
    double t = -r->start.z / r->dir.z;

    double x = r->start.x + t * r->dir.x;
    if (x > (floor_width / 2) || x < (-floor_width / 2))
        t = -1.0;

    double y = r->start.y + t * r->dir.y;
    if (y > (floor_width / 2) || y < (-floor_width / 2))
        t = -1.0;

    if (level == 0)
        return t;

    Vector intersection_point = r->start + t * r->dir;
    //std::cout <<  intersection_point.x << " " << intersection_point.y << " " << intersection_point.z << std::endl;

    double* intersection_color = GetColorAt(intersection_point);
    col[0] = intersection_color[0];
    col[1] = intersection_color[1];
    col[2] = intersection_color[2];

    Vector N(0,0,1);

    for (auto it = pointLights.begin(); it != pointLights.end(); it++)
    {
        Vector ray_direction = (intersection_point - it->pos).normalize();
        Ray point_ray(it->pos, ray_direction);

        double t_ray_min = (intersection_point.x - it->pos.x) / ray_direction.x;
        double dummyColor[3] = { 0,0,0 }, t_ray = -1.0;
        for (auto it_obj = objects.begin(); it_obj != objects.end(); it_obj++)
        {
            t_ray = (*it_obj)->intersect(&point_ray, dummyColor, 0);
            if ((t_ray > 0) && (t_ray < t_ray_min))
                break;
        }

        if ((t_ray > 0) && (t_ray < t_ray_min))
            continue;

        Vector L = (-1) * point_ray.dir;
        Vector R = 2 * (L.dot(N)) * N - L;
        Vector V = (-1) * r->dir;

        double lambertValue = L.dot(N);
        if (lambertValue < 0)
            lambertValue = 0;

        double phongValue = R.dot(V);
        if (phongValue < 0)
            phongValue = 0;

        col[0] += it->color[0] * coefficients[1] * lambertValue;
        col[1] += it->color[1] * coefficients[1] * lambertValue;
        col[2] += it->color[2] * coefficients[1] * lambertValue;

        col[0] += it->color[0] * coefficients[2] * pow(phongValue, shine);
        col[1] += it->color[1] * coefficients[2] * pow(phongValue, shine);
        col[2] += it->color[2] * coefficients[2] * pow(phongValue, shine);
    }

    for (auto it = spotLights.begin(); it != spotLights.end(); it++)
    {
        Vector ray_direction = (intersection_point - it->src.pos).normalize();
        Ray point_ray(it->src.pos, ray_direction);

        double ray_angle = acos(it->dir.dot(ray_direction)) * 180 / pi;
        if (ray_angle > it->cutoff_angle)
            continue;

        double t_ray_min = (intersection_point.x - it->src.pos.x) / ray_direction.x;
        double dummyColor[3] = { 0,0,0 }, t_ray = -1.0;
        for (auto it_obj = objects.begin(); it_obj != objects.end(); it_obj++)
        {
            t_ray = (*it_obj)->intersect(&point_ray, dummyColor, 0);
            if ((t_ray > 0) && (t_ray < t_ray_min))
                break;
        }

        if ((t_ray > 0) && (t_ray < t_ray_min))
            continue;

        Vector L = (-1) * point_ray.dir;
        Vector R = 2 * (L.dot(N)) * N - L;
        Vector V = (-1) * r->dir;

        double lambertValue = L.dot(N);
        if (lambertValue < 0)
            lambertValue = 0;

        double phongValue = R.dot(V);
        if (phongValue < 0)
            phongValue = 0;

        col[0] += it->src.color[0] * coefficients[1] * lambertValue;
        col[1] += it->src.color[1] * coefficients[1] * lambertValue;
        col[2] += it->src.color[2] * coefficients[1] * lambertValue;

        col[0] += it->src.color[0] * coefficients[2] * pow(phongValue, shine);
        col[1] += it->src.color[1] * coefficients[2] * pow(phongValue, shine);
        col[2] += it->src.color[2] * coefficients[2] * pow(phongValue, shine);
    }

    delete[] intersection_color;

    return t;
}


