
#include<iostream>
#include<fstream>
#include<cmath>
#include<cstdlib>
#include<ctime>
#include<stack>
#include "bitmap_image.hpp"

#define pi (2*acos(0.0))

using namespace std;

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
    res.x = vec.x*scale;
    res.y = vec.y*scale;
    res.z = vec.z*scale;
    return res;
}

class point
{
public:
    double x, y, z, w;

    point();
    point(double,double,double);
    point(double,double,double,double);
    void make_homogeneous();
    Vector toVector();
};

point::point()
{
        x = 0;
        y = 0;
        z = 0;
        w = 1;
}

point::point(double px, double py, double pz)
{
        x = px;
        y = py;
        z = pz;
        w = 1;
}

point::point(double px, double py, double pz, double pw)
{   
    x = px;
    y = py;
    z = pz;
    w = pw;
}

void point::make_homogeneous()
{
    x /= w;
    y /= w;
    z /= w;
    w = 1;
}

Vector point::toVector()
{
    return Vector(x,y,z);
}

class Matrix
{
public:
     int rows,cols;
     double** mat;

     Matrix() {}
     Matrix(int,int);
     Matrix(const Matrix&);
    ~ Matrix();
    void transpose();
    Matrix& operator=(const Matrix&);
    Matrix operator + (Matrix const&);
    Matrix operator - (Matrix const&);
    Matrix operator * (Matrix const&);
    point operator * (point const&);
};

 Matrix:: Matrix(int no_rows,int no_cols)
{
    rows = no_rows;
    cols = no_cols;

    mat = (double**) malloc(rows*sizeof(double*));
    for(int i=0;i<rows;i++)
        mat[i] = (double*) malloc(cols*sizeof(double));

    for(int i=0;i<rows;i++)
        for(int j=0;j<cols;j++)
            mat[i][j] = 0;
}

Matrix::Matrix(const Matrix& m)
{   
    rows = m.rows;
    cols = m.cols;

    mat = (double**) malloc(rows*sizeof(double*));
    for(int i=0;i<rows;i++)
        mat[i] = (double*) malloc(cols*sizeof(double));

    for(int i=0;i<rows;i++)
        for(int j=0;j<cols;j++)
            mat[i][j] = m.mat[i][j];    
}

 Matrix::~ Matrix()
{
    for(int i=0;i<rows;i++)
        free(mat[i]);
    free(mat);
}

void Matrix::transpose()
{
    double** arr = mat;

    int temp = rows;
    rows = cols;
    cols = temp;

    mat = (double**) malloc(rows*sizeof(double*));
    for(int i=0;i<rows;i++)
        mat[i] = (double*) malloc(cols*sizeof(double));

    for(int i=0;i<rows;i++)
        for(int j=0;j<cols;j++)
            mat[j][i] = arr[i][j];    

    for(int i=0;i<rows;i++)
        free(arr[i]);
    free(arr);
}

Matrix& Matrix::operator=(const Matrix& m)
{   
    for(int i=0;i<rows;i++)
        free(mat[i]);
    free(mat);

    rows = m.rows;
    cols = m.cols;

    mat = (double**) malloc(rows*sizeof(double*));
    for(int i=0;i<rows;i++)
        mat[i] = (double*) malloc(cols*sizeof(double));

    for(int i=0;i<rows;i++)
        for(int j=0;j<cols;j++)
            mat[i][j] = m.mat[i][j];    

    return *this;
}

Matrix Matrix::operator+(Matrix const& m)
{
    Matrix res(rows,cols);
    for(int i=0;i<rows;i++)
        for(int j=0;j<cols;j++)
            res.mat[i][j] = mat[i][j] + m.mat[i][j];   

    return res; 
}

Matrix Matrix::operator-(Matrix const& m)
{
    Matrix res(rows,cols);
    for(int i=0;i<rows;i++)
        for(int j=0;j<cols;j++)
            res.mat[i][j] = mat[i][j] - m.mat[i][j];   

    return res; 
}

Matrix Matrix::operator*(Matrix const& m)
{
    Matrix res(rows,m.cols);
    for(int i=0;i<rows;i++)
    {
        for(int j=0;j<m.cols;j++)
        {
            double cell = 0;
            for(int k=0;k<cols;k++)
                cell += mat[i][k]*m.mat[k][j];
            res.mat[i][j] = cell;
        }
    }

    return res; 
}

point Matrix::operator*(point const& p)
{
    if(cols!=4)
        return p;
    
    point res;
    res.x = mat[0][0]*p.x + mat[0][1]*p.y + mat[0][2]*p.z + mat[0][3]*p.w;
    res.y = mat[1][0]*p.x + mat[1][1]*p.y + mat[1][2]*p.z + mat[1][3]*p.w;
    res.z = mat[2][0]*p.x + mat[2][1]*p.y + mat[2][2]*p.z + mat[2][3]*p.w;
    res.w = mat[3][0]*p.x + mat[3][1]*p.y + mat[3][2]*p.z + mat[3][3]*p.w;

    return res;
}

Matrix operator*(double scale, Matrix const& mat)
{
    Matrix res(mat.rows,mat.cols);
    for(int i=0;i<mat.rows;i++)
        for(int j=0;j<mat.cols;j++)
            res.mat[i][j] = scale*mat.mat[i][j];
    return res;
}

Matrix operator*(Matrix const& mat, double scale)
{
    Matrix res(mat.rows,mat.cols);
    for(int i=0;i<mat.rows;i++)
        for(int j=0;j<mat.cols;j++)
            res.mat[i][j] = mat.mat[i][j]*scale;
    return res;
}

Vector RodriguesFormula(Vector unitVec, Vector a, double angle)
{
    double c = cos(angle*pi/180);
    double s = sin(angle*pi/180);

    return c*unitVec + (1-c)*(a.dot(unitVec))*a + s*(a.cross(unitVec));
}

Matrix createIdentityMatrx(int n)
{
    Matrix I(n,n);
    for(int i=0;i<n;i++)
        I.mat[i][i] = 1;
    return I;
}

Matrix createTranslationMatrix(double tx, double ty, double tz)
{
    Matrix T(4,4);
    T.mat[0][0] = 1;
    T.mat[1][1] = 1;
    T.mat[2][2] = 1;
    T.mat[3][3] = 1;
    T.mat[0][3] = tx;
    T.mat[1][3] = ty;
    T.mat[2][3] = tz;
    return T;
}

Matrix createRotationMatrix(double angle, double ax, double ay, double az)
{
    Matrix R(4,4);
    R.mat[3][3] = 1;

    Vector rotAxis(ax,ay,az);
    rotAxis = rotAxis.normalize();

    Vector c1 = RodriguesFormula(Vector(1,0,0),rotAxis,angle);
    Vector c2 = RodriguesFormula(Vector(0,1,0),rotAxis,angle);
    Vector c3 = RodriguesFormula(Vector(0,0,1),rotAxis,angle);

    R.mat[0][0] = c1.x;
    R.mat[1][0] = c1.y;
    R.mat[2][0] = c1.z;

    R.mat[0][1] = c2.x;
    R.mat[1][1] = c2.y;
    R.mat[2][1] = c2.z;

    R.mat[0][2] = c3.x;
    R.mat[1][2] = c3.y;
    R.mat[2][2] = c3.z;

    return R;
}

Matrix createScalingMatrix(double sx, double sy, double sz)
{
    Matrix S(4,4);
    S.mat[0][0] = sx;
    S.mat[1][1] = sy;
    S.mat[2][2] = sz;
    S.mat[3][3] = 1;
    return S;
}

struct Color
{
    double r,g,b;
};

struct  Triangle{
    point Points[3];
    Color color;
};

Vector eye, look, up;
double fovY, aspect, zNear, zFar;

double screen_width, screen_height;
double x_lim, y_lim;
double zFront, zRear;

double top_scanline(Triangle tri, double top_y, double del_y)
{
    double max_y = tri.Points[0].y;
    for(int i=1;i<3;i++)
    {
        if(tri.Points[i].y > max_y)
            max_y = tri.Points[i].y;
    }

    for(int i=0;i<screen_height;i++)
    {
        double cell_y = top_y - i*del_y;
        if(max_y > cell_y)
            return i;
    }

    return screen_height;
}

double bottom_scanline(Triangle tri, double bottom_y, double del_y)
{
    double min_y = tri.Points[0].y;
    for(int i=1;i<3;i++)
    {
        if(tri.Points[i].y < min_y)
            min_y = tri.Points[i].y;
    }

    for(int i=0;i<screen_height;i++)
    {
        double cell_y = bottom_y + i*del_y;
        if(min_y < cell_y)
            return screen_height - i - 1;
    }

    return -1;
}

double* left_right_scanline(Triangle tri, double row, double left_x, double right_x, double del_x)
{
    Vector v1(tri.Points[1].x - tri.Points[0].x,tri.Points[1].y - tri.Points[0].y,0);
    Vector v2(tri.Points[2].x - tri.Points[1].x,tri.Points[2].y - tri.Points[1].y,0);
    Vector v3(tri.Points[0].x - tri.Points[2].x,tri.Points[0].y - tri.Points[2].y,0);

    double t1 = (row - tri.Points[0].y)/v1.y;
    double t2 = (row - tri.Points[1].y)/v2.y;
    double t3 = (row - tri.Points[2].y)/v3.y;

    double max_x = left_x -1, min_x = right_x + 1;
    if(t1 > 0 && t1 < 1)
    {
        double x = tri.Points[0].x + t1*v1.x;
        if(x > max_x)
            max_x = x;
        if(x < min_x)
            min_x = x;
    }
    if(t2 > 0 && t2 < 1)
    {
        double x = tri.Points[1].x + t2*v2.x;
        if(x > max_x)
            max_x = x;
        if(x < min_x)
            min_x = x;
    }
    if(t3 > 0 && t3 < 1)
    {
        double x = tri.Points[2].x + t3*v3.x;
        if(x > max_x)
            max_x = x;
        if(x < min_x)
            min_x = x;
    }

    double* res = new double[2];
    res[0] = screen_width;
    res[1] = -1;

    for(int i=0;i<screen_width;i++)
    {
        double cell_x = left_x + i*del_x;
        if(min_x < cell_x){
            res[0] = i;
            break;
        }
    }

    for(int i=0;i<screen_width;i++)
    {
        double cell_x = right_x - i*del_x;
        if(max_x > cell_x){
            res[1] = screen_width - i - 1;
            break;
        }
    }

    return res;
}

void stage1()
{
    ifstream fin("scene.txt");
    ofstream fout("stage1.txt");

    fin >> eye.x >> eye.y >> eye.z;
    fin >> look.x >> look.y >> look.z;
    fin >> up.x >> up.y >> up.z;
    fin >> fovY >> aspect >> zNear >> zFar;

    stack<Matrix> matStack;
    Matrix modelMatrix = createIdentityMatrx(4);
    while (true)
    {
        string command;
        fin >> command;

        if(command=="triangle")
        {   
            double x,y,z;
            for(int i=0;i<3;i++){
                fin >> x >> y >> z;
                point p(x,y,z);
                p = modelMatrix*p;
                fout << p.x << " " << p.y << " " << p.z << endl;
            }
        }
        else if(command=="translate")
        {
            double x,y,z;
            fin >> x >> y >> z;
            Matrix T = createTranslationMatrix(x,y,z);
            modelMatrix = modelMatrix*T;
        }
        else if(command=="scale")
        {   
            double x,y,z;
            fin >> x >> y >> z;
            Matrix S = createScalingMatrix(x,y,z);
            modelMatrix = modelMatrix*S;
        }
        else if (command=="rotate")
        {
            double angle,ax,ay,az;
            fin >> angle >> ax >> ay >> az;
            Matrix R = createRotationMatrix(angle,ax,ay,az);
            modelMatrix = modelMatrix*R;
        }
        else if(command=="push")
        {
            matStack.push(modelMatrix);
        }
        else if(command=="pop")
        {
            modelMatrix = matStack.top();
            matStack.pop();
        }
        else if(command=="end")
        {
            break;
        }
        
    }

    fin.close();
    fout.close();
}

void stage2()
{
    ifstream fin("stage1.txt");
    ofstream fout("stage2.txt");

    Vector l = look - eye;
    l = l.normalize();
    Vector r = l.cross(up);
    r = r.normalize();
    Vector u = r.cross(l);

    Matrix T = createTranslationMatrix(-eye.x,-eye.y,-eye.z);

    Matrix R(4,4);
    R.mat[3][3] = 1;

    R.mat[0][0] = r.x;
    R.mat[0][1] = r.y;
    R.mat[0][2] = r.z;

    R.mat[1][0] = u.x;
    R.mat[1][1] = u.y;
    R.mat[1][2] = u.z;

    R.mat[2][0] = -l.x;
    R.mat[2][1] = -l.y;
    R.mat[2][2] = -l.z;

    Matrix V = R*T;

    double x,y,z;
    while (fin>>x>>y>>z)
    {
        point p(x,y,z);
        p = V*p;
        fout << p.x << " " << p.y << " " << p.z << endl;
    }

    fin.close();
    fout.close();
}

void stage3()
{
    ifstream fin("stage2.txt");
    ofstream fout("stage3.txt");

    double fovX = fovY*aspect;
    double t = zNear*tan(fovY*pi/360);
    double r = zNear*tan(fovX*pi/360);

    Matrix P(4,4);
    P.mat[0][0] = zNear/r;
    P.mat[1][1] = zNear/t;
    P.mat[2][2] = -(zFar+zNear)/(zFar-zNear);
    P.mat[2][3] = -(2*zFar*zNear)/(zFar-zNear);
    P.mat[3][2] = -1;

    double x,y,z;
    while (fin>>x>>y>>z)
    {
        point p(x,y,z);
        p = P*p;
        p.make_homogeneous();
        fout << p.x << " " << p.y << " " << p.z << endl;
    }

    fin.close();
    fout.close();
}

void stage4()
{   
    ifstream conf("config.txt");

    conf >> screen_width >> screen_height;
    conf >> x_lim;
    conf >> y_lim;
    conf >> zFront >> zRear;

    conf.close();

    ifstream fin("stage3.txt");
    ofstream zout("z_buffer.txt");
    
    double** zBuffer = (double**) malloc(screen_height*sizeof(double*));
    for(int i=0;i<screen_height;i++)
        zBuffer[i] = (double*) malloc(screen_width*sizeof(double));

    for(int i=0;i<screen_height;i++)
        for(int j=0;j<screen_width;j++)
            zBuffer[j][i] = zRear;    

    bitmap_image image(screen_width,screen_height);
    for(int i=0;i<screen_height;i++)
        for(int j=0;j<screen_width;j++)
            image.set_pixel(i,j,0,0,0);

    double dx = (-x_lim*2)/screen_width;
    double dy = (-y_lim*2)/screen_height;
    double Top_y = -y_lim - dy/2;
    double Left_x = x_lim + dx/2;

    while(!fin.eof())
    {
        Triangle t;
        for(int i=0;i<3;i++)
        {   
            double x,y,z;
            fin >> x >> y >> z;
            
            if(fin.eof())
                break;
            
            point p(x,y,z);
            t.Points[i] = p;
        }

        if(fin.eof())
            break;

        t.color.r = rand()%256;
        t.color.g = rand()%256;
        t.color.b = rand()%256;

        Vector n = (t.Points[1].toVector() - t.Points[0].toVector()).cross(t.Points[2].toVector() - t.Points[0].toVector());
        n = n.normalize();
        double d = n.dot(t.Points[0].toVector());

        double top_row = top_scanline(t,Top_y,dy);
        double bottom_row = bottom_scanline(t,-Top_y,dy);
        for(int i=top_row;i<=bottom_row;i++)
        {   
            double cell_y = Top_y - i*dy;
            double* cols = left_right_scanline(t,cell_y,Left_x,-Left_x,dx);

            for(int j=cols[0];j<=cols[1];j++)
            {   
                double cell_x = Left_x + j*dx;
                double cell_z = (d - n.x*cell_x - n.y*cell_y)/n.z;
                
                if(cell_z > zFront)
                {
                    if(cell_z < zBuffer[j][i])
                    {
                        zBuffer[j][i] = cell_z;
                        image.set_pixel(j,i,t.color.r,t.color.g,t.color.b);
                    }
                }
            }

            delete[] cols;
        }
    }

    for(int i=0;i<screen_height;i++)
    {
        for(int j=0;j<screen_width;j++)
            zout << zBuffer[i][j] << " ";
        zout<<endl;
    }

    image.save_image("output.bmp");

    for(int i=0;i<screen_height;i++)
        delete[] zBuffer[i];
    delete zBuffer;

    fin.close();
    zout.close();
}

int main()
{   
    //srand(time(NULL));

    stage1();
    stage2();
    stage3();
    stage4();
    return 0;
}




