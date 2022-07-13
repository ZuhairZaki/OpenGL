
#include<iostream>
#include<fstream>
#include<cmath>
#include<stack>

#define pi (2*acos(0.0))

using namespace std;

class point
{
public:
    double x, y, z, w;

    point();
    point(double,double,double);
    point(double,double,double,double);
    void make_homogeneous();
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

class Matrix
{
public:
     int rows,cols;
     double** mat;

     Matrix(int,int);
     Matrix(const Matrix&);
    ~ Matrix();
    void transpose();
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

    for(int i=0;i<rows;i++)
        free(mat[i]);
    free(mat);

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
    rotAxis.normalize();

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

Vector eye, look, up;
double fovY, aspect, zNear, zFar;

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
            modelMatrix = createIdentityMatrx(4);
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
    l.normalize();
    Vector r = l.cross(up);
    r.normalize();
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
}

int main()
{   
    return 0;
}




