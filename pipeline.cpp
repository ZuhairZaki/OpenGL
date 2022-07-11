#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stack>

#define pi (2*acos(0.0))

class point
{
public:
    double x, y, z, w;

    point();
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

point::point(double px, double py, double pz, double pw)
{   
    x = px;
    y = py;
    z = pz;

    if(pw)
        w = pw;
    else
        w = 1;
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
}

Matrix operator*(Matrix const& mat, double scale)
{
    Matrix res(mat.rows,mat.cols);
    for(int i=0;i<mat.rows;i++)
        for(int j=0;j<mat.cols;j++)
            res.mat[i][j] = mat.mat[i][j]*scale;
}

int main()
{
    return 0;
}




