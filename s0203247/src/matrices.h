//
// Created by User on 27/03/2021.
//

#ifndef ENGINE_MATRICES_H
#define ENGINE_MATRICES_H
Matrix eyePointTransformationMatrix(const double alpha, const double beta, const double r){
    Matrix m = Matrix();
    m(1,1) = -sin(alpha);
    m(1,2) = -cos(alpha)*cos(beta);
    m(1,3) = cos(alpha)*sin(beta);
    m(2,1) = cos(alpha);
    m(2,2) = -sin(alpha)*cos(beta);
    m(2,3) = sin(alpha)*sin(beta);
    m(3,2) = sin(beta);
    m(3,3) = cos(beta);
    m(4,3) = -r;
    return m;
}

Matrix scalingMatrix(const double scale){
    Matrix S = Matrix();
    for (int i = 1; i <= 3; i++)
        S(i,i) = scale;
    return S;
}
Matrix rotateX(const double alpha){
    Matrix m_x = Matrix();
    m_x(2,2) = cos(alpha);
    m_x(3,3) = cos(alpha);
    m_x(3,2) = -sin(alpha);
    m_x(2,3) = sin(alpha);
    return m_x;
}
Matrix rotateY(const double alpha){
    Matrix m_y = Matrix();
    m_y(1,1) = cos(alpha);
    m_y(1,3) = -sin(alpha);
    m_y(3,1) = sin(alpha);
    m_y(3,3) = cos(alpha);
    return m_y;
}

Matrix rotateZ(const double alpha){
    Matrix m_z = Matrix();
    m_z(1,1) = cos(alpha);
    m_z(2,2) = cos(alpha);
    m_z(2,1) = -sin(alpha);
    m_z(1,2) = sin(alpha);
    return m_z;
}

Matrix translate(const Vector3D &vector){
    Matrix m_t = Matrix();
    m_t(4,1)= vector.x;
    m_t(4,2)= vector.y;
    m_t(4,3)= vector.z;
    return m_t;
}

Matrix transformationMatrix(const double scale, const double alpha_x, const double alpha_y,const double alpha_z, const Vector3D &vector){
    return scalingMatrix(scale)*rotateX(alpha_x)*rotateY(alpha_y)*rotateZ(alpha_z)*translate(vector);
}
#endif //ENGINE_MATRICES_H
