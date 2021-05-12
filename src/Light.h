//
// Created by User on 7/05/2021.
//

#ifndef ENGINE_LIGHT_H
#define ENGINE_LIGHT_H
#include "vector3d.h"
#include "Color.h"
#include <list>
#include <cmath>

class Light
{
public:
    Light(const Color &ambientLight, const Color &diffuseLight, const Color &specularLight);

//de ambiente licht component
    Color ambientLight;
    //de diffuse licht component
    Color diffuseLight;
    //de diffuse licht component
    Color specularLight;
    virtual double getSpotAngle() const {};

    virtual bool isInfLight() const {
        return false;
    }
    virtual bool isPtLight() const {
        return false;
    }

    bool hasAmb(){
        return (!ambientLight.isZero());
    }
    bool hasDiff(){
        return (!diffuseLight.isZero());
    }
    bool hasSpec(){
        return (!specularLight.isZero());
    }
    virtual const Vector3D &getLdVector()const{};

    virtual const Vector3D &getLocation()const;
    virtual void applyTransformation(const Matrix & m){}
};

class InfLight: public Light
{
public:
    InfLight(const Color &ambientLight, const Color &diffuseLight, const Color &specularLight,
             const Vector3D& ldVector);

    //de richting waarin het
    //licht schijnt
    Vector3D ldVector;
    bool isInfLight() const override {
        return true;
    }
    void applyTransformation(const Matrix & m) override{
        ldVector*=m;
    }

    const Vector3D &getLdVector() const override;
};

class PointLight: public Light
{
public:
    PointLight(const Color &ambientLight, const Color &diffuseLight, const Color &specularLight,
               const Vector3D &location, const double &spotAngle);

    //de locatie van de puntbron
    Vector3D location;
    //de hoek van een spotlicht
    double spotAngle;
    ZBuffer shadow_mask;
    Matrix eye;
    double d,dx,dy;

    const Vector3D &getLocation() const override;
    void applyTransformation(const Matrix & m) override{
        location*=m;
    }
    bool isPtLight() const override{
        return true;
    }

    double getSpotAngle() const override;
};


typedef std::list<Light*> Lights3D;





#endif //ENGINE_LIGHT_H
