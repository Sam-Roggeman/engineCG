//
// Created by User on 7/05/2021.
//

#ifndef ENGINE_LIGHT_H
#define ENGINE_LIGHT_H
#include "vector3d.h"
#include "Color.h"
#include <list>

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

    virtual bool isInfLight() const {
        return false;
    }

    const virtual Vector3D &getLdVector()const{};
    virtual void applyTransformation(const Matrix & m){}
};

class InfLight: public Light
{
public:
    InfLight(const Color &ambientLight, const Color &diffuseLight, const Color &specularLight,
             Vector3D ldVector);

    //de richting waarin het
    //licht schijnt
    Vector3D ldVector;
    bool isInfLight() const override {
        return true;
    }
    virtual void applyTransformation(const Matrix & m){
        ldVector*=m;
    }

    const Vector3D &getLdVector() const override;
};

class PointLight: public Light
{
public:
    //de locatie van de puntbron
    Vector3D location;
    //de hoek van een spotlicht
    double spotAngle;
};
typedef std::list<Light*> Lights3D;





#endif //ENGINE_LIGHT_H
