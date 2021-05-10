//
// Created by User on 7/05/2021.
//

#include "Light.h"


Light::Light(const Color &ambientLight, const Color &diffuseLight, const Color &specularLight) : ambientLight(
        ambientLight), diffuseLight(diffuseLight), specularLight(specularLight) {}

const Vector3D &Light::getLocation() const {  }

InfLight::InfLight(const Color &ambientLight, const Color &diffuseLight, const Color &specularLight,
                   const Vector3D& ldVector) : Light(ambientLight, diffuseLight, specularLight), ldVector(ldVector){
}

const Vector3D &InfLight::getLdVector() const {
    return ldVector;
}

const Vector3D &PointLight::getLocation() const {
    return location;
}

PointLight::PointLight(const Color &ambientLight, const Color &diffuseLight, const Color &specularLight,
                       const Vector3D &location,const double &spotAngle):Light(ambientLight, diffuseLight, specularLight),location(location), spotAngle(M_PI *spotAngle/180){}

double PointLight::getSpotAngle() const {
    return spotAngle;
}
