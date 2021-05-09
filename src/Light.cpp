//
// Created by User on 7/05/2021.
//

#include "Light.h"


Light::Light(const Color &ambientLight, const Color &diffuseLight, const Color &specularLight) : ambientLight(
        ambientLight), diffuseLight(diffuseLight), specularLight(specularLight) {}

InfLight::InfLight(const Color &ambientLight, const Color &diffuseLight, const Color &specularLight,
                   Vector3D ldVector) : Light(ambientLight, diffuseLight, specularLight), ldVector(ldVector) {
}

const Vector3D &InfLight::getLdVector() const {
    return ldVector;
}
