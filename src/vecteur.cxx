/**
 * @file vecteur.cxx
 * @brief Implémentations de la classe vecteur et de l'opérateur de flux associé.
 */

#include "vecteur.hxx"
#include <cmath>
#include <iostream>

// Classe vecteur

vecteur::vecteur(double x, double y, double z) : x(x), y(y), z(z) {}

double vecteur::getX() const { return x; }
double vecteur::getY() const { return y; }
double vecteur::getZ() const { return z; }

void vecteur::setX(double dx) { x = dx; }
void vecteur::setY(double dy) { y = dy; }
void vecteur::setZ(double dz) { z = dz; }

vecteur vecteur::sum_Vectors(vecteur v) {
    return vecteur(x + v.x, y + v.y, z + v.z);
}

vecteur vecteur::sub_Vectors(vecteur v) {
    return vecteur(x - v.x, y - v.y, z - v.z);
}

vecteur vecteur::multV_par_lambda(double lambda) {
    return vecteur(x * lambda, y * lambda, z * lambda);
}

double vecteur::dot(vecteur v) {
    return x*v.x + y*v.y + z*v.z;
}

double vecteur::norme() {
    return std::sqrt(x*x + y*y + z*z);
}

void vecteur::print() {
    std::cout << "(" << x << ", "<< y<< ", "<< z<<")" << std::endl;
}

std::ostream& operator<<(std::ostream& os, const vecteur& v) {
    os << "Vecteur(" << v.getX() << ", " << v.getY() << ", " << v.getZ() << ")";
    return os;
}