/**
 * @file Particule.cpp
 * @brief Implémentations de la classe Particule et des fonctions associées.
 */

#include "Particule.hxx"
#include <iostream>
#include <cmath>

// ===============================
// Classe Particule
// ===============================

Particule::Particule(std::vector<double> pos,
                     std::vector<double> vitess,
                     std::vector<double> force,
                     double m,
                     int Id,
                     Categorie cat)
    : position(pos), vitesse(vitess), force(force), m(m), Id(Id), Cat(cat) {}

int Particule::getDim() const { return position.size(); }

double Particule::getPosition(int i) const { return position[i]; }
double Particule::getVitesse(int i) const  { return vitesse[i]; }
double Particule::getForce(int i) const    { return force[i]; }

void Particule::setPosition(int i, double val) { position[i] = val; }
void Particule::setVitesse(int i, double val)  { vitesse[i] = val; }
void Particule::setForce(int i, double val)    { force[i] = val; }

double Particule::getMas() const { return m; }
int    Particule::getId()  const { return Id; }

// La modification repose sur la 3ème loi de Newton :
// si la particule i exerce une force sur j, alors j exerce la même force sur i, mais dans le sens opposé
void Particule::Fij(Particule& p2){
    double rij = 0;
    double m2 =  p2.getMas();
    int dim = getDim();
    for(int i = 0; i < dim; i++){
        double res_int = std::pow(position[i]- p2.getPosition(i), 2);
        rij+=res_int;
    }
    if (rij < 1e-9) return;
    double dist  = std::sqrt(rij);
    double res = (m2*m)/pow(dist,3);
    for (int i = 0; i < dim; i++) {
        double f_comp = res * (p2.getPosition(i) - position[i]);
        force[i] += f_comp;
        p2.setForce(i, p2.getForce(i) - f_comp);
    }
}


double dist(const Particule& p1, const Particule& p2) {
    double sum = 0.0;
    for (int i = 0; i < p1.getDim(); ++i) {
        double d = p2.getPosition(i) - p1.getPosition(i);
        sum += d * d;
    }
    return std::sqrt(sum);
}

void CalculeF(std::vector<Particule>& particules) {
    // Remise à zéro des forces
    for (auto& p : particules) {
        for (int k = 0; k < p.getDim(); ++k) {
            p.setForce(k, 0.0);
        }
    }

    for (size_t i = 0; i < particules.size(); ++i) {
        for (size_t j = i + 1; j < particules.size(); ++j) {
            particules[i].Fij(particules[j]);
        }
    }
}

void Stormer_Verlet(std::vector<Particule>& particules, std::vector<std::vector<double>>& Fold) {
    double tend = 468.5;
    double dt   = 0.015;
    double t    = 0.0;
    int    dim  = particules[0].getDim();

    CalculeF(particules);

    while (t < tend) {
        t += dt;

        for (size_t i = 0; i < particules.size(); ++i) {
            Particule& p = particules[i];
            double mi = p.getMas();
            for (int k = 0; k < dim; ++k) {
                double xk  = p.getPosition(k);
                double vk  = p.getVitesse(k);
                double Fik = p.getForce(k);
                p.setPosition(k, xk + dt * (vk + 0.5 / mi * Fik * dt));
                Fold[i][k] = Fik;
            }
        }

        CalculeF(particules);

        for (size_t i = 0; i < particules.size(); ++i) {
            Particule& p = particules[i];
            double mi = p.getMas();
            for (int k = 0; k < dim; ++k) {
                double vk  = p.getVitesse(k);
                double Fik = p.getForce(k);
                p.setVitesse(k, vk + dt * 0.5 / mi * (Fik + Fold[i][k]));
            }
        }
    }

    for (size_t i = 0; i < particules.size(); ++i) {
        Particule& p = particules[i];
        std::cout << t;
        for (int k = 0; k < dim; ++k)
            std::cout << " " << p.getPosition(k);
        for (int k = 0; k < dim; ++k)
            std::cout << " " << p.getVitesse(k);
        std::cout << "\n";
    }
}