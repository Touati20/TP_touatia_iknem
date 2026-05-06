/**
 * @file univers.cxx
 * @brief Implémentations des classes Cellule, Grille et Univers.
 */

#include "Univers.hxx"
#include <iostream>
#include <cmath>
#include <algorithm>

// classe cellule

Cellule::Cellule() : cx(0.0), cy(0.0) {}
Cellule::Cellule(double x, double y) : cx(x), cy(y) {}

double Cellule::getCx() const         { return cx; }
double Cellule::getCy() const         { return cy; }
int    Cellule::getNbVoisins() const  { return nb_voisins; }
int    Cellule::nbParticules() const  { return particules_internes.size(); }

void Cellule::setNbVoisin(int dim) {
    nb_voisins = static_cast<int>(std::pow(3, dim));
}

void Cellule::ajouterVoisin(int index) {
    indices_voisines.push_back(index);
}

void Cellule::ajouterParticule(const Particule& p) {
    particules_internes.push_back(p);
}

void Cellule::vider() {
    particules_internes.clear();
}

std::vector<Particule>& Cellule::getParticules() {
    return particules_internes;
}

const std::vector<Particule>& Cellule::getParticules() const {
    return particules_internes;
}

const std::vector<int>& Cellule::getIndicesVoisines() const {
    return indices_voisines;
}

Particule& Cellule::operator[](int i) {
    return particules_internes[i];
}

// classe Grille

Grille::Grille() : nx(1), ny(1), taille_cellule(3.0) {}

Grille::Grille(int x, int y, double taille)
    : nx(x), ny(y), taille_cellule(taille) {
    cellules.resize(x * y);
    for (int j = 0; j < y; ++j){
        for (int i = 0; i < x; ++i){
            cellules[i + j * x] = Cellule((i + 0.5) * taille, (j + 0.5) * taille);
        }
    }
}

int Grille::getNx() const { return nx; }
int Grille::getNy() const { return ny; }

int Grille::getIndex(int x, int y) const {
    return x + nx * y;
}

Cellule& Grille::getCellule(int x, int y) {
    return cellules[getIndex(x, y)];
}

std::vector<Cellule>& Grille::getCellules() {
    return cellules;
}

const std::vector<Cellule>& Grille::getCellules() const {
    return cellules;
}

void Grille::initialiserVoisins() {
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            int currentIdx = x + y * nx;
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dx = -1; dx <= 1; ++dx) {
                    int nx_voisin = x + dx;
                    int ny_voisin = y + dy;
                    if (nx_voisin >= 0 && nx_voisin < nx &&
                        ny_voisin >= 0 && ny_voisin < ny) {
                        int voisinIdx = nx_voisin + ny_voisin * nx;
                        cellules[currentIdx].ajouterVoisin(voisinIdx);
                    }
                }
            }
        }
    }
}

// classe Univers

Univers::Univers(int dim, int nb_part, double Lx, double Ly, double rcut)
    : dim(dim), nb_particule(0), Lx(Lx), Ly(Ly), rcut(rcut) {
    int ncd_x = std::max(1, static_cast<int>(Lx / rcut));
    int ncd_y = std::max(1, static_cast<int>(Ly / rcut));
    grille = Grille(ncd_x, ncd_y, rcut);
    grille.initialiserVoisins();
}


int Univers::getNbPart() const { return nb_particule; }

std::vector<Cellule>& Univers::getCellules() {
    return grille.getCellules();
}

const std::vector<Cellule>& Univers::getCellules() const {
    return grille.getCellules();
}

void Univers::ajouterParticule(const Particule& p) {
    double cell_dx = Lx / grille.getNx();
    double cell_dy = Ly / grille.getNy();
    int cx = std::max(0, std::min(static_cast<int>(p.getPosition(0) / cell_dx),
                                  grille.getNx() - 1));
    int cy = std::max(0, std::min(static_cast<int>(p.getPosition(1) / cell_dy),
                                  grille.getNy() - 1));
    grille.getCellule(cx, cy).ajouterParticule(p);
    ++nb_particule;
}

void Univers::maj_cellules() {
    double cell_dx = Lx / grille.getNx();
    double cell_dy = Ly / grille.getNy();
    auto& cellules = grille.getCellules();

    std::vector<Particule> toutes;
    for (Cellule& c : cellules){
        for (Particule& p : c){
            toutes.push_back(p);
        }
    }

    for (Cellule& c : cellules){
        c.vider();
    }

    for (auto& p : toutes) {
        int cx = std::max(0, std::min(static_cast<int>(p.getPosition(0) / cell_dx),
                                      grille.getNx() - 1));
        int cy = std::max(0, std::min(static_cast<int>(p.getPosition(1) / cell_dy),
                                      grille.getNy() - 1));
        grille.getCellule(cx, cy).ajouterParticule(p);
    }
}


void Univers::avancer_parts(double dt) {
    for (auto& c : grille.getCellules()){
        for (auto& p : c){
            for (int d = 0; d < dim; ++d){
                p.setPosition(d, p.getPosition(d) + p.getVitesse(d) * dt);
            }
        }
    }
}

void Univers::maj_vitesse(double dt) {
    for (Cellule& c : grille.getCellules()){
        for (auto& p : c){
            for (int d = 0; d < dim; ++d) {
                double a = p.getForce(d) / p.getMas();
                p.setVitesse(d, p.getVitesse(d) + a * dt);
            }
        }
    }
}

double Univers::P_Lennard_Jones(double rij) const {
    double sr6 = std::pow(sigma / rij, 6);
    return 4.0 * sr6 * (sr6 - 1.0);
}


double Univers::P_limites_reflex(double r) const {
    if (r < 1e-10) return 0.0;
    double temp = sigma / (2.0 * r);
    double t6   = std::pow(temp, 6);
    return -24.0 * epsilon * (1.0 / (2.0 * r)) * t6 * (1.0 - 2.0 * t6);
}


void Univers::all_forces() {
    auto& cellules = grille.getCellules();

    // Remise à zéro des forces
    for (Cellule& c : cellules){
        for (Particule& p : c){
            for (int d = 0; d < dim; ++d){
                p.setForce(d, 0.0);
            }
        }
    }
    // Lennard-Jones 
    for (Cellule& celluleA : cellules) {
        for (Particule& pI : celluleA) {
            for (int idxVoisin : celluleA.getIndicesVoisines()) {
                Cellule& cellVoisine = cellules[idxVoisin];

                double dxC = pI.getPosition(0) - cellVoisine.getCx();
                double dyC = pI.getPosition(1) - cellVoisine.getCy();
                double distCentre = std::sqrt(dxC * dxC + dyC * dyC);

                if (distCentre <= rcut) {
                    for (Particule& pJ : cellVoisine) {
                        if (&pI == &pJ) continue;

                        double dx = pI.getPosition(0) - pJ.getPosition(0);
                        double dy = pI.getPosition(1) - pJ.getPosition(1);
                        double r2 = dx * dx + dy * dy;
                        double r  = std::sqrt(r2);

                        if (r <= rcut) {
                            double s_r6 = std::pow(sigma / r, 6);
                            double intensite = (24.0 * epsilon / r2) * s_r6
                                             * (2.0 * s_r6 - 1.0);
                            if (intensite >  5000.0) intensite =  5000.0;
                            if (intensite < -5000.0) intensite = -5000.0;
                            pI.setForce(0, pI.getForce(0) + intensite * dx);
                            pI.setForce(1, pI.getForce(1) + intensite * dy);
                        }
                    }
                }
            }
        }
    }
}

void Univers::applique_Potentiel_gravitationnel() {
    for (Cellule& c : grille.getCellules()){
        for (Particule& p : c) {
            double df = p.getMas() * G;
            p.setForce(1, p.getForce(1) + df);
        }
    }
}


void Univers::appliquer_reflection() {
    for (Cellule& c : grille.getCellules()) {
        for (Particule& p : c) {
            for (int d = 0; d < dim; ++d) {
                double pos = p.getPosition(d);
                double L   = (d == 0) ? Lx : Ly;
                if (pos < 0) {
                    p.setPosition(d, -pos);
                    p.setVitesse(d, -p.getVitesse(d));
                } else if (pos >= L) {
                    p.setPosition(d, 2 * L - pos);
                    p.setVitesse(d, -p.getVitesse(d));
                }
            }
        }
    }
}


void Univers::appliquer_absorption() {
    for (Cellule& c : grille.getCellules()) {
        auto& parts = c.getParticules();
        for (auto it = parts.begin(); it != parts.end(); ) {
            bool est_absorbee = false;
            for (int d = 0; d < dim; ++d) {
                double pos = it->getPosition(d);
                double L   = (d == 0) ? Lx : Ly;
                if (pos < 0 || pos >= L) {
                    est_absorbee = true;
                    break;
                }
            }
            if (est_absorbee) {
                it = parts.erase(it);
                --nb_particule;
            } else {
                ++it;
            }
        }
    }
}

void Univers::appliquer_periodique() {
    for (Cellule& c : grille.getCellules()){
        for (Particule& p : c){
            for (int d = 0; d < dim; ++d) {
                double L   = (d == 0) ? Lx : Ly;
                double pos = p.getPosition(d);
                if (pos < 0)  p.setPosition(d, pos + L);
                if (pos >= L) p.setPosition(d, pos - L);
            }
        }
    }
}

void Univers::potentiel_reflex() {
    double rcut_bis = sigma * std::pow(2.0, 1.0 / 6.0);
    for (Cellule& c : grille.getCellules()) {
        for (Particule& p : c) {
            for (int d = 0; d < dim; ++d) {
                double pos = p.getPosition(d);
                double L   = (d == 0) ? Lx : Ly;

                double r = pos;
                if (r > 0 && r < rcut_bis) {
                    double df = P_limites_reflex(r);
                    p.setForce(d, p.getForce(d) - df);
                }
                r = L - pos;
                if (r > 0 && r < rcut_bis) {
                    double df = P_limites_reflex(r);
                    p.setForce(d, p.getForce(d) + df);
                }
            }
        }
    }
}


void Univers::univers_state(int step) const {
    if (step >= 0) std::cout << "--- Step " << step << " ---\n";
    int i = 1;
    for (const auto& c : grille.getCellules()) {
        for (const auto& p : c) {
            std::cout << "P" << i << " pos=("
                      << p.getPosition(0) << ", " << p.getPosition(1) << ")"
                      << " vit=("
                      << p.getVitesse(0) << ", " << p.getVitesse(1) << ")\n";
            ++i;
        }
    }
    if (i == 1) std::cout << "  (aucune particule)\n";
}

double Univers::calculer_Ec() const {
    double Ec = 0.0;
    for (const auto& c : grille.getCellules()){
        for (const auto& p : c){
            for (int d = 0; d < dim; ++d){
                Ec += 0.5 * p.getMas() * p.getVitesse(d) * p.getVitesse(d);
            }
        }
    }
    return Ec;
}

void Univers::rescaleV(double Ec_cible) {
    double beta = std::sqrt(Ec_cible / calculer_Ec());
    for (Cellule& c : grille.getCellules()){
        for (Particule& p : c){
            for (int d = 0; d < dim; ++d){
                p.setVitesse(d, beta * p.getVitesse(d));
            }
        }
    }
}