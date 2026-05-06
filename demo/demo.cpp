/**
 * @file demo.cpp
 * @brief Simulation d'impact :
 *        un disque de N1 = 395 particules tombe sur une membrane
 *        de N2 = 17227 particules.
 *
 * Paramètres :
 *   L1 = 250, L2 = 180
 *   ε = 1, σ = 1, m = 1, v = (0, 10)
 *   N1 = 395, N2 = 17227
 *   rcut = 2.5σ, δt = 0.0005
 *   G = -12, Ec_D = 0.005 * (N1 + N2)
 */

#include "Particule.hxx"
#include "Univers.hxx"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <chrono>
#include <stdexcept>
#include <filesystem>

std::vector<Particule*> collecterToutesLesParticules(Univers& univers) {
    std::vector<Particule*> liste_pointeurs;
    for (Cellule& cell : univers.getCellules())
        for (Particule& p : cell)
            liste_pointeurs.push_back(&p);
    return liste_pointeurs;
}

void sauvegarderVTK(const std::string& nom, Univers& univ) {
    std::ofstream out(nom);
    if (!out.is_open())
        throw std::runtime_error("Erreur : impossible d'ouvrir " + nom);

    std::vector<Particule*> ptrs = collecterToutesLesParticules(univ);
    int nbPart = ptrs.size();

    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
        << "  <UnstructuredGrid>\n"
        << "    <Piece NumberOfPoints=\"" << nbPart << "\" NumberOfCells=\"0\">\n";

    out << "      <Points>\n"
        << "        <DataArray Name=\"Points\" type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (Particule* p : ptrs)
        out << p->getPosition(0) << " " << p->getPosition(1) << " 0 ";
    out << "\n        </DataArray>\n      </Points>\n";

    out << "      <PointData Vectors=\"Velocity\">\n"
        << "        <DataArray Name=\"Velocity\" type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (Particule* p : ptrs)
        out << p->getVitesse(0) << " " << p->getVitesse(1) << " 0 ";
    out << "\n        </DataArray>\n";

    out << "        <DataArray Name=\"Masse\" type=\"Float32\" format=\"ascii\">\n";
    for (Particule* p : ptrs)
        out << p->getMas() << " ";
    out << "\n        </DataArray>\n      </PointData>\n";

    out << "      <Cells>\n"
        << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\"></DataArray>\n"
        << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\"></DataArray>\n"
        << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\"></DataArray>\n"
        << "      </Cells>\n";

    out << "    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n";
}

void Stormer_Verlet(Univers& univers, int nb_part, int dim, double tend, double dt, double Ec_cible) {

    std::vector<double> Fold(nb_part * dim, 0.0);

    univers.all_forces();
    univers.applique_Potentiel_gravitationnel();
    univers.appliquer_reflection();

    double t = 0.0;
    int    step = 0;

    sauvegarderVTK("frames/frame_0.vtu", univers);

    while (t < tend) {
        t += dt;
        ++step;

        for (Cellule& cell : univers.getCellules()) {
            for (Particule& p : cell) {
                int    i  = p.getId();
                double mi = p.getMas();
                for (int k = 0; k < dim; ++k) {
                    double xi = p.getPosition(k);
                    double vi = p.getVitesse(k);
                    double Fi = p.getForce(k);
                    p.setPosition(k, xi + dt * (vi + 0.5 / mi * Fi * dt));
                    Fold[i * dim + k] = Fi;
                }
            }
        }

        univers.maj_cellules();
        univers.all_forces();
        univers.applique_Potentiel_gravitationnel();
        univers.appliquer_reflection();

        for (Cellule& cell : univers.getCellules()) {
            for (Particule& p : cell) {
                int    i  = p.getId();
                double mi = p.getMas();
                for (int k = 0; k < dim; ++k) {
                    double vi     = p.getVitesse(k);
                    double Fi     = p.getForce(k);
                    double Fold_i = Fold[i * dim + k];
                    p.setVitesse(k, vi + dt * 0.5 / mi * (Fi + Fold_i));
                }
            }
        }

        if (step % 1000 == 0) {
            std::cout << "Progression t = " << t << " / " << tend
                      << "   step = " << step << "\r" << std::flush;
            for (Cellule& cell : univers.getCellules()) {
                for (Particule& p : cell) {
                    if (std::isnan(p.getVitesse(0)) || std::isnan(p.getPosition(0))) {
                        std::cout << "\nNaN détecté à step=" << step
                                  << " id=" << p.getId() << "\n";
                        return;
                    }
                }
            }
            sauvegarderVTK("frames/frame_" + std::to_string(step) + ".vtu",
                           univers);
        }

        if (step % 1000 == 0)
            univers.rescaleV(Ec_cible);
    }

    sauvegarderVTK("etat_final.vtu", univers);
    std::cout << "\nSimulation terminée à t = " << t << "\n";
}

int main() {
    try {
        const int    dim   = 2;
        const double Lx    = 250.0;     ///< L1.
        const double Ly    = 180.0;     ///< L2.
        const double rcut  = 2.5;       ///< 2.5 σ.
        const double dt    = 0.0005;    ///< δt.
        const double tend  = 29.5;      ///< t_end.
        const double sigma = 1.0;
        const double d0    = std::pow(2.0, 1.0/6.0) * sigma;

        std::filesystem::create_directory("frames");

        Univers univers(dim, 20000, Lx, Ly, rcut);

        int id = 0;
        const int    M_disque = 23;
        const double R        = 11.21 * d0;
        const double cx       = Lx / 2.0;
        const double cy       = 120.0;      

        for (int j = 0; j < M_disque; ++j) {
            for (int i = 0; i < M_disque; ++i) {
                double x  = cx + i*d0 - (M_disque - 1)/2.0 * d0;
                double y  = cy + j*d0 - (M_disque - 1)/2.0 * d0;
                double dx = x - cx;
                double dy = y - cy;
                if (dx*dx + dy*dy <= R*R) {
                    Particule p({x, y}, {0.0, 10.0}, {0.0, 0.0},
                                1.0, id++, Categorie::Proton);
                    univers.ajouterParticule(p);
                }
            }
        }
        const int N1 = id;
        std::cout << "Disque : " << N1 << " particules placees (cible : 395)\n";

        const int n_col = static_cast<int>(Lx / d0);     
        const int n_lig = 78;                            

        for (int j = 1; j <= n_lig; ++j) {
            for (int i = 0; i < n_col; ++i) {
                Particule p({i*d0, j*d0}, {0.0, 0.0}, {0.0, 0.0},
                            1.001, id++, Categorie::Proton);
                univers.ajouterParticule(p);
            }
        }
        const int N2 = id - N1;
        std::cout << "Membrane : " << N2 << " particules placees (cible : 17227)\n";
        std::cout << "Total : " << id << " particules\n";

        const double Ec_D    = 0.005 * (N1 + N2);
        const int    nb_part = id;

        std::cout << "Ec cible = " << Ec_D << "\n";
        std::cout << "Lancement de la simulation...\n";

        auto start = std::chrono::steady_clock::now();
        Stormer_Verlet(univers, nb_part, dim, tend, dt, Ec_D);
        auto end_t = std::chrono::steady_clock::now();

        std::cout << "Temps ecoule : "
                  << std::chrono::duration<double>(end_t - start).count()
                  << " s\n";
    }
    catch (const std::exception& e) {
        std::cerr << "\n[ ARRET DE LA SIMULATION ]\n"
                  << "Raison : " << e.what() << "\n";
        return 1;
    }

    return 0;
}