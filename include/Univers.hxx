/**
 * @file Univers.hxx
 * @brief Déclarations des classes Cellule, Grille et Univers
 *        pour la simulation de particules par méthode des cellules liées.
 */

#ifndef UNIVERS_HXX
#define UNIVERS_HXX

#include "Particule.hxx"
#include <vector>

/**
 * @class Cellule
 * @brief Cellule spatiale contenant un sous-ensemble de particules,
 *        et connaissant les indices de ses cellules voisines dans la grille.
 */
class Cellule {
private:
    double cx, cy;                                ///< Centre géométrique de la cellule.
    int    nb_voisins = 9;                        ///< Nombre de voisines (9 en 2D).
    std::vector<Particule> particules_internes;   ///< Particules contenues.
    std::vector<int>       indices_voisines;      ///< Indices des cellules voisines.

public:
    Cellule();
    Cellule(double x, double y);

    double getCx() const;
    double getCy() const;
    int    getNbVoisins() const;
    int    nbParticules() const;

    void setNbVoisin(int dim);
    void ajouterVoisin(int index);
    void ajouterParticule(const Particule& p);
    void vider();

    std::vector<Particule>&       getParticules();
    const std::vector<Particule>& getParticules() const;
    const std::vector<int>&       getIndicesVoisines() const;

    Particule& operator[](int i);

    /// Itérateurs (permettent la boucle range-based for sur la cellule).
    auto begin()       { return particules_internes.begin(); }
    auto end()         { return particules_internes.end();   }
    auto begin() const { return particules_internes.begin(); }
    auto end()   const { return particules_internes.end();   }
};

/**
 * @class Grille
 * @brief Grille 2D de cellules pour partitionner l'espace de simulation.
 *
 * Stockage row-major (`index = x + nx * y`).
 * Voisinage de Moore [-1, +1] × [-1, +1] avec gestion correcte des bords.
 */
class Grille {
private:
    int    nx, ny;                  ///< Nombre de cellules selon x et y.
    double taille_cellule;          ///< Côté d'une cellule.
    std::vector<Cellule> cellules;  ///< Tableau linéaire des cellules.

public:
    Grille();
    Grille(int x, int y, double taille);

    int  getNx() const;
    int  getNy() const;
    int  getIndex(int x, int y) const;

    Cellule&                    getCellule(int x, int y);
    std::vector<Cellule>&       getCellules();
    const std::vector<Cellule>& getCellules() const;

    /** @brief Initialise les voisinages de chaque cellule (Moore 3×3 borné). */
    void initialiserVoisins();
};

/**
 * @class Univers
 * @brief Univers de simulation : grille spatiale, forces de Lennard-Jones,
 *        conditions aux bords, énergie cinétique.
 *
 * Le potentiel de Lennard-Jones est défini par
 * \f[
 *   U_{LJ}(r) = 4\,\varepsilon\left[\left(\frac{\sigma}{r}\right)^{12}
 *                                   - \left(\frac{\sigma}{r}\right)^{6}\right].
 * \f]
 *
 * @note La méthode @c all_forces() ne calcule QUE les forces inter-particules.
 *       La gravité et les conditions aux bords doivent être appelées séparément.
 */
class Univers {
private:
    double G       = -12.0;     ///< Accélération gravitationnelle suivant -y.
    int    dim     = 2;         ///< Dimension de l'espace.
    double sigma   = 1.0;       ///< Paramètre σ du potentiel de Lennard-Jones.
    double epsilon = 1.0;       ///< Paramètre ε du potentiel de Lennard-Jones.
    int    nb_particule = 0;    ///< Nombre courant de particules dans l'univers.
    double Lx, Ly;              ///< Dimensions du domaine.
    double rcut;                ///< Rayon de coupure.
    Grille grille;              ///< Partitionnement spatial.

public:
    Univers(int dim, int nb_part, double Lx, double Ly, double rcut);
    Univers() = default;

    // ----- Accès -----
    int  getNbPart() const;
    std::vector<Cellule>&       getCellules();
    const std::vector<Cellule>& getCellules() const;

    // ----- Gestion des particules -----
    void ajouterParticule(const Particule& p);
    void maj_cellules();

    // ----- Intégration en temps -----
    void avancer_parts(double dt);
    void maj_vitesse(double dt);

    // ----- Forces -----
    /** @brief Calcule les forces inter-particules (Lennard-Jones tronqué à rcut). */
    void all_forces();

    /** @brief Ajoute la force gravitationnelle constante \f$ m\,G \f$ suivant y. */
    void applique_Potentiel_gravitationnel();

    double P_Lennard_Jones(double rij) const;
    double P_limites_reflex(double r) const;

    // ----- Conditions aux bords -----
    void appliquer_reflection();
    void appliquer_absorption();
    void appliquer_periodique();
    void potentiel_reflex();

    // ----- Diagnostic -----
    void   univers_state(int step = -1) const;
    double calculer_Ec() const;
    void   rescaleV(double Ec_cible);
};

#endif