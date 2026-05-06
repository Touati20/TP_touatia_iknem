/**
 * @file Particule.hxx
 * @brief Déclarations de la classe Particule, de l'énumération Categorie
 *        et des fonctions libres associées (calcul de force, intégration
 *        en temps de Störmer–Verlet).
 */

#ifndef PARTICULE_HXX
#define PARTICULE_HXX

#include <vector>

/**
 * @enum Categorie
 * @brief Catégorie physique d'une particule.
 */
enum class Categorie {
    Proton,    ///< Particule de type proton.
    Electron,  ///< Particule de type électron.
    Neutron    ///< Particule de type neutron.
};

/**
 * @class Particule
 * @brief Représente une particule physique dans un espace de dimension quelconque.
 *
 * La classe encapsule la position, la vitesse, la force, la masse,
 * un identifiant et une catégorie. Elle fournit les accesseurs et
 * mutateurs nécessaires, ainsi que la méthode Fij qui calcule
 * l'interaction gravitationnelle entre deux particules en exploitant
 * la troisième loi de Newton.
 */
class Particule {
    private:
        std::vector<double> position;       ///< Vecteur position.
        std::vector<double> vitesse;        ///< Vecteur vitesse.
        std::vector<double> force;          ///< Vecteur force totale.
        double m = 1.0;                     ///< Masse de la particule.
        int Id = 0;                         ///< Identifiant unique.
        Categorie Cat = Categorie::Proton;  ///< Catégorie de la particule.

    public:
        /**
         * @brief Constructeur paramétré.
         * @param pos    Vecteur position initial.
         * @param vitess Vecteur vitesse initial.
         * @param force  Vecteur force initial.
         * @param m      Masse.
         * @param Id     Identifiant unique.
         * @param cat    Catégorie.
         */
        Particule(std::vector<double> pos,
                  std::vector<double> vitess,
                  std::vector<double> force,
                  double m,
                  int Id,
                  Categorie cat);

        /**
         * @brief Constructeur par défaut.
         */
        Particule() = default;

        /**
         * @brief Retourne la dimension de l'espace.
         * @return Nombre de composantes du vecteur position.
         */
        int getDim() const;

        /**
         * @brief Accède à la i-ème composante de la position.
         * @param i Indice de la composante.
         */
        double getPosition(int i) const;

        /**
         * @brief Accède à la i-ème composante de la vitesse.
         * @param i Indice de la composante.
         */
        double getVitesse(int i) const;

        /**
         * @brief Accède à la i-ème composante de la force.
         * @param i Indice de la composante.
         */
        double getForce(int i) const;

        /**
         * @brief Modifie la i-ème composante de la position.
         */
        void setPosition(int i, double val);

        /**
         * @brief Modifie la i-ème composante de la vitesse.
         */
        void setVitesse(int i, double val);

        /**
         * @brief Modifie la i-ème composante de la force.
         */
        void setForce(int i, double val);

        /**
         * @brief Retourne la masse de la particule.
         */
        double getMas() const;

        /**
         * @brief Retourne l'identifiant de la particule.
         */
        int getId() const;

        /**
         * @brief Calcule la force d'interaction gravitationnelle entre la
         *        particule courante et @p p2, et met à jour les forces des deux particules.
         *
         * Exploite la 3ème loi de Newton pour ne calculer qu'une seule fois
         * l'interaction par paire.
         *
         * @param p2 Référence vers la seconde particule.
         */
        void Fij(Particule& p2);
};

/**
 * @brief Calcule la distance euclidienne entre deux particules.
 * @param p1 Première particule.
 * @param p2 Seconde particule.
 * @return Distance euclidienne entre p1 et p2.
 */
double dist(const Particule& p1, const Particule& p2);

/**
 * @brief Met à jour les forces de toutes les particules en parcourant chaque paire.
 *
 * La force de chaque particule est d'abord remise à zéro, puis chaque paire
 * (i, j) avec i < j contribue via Fij.
 *
 * @param particules Vecteur des particules (modifié).
 */
void CalculeF(std::vector<Particule>& particules);

/**
 * @brief Intègre la dynamique du système par le schéma de Störmer–Verlet.
 *
 * Schéma symplectique adapté à l'intégration en temps des équations
 * différentielles issues de la mécanique newtonienne. Met à jour les
 * positions, vitesses et forces des particules de t = 0 à t = tend
 * avec un pas dt fixé.
 *
 * @param particules Vecteur des particules (modifié).
 * @param Fold       Tableau auxiliaire pour stocker les forces du pas précédent.
 */
void Stormer_Verlet(std::vector<Particule>& particules,
                    std::vector<std::vector<double>>& Fold);

#endif // PARTICULE_HPP