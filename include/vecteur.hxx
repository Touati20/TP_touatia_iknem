/**
 * @file vecteur.hxx
 * @brief Déclarations de la classe vecteur (vecteur 3D) et de l'opérateur
 *        de flux associé.
 */

#ifndef VECTEUR_HXX
#define VECTEUR_HXX

#include <ostream>

/**
 * @class vecteur
 * @brief Représente un vecteur dans l'espace euclidien à trois dimensions.
 *
 * La classe encapsule trois composantes réelles (x, y, z) et fournit les
 * opérations vectorielles courantes : addition, soustraction, multiplication
 * par un scalaire, produit scalaire, norme euclidienne et affichage.
 */
class vecteur {

private:
    /** @brief Composantes du vecteur 3D : x, y et z (initialisées à 0.0). */
    double x = 0.0, y = 0.0, z = 0.0;

public:

    /**
     * @brief Constructeur paramétré.
     * @param x Composante selon l'axe x.
     * @param y Composante selon l'axe y.
     * @param z Composante selon l'axe z.
     */
    vecteur(double x, double y, double z);

    /**
     * @brief Constructeur par défaut.
     *
     * Construit un vecteur nul (0, 0, 0).
     */
    vecteur() = default;

    /**
     * @brief Retourne la composante x du vecteur.
     * @return Valeur de x.
     */
    double getX() const;

    /**
     * @brief Retourne la composante y du vecteur.
     * @return Valeur de y.
     */
    double getY() const;

    /**
     * @brief Retourne la composante z du vecteur.
     * @return Valeur de z.
     */
    double getZ() const;

    /**
     * @brief Modifie la composante x du vecteur.
     * @param dx Nouvelle valeur de x.
     */
    void setX(double dx);

    /**
     * @brief Modifie la composante y du vecteur.
     * @param dy Nouvelle valeur de y.
     */
    void setY(double dy);

    /**
     * @brief Modifie la composante z du vecteur.
     * @param dz Nouvelle valeur de z.
     */
    void setZ(double dz);

    /**
     * @brief Calcule la somme du vecteur courant et d'un autre vecteur.
     * @param v Vecteur à additionner.
     * @return Nouveau vecteur résultat de la somme.
     */
    vecteur sum_Vectors(vecteur v);

    /**
     * @brief Calcule la différence entre le vecteur courant et un autre vecteur.
     * @param v Vecteur à soustraire.
     * @return Nouveau vecteur résultat de la différence (this - v).
     */
    vecteur sub_Vectors(vecteur v);

    /**
     * @brief Multiplie le vecteur courant par un scalaire.
     * @param lambda Scalaire multiplicateur.
     * @return Nouveau vecteur dont chaque composante est multipliée par lambda.
     */
    vecteur multV_par_lambda(double lambda);

    /**
     * @brief Calcule le produit scalaire entre le vecteur courant et un autre vecteur.
     * @param v Second vecteur du produit scalaire.
     * @return Valeur du produit scalaire \f$ x \cdot v_x + y \cdot v_y + z \cdot v_z \f$.
     */
    double dot(vecteur v);

    /**
     * @brief Calcule la norme euclidienne du vecteur.
     * @return Norme du vecteur, soit \f$ \sqrt{x^2 + y^2 + z^2} \f$.
     */
    double norme();

    /**
     * @brief Affiche le vecteur sur la sortie standard sous la forme (x, y, z).
     */
    void print();
};

/**
 * @brief Opérateur de flux pour afficher un vecteur sur un flux de sortie.
 *
 * Affiche le vecteur sous la forme "Vecteur(x, y, z)".
 *
 * @param os Flux de sortie.
 * @param v  Vecteur à afficher.
 * @return Référence vers le flux de sortie (pour permettre le chaînage).
 */
std::ostream& operator<<(std::ostream& os, const vecteur& v);

#endif 