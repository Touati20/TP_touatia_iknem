/**
 * @file test_univers.cxx
 * @brief Tests unitaires pour les classes Cellule, Grille et Univers.
 */

#include <gtest/gtest.h>
#include "Univers.hxx"
#include <cmath>


TEST(CelluleTest, GestionParticules) {
    Cellule c(0.0, 0.0);
    EXPECT_EQ(c.nbParticules(), 0);

    Particule p({1.0, 1.0}, {0.0, 0.0}, {0.0, 0.0}, 1.0, 1, Categorie::Proton);
    c.ajouterParticule(p);
    EXPECT_EQ(c.nbParticules(), 1);

    EXPECT_DOUBLE_EQ(c.getCx(), 0.0);
    EXPECT_DOUBLE_EQ(c.getCy(), 0.0);

    c.vider();
    EXPECT_EQ(c.nbParticules(), 0);
}

TEST(CelluleTest, AjouterEtAccesParticule) {
    Cellule c(2.5, 3.5);

    Particule p1({1.0, 1.0}, {0.0, 0.0}, {0.0, 0.0}, 1.0, 0, Categorie::Proton);
    Particule p2({2.0, 2.0}, {0.5, 0.5}, {0.0, 0.0}, 1.5, 1, Categorie::Electron);
    c.ajouterParticule(p1);
    c.ajouterParticule(p2);

    EXPECT_EQ(c.nbParticules(), 2);
    EXPECT_EQ(c[0].getId(), 0);
    EXPECT_EQ(c[1].getId(), 1);
    EXPECT_DOUBLE_EQ(c[1].getMas(), 1.5);
}

TEST(CelluleTest, GestionVoisins) {
    Cellule c;
    EXPECT_EQ(c.getNbVoisins(), 9);  // valeur par défaut

    c.ajouterVoisin(5);
    c.ajouterVoisin(10);
    EXPECT_EQ(c.getIndicesVoisines().size(), 2);
    EXPECT_EQ(c.getIndicesVoisines()[0], 5);
    EXPECT_EQ(c.getIndicesVoisines()[1], 10);
}


TEST(GrilleTest, ConstructionEtIndexation) {
    Grille g(4, 3, 2.5);
    EXPECT_EQ(g.getNx(), 4);
    EXPECT_EQ(g.getNy(), 3);

    // Indexage row-major : index = x + nx * y
    EXPECT_EQ(g.getIndex(0, 0), 0);
    EXPECT_EQ(g.getIndex(3, 0), 3);
    EXPECT_EQ(g.getIndex(0, 1), 4);
    EXPECT_EQ(g.getIndex(2, 2), 10);
}

TEST(GrilleTest, CentresDesCellules) {
    Grille g(3, 2, 1.0);
    // Cellule (0,0) doit avoir son centre en (0.5, 0.5)
    EXPECT_DOUBLE_EQ(g.getCellule(0, 0).getCx(), 0.5);
    EXPECT_DOUBLE_EQ(g.getCellule(0, 0).getCy(), 0.5);
    // Cellule (2,1) doit avoir son centre en (2.5, 1.5)
    EXPECT_DOUBLE_EQ(g.getCellule(2, 1).getCx(), 2.5);
    EXPECT_DOUBLE_EQ(g.getCellule(2, 1).getCy(), 1.5);
}

TEST(GrilleTest, InitialisationVoisins) {
    Grille g(3, 3, 1.0);
    g.initialiserVoisins();

    // La cellule centrale (1,1) doit avoir 9 voisines (elle-même incluse)
    Cellule& centre = g.getCellule(1, 1);
    EXPECT_EQ(centre.getIndicesVoisines().size(), 9);

    // La cellule en coin (0,0) ne doit avoir que 4 voisines
    Cellule& coin = g.getCellule(0, 0);
    EXPECT_EQ(coin.getIndicesVoisines().size(), 4);

    // Une cellule de bord (1,0) doit avoir 6 voisines
    Cellule& bord = g.getCellule(1, 0);
    EXPECT_EQ(bord.getIndicesVoisines().size(), 6);
}

TEST(UniversTest, MiseAJourVitesse) {
    Univers u(2, 1, 10.0, 10.0, 2.5);
    Particule p({0.0, 0.0}, {0.0, 0.0}, {4.0, 0.0},
                2.0, 0, Categorie::Proton);
    u.ajouterParticule(p);

    // a = F/m = 4/2 = 2, v_new = 0 + 2 * 0.5 = 1
    u.maj_vitesse(0.5);

    // Récupération de la particule mise à jour (elle est dans la cellule (0,0))
    bool trouvee = false;
    for (Cellule& c : u.getCellules()) {
        for (Particule& part : c) {
            EXPECT_DOUBLE_EQ(part.getVitesse(0), 1.0);
            EXPECT_DOUBLE_EQ(part.getVitesse(1), 0.0);
            trouvee = true;
        }
    }
    EXPECT_TRUE(trouvee);
}


TEST(UniversTest, AvancerParticules) {
    Univers u(2, 1, 10.0, 10.0, 2.5);
    Particule p({1.0, 2.0}, {3.0, -1.0}, {0.0, 0.0},
                1.0, 0, Categorie::Proton);
    u.ajouterParticule(p);

    // x_new = x + v*dt = (1 + 3*2, 2 + (-1)*2) = (7, 0)
    u.avancer_parts(2.0);

    bool trouvee = false;
    for (Cellule& c : u.getCellules()) {
        for (Particule& part : c) {
            EXPECT_DOUBLE_EQ(part.getPosition(0), 7.0);
            EXPECT_DOUBLE_EQ(part.getPosition(1), 0.0);
            trouvee = true;
        }
    }
    EXPECT_TRUE(trouvee);
}

TEST(UniversTest, AjouterParticuleEtNbPart) {
    Univers u(2, 0, 10.0, 10.0, 2.5);
    EXPECT_EQ(u.getNbPart(), 0);

    u.ajouterParticule(Particule({1.0, 1.0}, {0.0, 0.0}, {0.0, 0.0},
                                 1.0, 0, Categorie::Proton));
    EXPECT_EQ(u.getNbPart(), 1);

    u.ajouterParticule(Particule({5.0, 5.0}, {0.0, 0.0}, {0.0, 0.0},
                                 1.0, 1, Categorie::Proton));
    EXPECT_EQ(u.getNbPart(), 2);
}

TEST(UniversTest, ReflectionGeometrique) {
    Univers u(2, 1, 10.0, 10.0, 2.5);
    // Particule à 9.5 avec vitesse +1 → après dt=1, à 10.5, hors domaine
    Particule p({9.5, 5.0}, {1.0, 0.0}, {0.0, 0.0},
                1.0, 0, Categorie::Proton);
    u.ajouterParticule(p);

    u.avancer_parts(1.0);
    u.appliquer_reflection();
    u.maj_cellules();

    // Position attendue : 2*L - pos_dépassée = 2*10 - 10.5 = 9.5
    // Vitesse inversée : -1.0
    for (Cellule& c : u.getCellules()) {
        for (Particule& part : c) {
            EXPECT_DOUBLE_EQ(part.getPosition(0), 9.5);
            EXPECT_DOUBLE_EQ(part.getVitesse(0), -1.0);
        }
    }
}

TEST(UniversTest, EnergieCinetique) {
    Univers u(2, 0, 10.0, 10.0, 2.5);

    // Deux particules : Ec = 0.5 * m * v^2 pour chaque
    u.ajouterParticule(Particule({1.0, 1.0}, {3.0, 0.0}, {0.0, 0.0},
                                 2.0, 0, Categorie::Proton));
    u.ajouterParticule(Particule({5.0, 5.0}, {0.0, 4.0}, {0.0, 0.0},
                                 1.0, 1, Categorie::Proton));

    // Ec_total = 0.5*2*9 + 0.5*1*16 = 9 + 8 = 17
    EXPECT_DOUBLE_EQ(u.calculer_Ec(), 17.0);
}

TEST(UniversTest, RescaleV) {
    Univers u(2, 0, 10.0, 10.0, 2.5);
    u.ajouterParticule(Particule({1.0, 1.0}, {2.0, 0.0}, {0.0, 0.0}, 1.0, 0, Categorie::Proton));
    u.rescaleV(8.0);

    for (Cellule& c : u.getCellules()) {
        for (Particule& part : c) {
            EXPECT_DOUBLE_EQ(part.getVitesse(0), 4.0);
        }
    }
    EXPECT_DOUBLE_EQ(u.calculer_Ec(), 8.0);
}

TEST(UniversTest, GraviteAjouteForce) {
    Univers u(2, 0, 10.0, 10.0, 2.5);
    u.ajouterParticule(Particule({5.0, 5.0}, {0.0, 0.0}, {0.0, 0.0}, 2.0, 0, Categorie::Proton));

    u.applique_Potentiel_gravitationnel();

    for (Cellule& c : u.getCellules()) {
        for (Particule& part : c) {
            EXPECT_DOUBLE_EQ(part.getForce(0), 0.0);
            EXPECT_DOUBLE_EQ(part.getForce(1), -24.0);
        }
    }
}