#include <gtest/gtest.h>
#include "Particule.hxx" 

// Suite de tests pour la classe Particule


TEST(ParticuleTest, ConstructeurEtGetters) {
    std::vector<double> pos = {1.0, 2.0, 3.0};
    std::vector<double> vit = {0.0, 0.0, 0.0};
    std::vector<double> force = {0.0, 0.0, 0.0};
    
    Particule p(pos, vit, force, 2.5, 42, Categorie::Proton);

    EXPECT_EQ(p.getDim(), 3);
    EXPECT_EQ(p.getMas(), 2.5);
    EXPECT_EQ(p.getId(), 42);
    
    EXPECT_DOUBLE_EQ(p.getPosition(0), 1.0);
    EXPECT_DOUBLE_EQ(p.getPosition(1), 2.0);
    EXPECT_DOUBLE_EQ(p.getPosition(2), 3.0);
}

TEST(ParticuleTest, Setters) {
    std::vector<double> pos = {0.0, 0.0};
    std::vector<double> vit = {0.0, 0.0};
    std::vector<double> force = {0.0, 0.0};
    Particule p(pos, vit, force, 1.0, 1, Categorie::Electron);

    p.setPosition(0, 5.5);
    p.setVitesse(1, -2.3);
    p.setForce(0, 10.0);

    EXPECT_DOUBLE_EQ(p.getPosition(0), 5.5);
    EXPECT_DOUBLE_EQ(p.getVitesse(1), -2.3);
    EXPECT_DOUBLE_EQ(p.getForce(0), 10.0);
}

TEST(ParticuleTest, CalculForceFij) {
    // Particule 1 à l'origine (0, 0)
    std::vector<double> pos1 = {0.0, 0.0};
    std::vector<double> vit1 = {0.0, 0.0};
    std::vector<double> f1 = {0.0, 0.0};
    Particule p1(pos1, vit1, f1, 1.0, 1, Categorie::Proton);

    // Particule 2 en (2, 0) avec une masse de 4.0
    std::vector<double> pos2 = {2.0, 0.0};
    std::vector<double> vit2 = {0.0, 0.0};
    std::vector<double> f2 = {0.0, 0.0};
    Particule p2(pos2, vit2, f2, 4.0, 2, Categorie::Proton);

    // Calcul de la force
    p1.Fij(p2);

    EXPECT_DOUBLE_EQ(p1.getForce(0), 1.0);
    EXPECT_DOUBLE_EQ(p1.getForce(1), 0.0);
    
    EXPECT_DOUBLE_EQ(p2.getForce(0), -1.0);
    EXPECT_DOUBLE_EQ(p2.getForce(1), 0.0);
}

TEST(ParticuleTest, CalculForceFijTropProche) {
    std::vector<double> pos1 = {1.0, 1.0};
    std::vector<double> vit = {0.0, 0.0};
    std::vector<double> f = {0.0, 0.0};
    
    Particule p1(pos1, vit, f, 1.0, 1, Categorie::Proton);
    Particule p2(pos1, vit, f, 1.0, 2, Categorie::Proton);

    p1.Fij(p2);

    EXPECT_DOUBLE_EQ(p1.getForce(0), 0.0);
    EXPECT_DOUBLE_EQ(p2.getForce(0), 0.0);
}