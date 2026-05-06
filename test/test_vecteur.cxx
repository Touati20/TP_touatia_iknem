#include <gtest/gtest.h>
#include "vecteur.hxx" 

TEST(VecteurTest, ConstructeurEtGetters) {
    vecteur v(1.5, -2.0, 3.14);
    EXPECT_DOUBLE_EQ(v.getX(), 1.5);
    EXPECT_DOUBLE_EQ(v.getY(), -2.0);
    EXPECT_DOUBLE_EQ(v.getZ(), 3.14);
}

TEST(VecteurTest, Setters) {
    vecteur v;
    v.setX(10.0);
    v.setY(20.0);
    v.setZ(30.0);
    EXPECT_DOUBLE_EQ(v.getX(), 10.0);
    EXPECT_DOUBLE_EQ(v.getY(), 20.0);
    EXPECT_DOUBLE_EQ(v.getZ(), 30.0);
}

TEST(VecteurTest, OperationsMathematiques) {
    vecteur v1(1.0, 2.0, 3.0);
    vecteur v2(4.0, 5.0, 6.0);

    // Test de l'addition
    vecteur v_sum = v1.sum_Vectors(v2);
    EXPECT_DOUBLE_EQ(v_sum.getX(), 5.0);
    EXPECT_DOUBLE_EQ(v_sum.getY(), 7.0);
    EXPECT_DOUBLE_EQ(v_sum.getZ(), 9.0);

    // Test de la soustraction
    vecteur v_sub = v2.sub_Vectors(v1);
    EXPECT_DOUBLE_EQ(v_sub.getX(), 3.0);
    EXPECT_DOUBLE_EQ(v_sub.getY(), 3.0);
    EXPECT_DOUBLE_EQ(v_sub.getZ(), 3.0);

    vecteur v_mult = v1.multV_par_lambda(2.5);
    EXPECT_DOUBLE_EQ(v_mult.getX(), 2.5);
    EXPECT_DOUBLE_EQ(v_mult.getY(), 5.0);
    EXPECT_DOUBLE_EQ(v_mult.getZ(), 7.5);
}

TEST(VecteurTest, ProduitScalaireEtNorme) {
    vecteur v1(1.0, 2.0, 3.0);
    vecteur v2(4.0, -5.0, 6.0);

    
    EXPECT_DOUBLE_EQ(v1.dot(v2), 12.0);

    vecteur v3(3.0, 4.0, 0.0);
    EXPECT_DOUBLE_EQ(v3.norme(), 5.0);
}