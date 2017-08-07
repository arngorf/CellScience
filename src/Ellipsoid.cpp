#include "Ellipsoid.hpp"

Ellipsoid::Ellipsoid(Vector const &center, Vector const &radii, Mat const &evecs, Vector const &v)
                    : center(center), radii(radii), R(evecs.transpose())
                    , D(radii.asDiagonal()), v(v) {
    setDiagonalMatrix(radii);
    sortDiagonalEntries();
}

Ellipsoid::Ellipsoid(Vector const &center, Vector const &radii, Mat const &R)
                    : center(center), radii(radii), R(R)
                    , D(radii.asDiagonal()), v(Vector::Zero(10)) {
    setDiagonalMatrix(radii);
    sortDiagonalEntries();
    // Missing calculation of v for completion.
}

void Ellipsoid::setDiagonalMatrix(Vector const &radii) {
    D = Mat::Zero(3,3);
    for (int i = 0; i < 3; ++i) {
        D(i,i) = 1.0/std::pow(radii(i), 2);
    }
}

void Ellipsoid::sortDiagonalEntries() {
    for (int i = 3; i > 1; --i) {
        for (int j = 0; j < i - 1; ++j) {
            if (radii(j+1) > radii(j)) {
                R.row(j).swap(R.row(j+1));
                std::swap(radii(j),radii(j+1));
            }
        }
    }

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (R(i,j) > 0) {
                break;
            } else if (R(i,j) < 0) {
                R.row(i) = -R.row(i);
                break;
            }
        }
    }
}