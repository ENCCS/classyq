#pragma once

#include <Eigen/Core>

#include "cavity/TsLess.hpp"

namespace classyq {
/** Form \f$\mathbf{S}\f$ matrix for a generic environment.
 *
 * @tparam Scalar underlying type of matrix, defaults to `double`.
 * @param cavity the discretization of the surface.
 * @param gf Green's function.
 * @returns The matrix representation, by collocation, of the \f$\mathcal{S}\f$ boundary integral operator.
 */
template<typename Scalar = double>
[[nodiscard]] auto
form_S(const TsLess& cavity, int gf) -> Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

/** Form \f$\mathbf{D}\f$ matrix for a generic environment.
 *
 * @tparam Scalar underlying type of matrix, defaults to `double`.
 * @param cavity the discretization of the surface.
 * @param gf Green's function.
 * @returns The matrix representation, by collocation, of the \f$\mathcal{D}\f$ boundary integral operator.
 */
template<typename Scalar = double>
[[nodiscard]] auto
form_D(const TsLess& cavity, int gf) -> Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

/** Form \f$\mathbf{S}\f$ matrix for vacuum.
 *
 * @param cavity the discretization of the surface.
 * @returns The matrix representation, by collocation, of the \f$\mathcal{S}\f$ boundary integral operator.
 *
 * This function forms the matrix representation, by collocation, of the vacuum \f$\mathcal{S}\f$ boundary integral
 * operator, whose kernel is:
 *
 * \f[
 *   \mathcal{S}f
 *   =
 *   \int \mathrm{d}\mathbf{r}^{\prime}
 *   \frac{f(\mathbf{r}^{\prime})}{| \mathbf{r} - \mathbf{r}^{\prime} |}
 * \f]
 *
 * Given the set of collocation points \f$\lbrace\mathbf{s}_{i}\rbrace\f$, the off-diagonal is straightforward:
 *
 * \f[
 *   S_{ij}
 *   =
 *   \frac{1}{| \mathbf{s}_{i} - \mathbf{s}_{j} |}
 * \f]
 *
 * The diagonal is computed using the formulas reported by \citet{Scalmani2010-tw} for the
 * *self-potential factors*:
 *
 * \f[
 *   f^{\alpha}_{i} = 4\pi -  \sum_{j \neq i} \omega_{\alpha j}\frac{\erf{\zeta_{ij}s_{ij}}}{s_{ij}},
 * \f]
 *
 * these are exact for the centroid collocation of the \f$\mathcal{S}\f$ operator for a single sphere and hence
 * reproduce the Born energy of a ion exactly.
 * Here \f$s_{ij} = |\mathbf{s}_{i} - \mathbf{s}_{j} |\f$ are the distances between collocation points on the unit
 * sphere.
 *
 * Since we use the EQ partition to build the partition of the cavity spheres, the weights are the same for all points
 * and there is no Gaussian smearing of the point charges:
 *
 * \f[
 *   f^{\alpha}_{i} = 4\pi -  \omega_{\alpha}\sum_{j \neq i} \frac{1}{s_{ij}},
 * \f]
 *
 * The self-potential factors for the non-unit sphere give the diagonal are:
 *
 * \f[
 *   S_{ii}(R_{\alpha})
 *   =
 *   R_{\alpha}f^{\alpha}_{i}\
 * \f].
 */
[[nodiscard]] auto
form_S(const TsLess& cavity) -> Eigen::MatrixXd;

/** Form \f$\mathbf{D}\f$ matrix for vacuum.
 *
 * @param cavity the discretization of the surface.
 * @returns The matrix representation, by collocation, of the \f$\mathcal{D}\f$ boundary integral operator.
 *
 * This function forms the matrix representation, by collocation, of the vacuum \f$\mathcal{D}\f$ boundary integral
 * operator, whose kernel is:
 *
 * \f[
 *   \mathcal{D}f
 *   =
 *   \int \mathrm{d}\mathbf{r}^{\prime}
 *   f(\mathbf{r}^{\prime})
 *   \frac{(\mathbf{r} - \mathbf{r}^{\prime})\cdot \mathbf{n}^{\prime}}{| \mathbf{r} - \mathbf{r}^{\prime} |^{3}}
 * \f]
 *
 * Given the set of collocation points \f$\lbrace\mathbf{s}_{i}\rbrace\f$, the off-diagonal is straightforward:
 *
 * \f[
 *   D_{ij}
 *   =
 *   \frac{(\mathbf{s}_{i} - \mathbf{s}_{j})\cdot \mathbf{n}_{j}}{| \mathbf{s}_{i} - \mathbf{s}_{j} |^{3}}
 * \f]
 *
 * The formula is given by \citet{Scalmani2010-tw} for the *self-field factors*:
 *
 * \f[
 *   g^{\alpha}_{i}
 *   =
 *   - 2\pi
 *   - \sum_{j \neq i}
 *      \omega_{\alpha j}
 *      \mathbf{n}_{j} \cdot \frac{\partial }{\partial \mathbf{s}_{j}}
 *      \frac{\erf{\zeta_{ij}s_{ij}}}{s_{ij}}
 * \f]
 *
 * It is derived from Gauss' theorem and was first published by \citet{Purisima1995-kc}
 * These are exact for the centroid collocation of the \f$\mathcal{D}\f$ operator for a single sphere and hence
 * reproduce the Born energy of a ion exactly.
 * Here \f$s_{ij} = |\mathbf{s}_{i} - \mathbf{s}_{j} |\f$ are the distances between collocation points on the unit
 * sphere.
 *
 * Since we use the EQ partition to build the partition of the cavity spheres, the weights are the same for all points
 * and there is no Gaussian smearing of the point charges:
 *
 * \f[
 *   g^{\alpha}_{i}
 *   =
 *   - 2\pi
 *   - \omega_{\alpha}\sum_{j \neq i} \frac{(\mathbf{s}_i - \mathbf{s}_{j})\cdot
 * \hat{\mathbf{n}}_{j}}{s^{3}_{ij}},
 * \f]
 *
 * \f$\hat{\mathbf{n}}_{i}\f$ is the normal vector at EQ point \f$i\f$. Since these are computed
 * for the unit sphere at the origin, the EQ points **are** also the
 * normal vectors.
 *
 * The self-potential factors for the non-unit sphere give the diagonal as:
 *
 * \f[
 *   D_{ii}(R_{\alpha})
 *   =
 *   g^{\alpha}_{i}
 * \f].
 */
[[nodiscard]] auto
form_D(const TsLess& cavity) -> Eigen::MatrixXd;

/** Form \f$\mathbf{R}_{\varepsilon}\f$ matrix.
 *
 * @param cavity the discretization of the surface.
 * @param epsilon permittivity of the uniform dielectric, with \f$\varepsilon > 1\f$.
 * @returns The matrix representation, by collocation, of the \f$\mathcal{R}_{\varepsilon}\f$ boundary integral
 * operator:
 *
 * \f[
 *   \mathbf{R}_{\varepsilon}
 *   =
 *   2\pi\left( \frac{\varepsilon+1}{\varepsilon-1}\right) \mathbf{I} - \mathbf{D}\mathbf{W}
 * \f]
 *
 * @note This is used in the implementation the isotropic IEF solver.
 */
[[nodiscard]] auto
form_R_epsilon(const TsLess& cavity, double epsilon) -> Eigen::MatrixXd;

/** Form \f$\mathbf{R}_{\infty}\f$ matrix.
 *
 * @param cavity the discretization of the surface.
 * @returns The matrix representation, by collocation, of the \f$\mathcal{R}_{\infty}\f$ boundary integral
 * operator:
 *
 * \f[
 *   \mathbf{R}_{\infty}
 *   =
 *   2\pi \mathbf{I} - \mathbf{D}\mathbf{W}
 * \f]
 *
 * @note This is used in the implementation of the isotropic IEF solver.
 */
[[nodiscard]] auto
form_R_infinity(const TsLess& cavity) -> Eigen::MatrixXd;
} // namespace classyq
