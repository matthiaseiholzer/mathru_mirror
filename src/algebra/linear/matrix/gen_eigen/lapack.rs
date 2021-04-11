use crate::{
    algebra::{
        abstr::{Field, Scalar},
        linear::{matrix::GenEigenDec, Matrix, Vector},
    },
    elementary::Power,
};

impl<T> Matrix<T> where T: Field + Scalar + Power
{
    /// Computes for a pair of N-by-N real nonsymmetric matrices (A,B) the generalized eigenvalues and the left and right
    /// generalized eigenvectors.
    ///
    /// # Arguments
    ///
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::{matrix::GenEigenDec, Matrix, Vector};
    ///
    /// let a: Matrix<f64> = Matrix::new(3, 3, vec![1.0, -3.0, 3.0, 3.0, -5.0, 3.0, 6.0, -6.0, 4.0]);
    /// let b: Matrix<f64> = Matrix::new(3, 3, vec![1.0, -3.0, 3.0, 3.0, -5.0, 3.0, 6.0, -6.0, 4.0]);
    /// //let eigen: GenEigenDec<f64> = a.dec_gen_eigen(b);
    /// ```
    pub fn dec_gen_eigen(self: Self, b: Self) //-> GenEigenDec<T>
    {

        // return GenEigenDec::new(Vector::new_column(m, wr), Matrix::new(n, n, temp2));
    }
}
