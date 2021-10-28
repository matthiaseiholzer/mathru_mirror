use crate::algebra::{
    abstr::{Field, Scalar},
    linear::Matrix,
};
use std::ops::Sub;

impl<T> Sub for Matrix<T> where T: Field + Scalar
{
    type Output = Matrix<T>;

    /// Subtracts two matrices
    ///
    /// A = (a_{ij}) \in T^{m \times n}
    /// B = (b_{ij}) \in T^{m \times n}
    /// A - B = ( a_{ij} - b_{ij} )
    ///
    /// # Arguments
    ///
    /// rhs:
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::Matrix;
    ///
    /// let a: Matrix<f64> = Matrix::new(2, 2, vec![1.0, 0.0, 3.0, -7.0]);
    /// let b: Matrix<f64> = Matrix::new(2, 2, vec![1.0, 0.0, 3.0, -7.0]);
    ///
    /// let c: Matrix<f64> = a - b;
    /// ```
    fn sub(self: Self, rhs: Self) -> Self::Output
    {
        (&self).sub(&rhs)
    }
}

impl<'a, 'b, T> Sub<&'b Matrix<T>> for &'a Matrix<T> where T: Field + Scalar
{
    type Output = Matrix<T>;

    fn sub(self: Self, rhs: &'b Matrix<T>) -> Self::Output
    {
        assert_eq!(self.dim(), rhs.dim());
        let (m, n) = rhs.dim();
        Matrix { m,
                 n,
                 data: self.data
                           .iter()
                           .zip(rhs.data.iter())
                           .map(|(x, y)| *x - *y)
                           .collect::<Vec<T>>() }
    }
}

///
/// Subtracts scalar from all matrix elements
impl<'a, 'b, T> Sub<&'b T> for &'a Matrix<T> where T: Field + Scalar
{
    type Output = Matrix<T>;

    /// Subtracts a scalar value from all matrix elements
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::Matrix;
    ///
    /// let a: Matrix<f64> = Matrix::new(2, 2, vec![1.0, 0.0, 3.0, -7.0]);
    /// let b: Matrix<f64> = &a - &-4.0;
    /// ```
    fn sub(self: Self, rhs: &T) -> Self::Output
    {
        return self.apply(&|x: &T| -> T { *x - *rhs });
    }
}

///
/// Subtracts scalar from all matrix elements
impl<'a, 'b, T> Sub<&'b T> for &'a mut Matrix<T> where T: Field + Scalar
{
    type Output = &'a mut Matrix<T>;

    /// Subtracts a scalar value from all matrix elements
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::Matrix;
    ///
    /// let a: Matrix<f64> = Matrix::new(2, 2, vec![1.0, 0.0, 3.0, -7.0]);
    /// let b: Matrix<f64> = &a - &-4.0;
    /// ```
    fn sub(self: Self, rhs: &T) -> Self::Output
    {
        self.mut_apply(&|x: &mut T| -> T { *x - *rhs });
        return self;
    }
}

impl<T> Sub<T> for Matrix<T> where T: Field + Scalar
{
    type Output = Matrix<T>;

    /// Subtracts a scalar from all matrix elements
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::Matrix;
    ///
    /// let a: Matrix<f64> = Matrix::new(2, 2, vec![1.0, 0.0, 3.0, -7.0]);
    /// let b: Matrix<f64> = a - -4.0;
    /// ```
    fn sub(self: Self, rhs: T) -> Self::Output
    {
        return (&self).sub(&rhs);
    }
}
