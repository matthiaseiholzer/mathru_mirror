use crate::algebra::{
    abstr::{Field, Scalar},
    linear::Matrix,
};
use std::ops::{Add, AddAssign};

impl<T> Add<Self> for Matrix<T> where T: Field + Scalar
{
    type Output = Matrix<T>;

    /// Adds two matrices
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::Matrix;
    ///
    /// let a: Matrix<f64> = Matrix::new(2, 2, vec![1.0, 0.0, 3.0, -7.0]);
    /// let b: Matrix<f64> = Matrix::new(2, 2, vec![1.0, 0.0, 3.0, -7.0]);
    ///
    /// let c: Matrix<f64> = a + b;
    /// ```
    fn add(self: Self, rhs: Self) -> Self::Output
    {
        (&self).add(&rhs)
    }
}

///Adds two matrices
impl<'a, 'b, T> Add<&'b Matrix<T>> for &'a Matrix<T> where T: Field + Scalar
{
    type Output = Matrix<T>;

    /// Adds two matrices
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::Matrix;
    ///
    /// let a: Matrix<f64> = Matrix::new(2, 2, vec![1.0, 0.0, 3.0, -7.0]);
    /// let b: Matrix<f64> = Matrix::new(2, 2, vec![1.0, 0.0, 3.0, -7.0]);
    ///
    /// let c: Matrix<f64> = &b + &a;
    /// ```
    fn add(self: Self, rhs: &'b Matrix<T>) -> Self::Output
    {
        // assert_eq!(self.dim(), rhs.dim());
        let (m, n) = self.dim();
        Matrix { m,
                 n,
                 data: self.data
                           .iter()
                           .zip(rhs.data.iter())
                           .map(|(x, y)| *x + *y)
                           .collect::<Vec<T>>() }
    }
}

///
/// Add scalar to matrix
impl<'a, 'b, T> Add<&'a T> for &'a mut Matrix<T>
    where T: Field + Scalar
{
    type Output = &'a mut Matrix<T>;

    /// Add a scalar to the matrix
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::Matrix;
    ///
    /// let a: Matrix<f64> = Matrix::new(2, 2, vec![1.0, 0.0, 3.0, -7.0]);
    /// let b: Matrix<f64> = &a + &-4.0;
    /// ```
    fn add(self: Self, rhs: &'a T) -> Self::Output
    {
        self.data.iter_mut().for_each(&|x: &mut T| *x += *rhs);
        self
    }
}

///
/// Add scalar to matrix
impl<T> Add<T> for Matrix<T>
    where T: Field + Scalar
{
    type Output = Matrix<T>;

    /// Add a scalar to the matrix
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::Matrix;
    ///
    /// let a: Matrix<f64> = Matrix::new(2, 2, vec![1.0, 0.0, 3.0, -7.0]);
    /// let b: Matrix<f64> = a + -4.0;
    /// ```
    fn add(mut self: Self, rhs: T) -> Self::Output
    {
        (&mut self).add(&rhs);
        return self;
    }
}

// Add scalar to matrix
impl<T> AddAssign<Matrix<T>> for Matrix<T>
    where T: Field + Scalar
{
    /// Add a scalar to the matrix
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::Matrix;
    ///
    /// let mut a: Matrix<f64> = Matrix::new(2, 2, vec![1.0, 0.0, 3.0, -7.0]);
    /// let b: Matrix<f64> = Matrix::new(2, 2, vec![2.0, 3.0, -5.0, 2.0]);
    /// a += b;
    /// ```
    fn add_assign(&mut self, rhs: Matrix<T>)
    {
        self.data.iter_mut().zip(rhs.data.iter()).for_each(|(a, b)|
            *a += *b
        )
    }
}

// Add scalar to matrix
impl<T> AddAssign<T> for Matrix<T>
    where T: Field + Scalar
{
    /// Add a scalar to the matrix
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::Matrix;
    ///
    /// let mut a: Matrix<f64> = Matrix::new(2, 2, vec![1.0, 0.0, 3.0, -7.0]);
    /// a += -4.0;
    /// ```
    fn add_assign(&mut self, rhs: T)
    {
        self.data.iter_mut().for_each(|a: &mut T|
            *a += rhs
        );
    }
}

