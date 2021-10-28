use crate::algebra::{
    abstr::{Field, Scalar},
    linear::Matrix,
};
use std::ops::Add;

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
    /// ```
    fn add(self: Self, rhs: Self) -> Self::Output
    {
        assert_eq!(self.dim(), rhs.dim());

        let (m, n): (usize, usize) = rhs.dim();

        let mut c: Matrix<T> = rhs;

        T::xaxpy((m * n) as i32, T::one(), &self.data[..], 1, &mut c.data[..], 1);

        return c;
    }
}

///
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
    /// ```
    fn add(self: Self, rhs: &'b Matrix<T>) -> Self::Output
    {
        assert_eq!(self.dim(), rhs.dim());

        let (m, n): (usize, usize) = rhs.dim();

        let mut c: Matrix<T> = rhs.clone();

        T::xaxpy((m * n) as i32, T::one(), &self.data[..], 1, &mut c.data[..], 1);

        return c;
    }
}

impl<T> Matrix<T> where T: Field + Scalar
{
    pub fn add_func(self: &Self, rhs: &Matrix<T>) -> Matrix<T>
    {
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
impl<'a, 'b, T> Add<&'b T> for &'a Matrix<T> where T: Field + Scalar
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
    /// let b: Matrix<f64> = Matrix::new(2, 2, vec![-3.0, -4.0, -1.0, -11.0]);
    /// ```
    fn add(self: Self, rhs: &T) -> Self::Output
    {
        self.apply_mut(&|x: &T| -> T { *x + *rhs });

        return *self;
    }
}

///
/// Add scalar to matrix
impl<T> Add<T> for Matrix<T> where T: Field + Scalar
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
    /// let b: Matrix<f64> = Matrix::new(2, 2, vec![-3.0, -4.0, -1.0, -11.0]);
    ///
    /// ```
    fn add(self: Self, rhs: T) -> Self::Output
    {
        return (&self).add(&rhs);
    }
}
