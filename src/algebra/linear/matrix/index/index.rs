use std::ops::Index;
use crate::algebra::linear::Matrix;

impl<T> Index<[usize; 2]> for Matrix<T>
{
    type Output = T;

    fn index(&self, index: [usize; 2]) -> &Self::Output
    {
        return &self.data[index[1] * self.m + index[0]];
    }
}

impl<T> Matrix<T>
{
    /// Returns the mutual element a_ij from the matrix
    ///
    /// # Example
    ///
    /// ```
    /// # #[macro_use]
        /// # extern crate mathru;
    /// # fn main()
    /// # {
    /// use mathru::algebra::linear::{Matrix};
    ///
    /// let mut a: Matrix<f64> = matrix![1.0, 0.0; 3.0, -7.0];
    /// *a.get_mut(1, 0) = -8.0;
    ///
    /// let a_updated: Matrix<f64> = matrix![1.0, 0.0; -8.0, -7.0];
    /// assert_eq!(a_updated, a);
    /// # }
    /// ```
    pub fn get_mut(self: &mut Self, i: usize, j: usize) -> &mut T
    {
        &mut (self.data[j * self.m + i])
    }

    /// Returns the element a_ij from the matrix
    ///
    /// # Example
    ///
    /// ```
    /// # #[macro_use]
    /// # extern crate mathru;
    /// # fn main()
    /// # {
    /// use mathru::algebra::linear::{Matrix};
    ///
    /// let a: Matrix<f64> = matrix![   1.0, 0.0;
    ///                                 3.0, -7.0];
    ///
    /// let a_ref: f64 = 3.0;
    /// let element: f64 = *a.get(1, 0);
    ///
    /// assert_eq!(a_ref, element);
    /// # }
    /// ```
    pub fn get(self: &Self, i: usize, j: usize) -> &T
    {
        assert!(i < self.m);
        assert!(j < self.n);

        return &self.data[j * self.m + i];
    }
}