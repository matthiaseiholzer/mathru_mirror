use std::ops::{Index, IndexMut};
use crate::algebra::linear::Vector;

impl<T> Index<usize> for Vector<T>
{
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output
    {
        let (m, n): (usize, usize) = self.data.dim();

        if m == 1
        {
            return self.data.get(0, index)
        }

        if n == 1
        {
            return self.data.get(index, 0)
        }
        panic!("")
    }
}

impl<T> Vector<T>
{
    /// Returns the component
    ///
    /// # Panics
    ///
    /// if i out of bounds
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::Vector;
    ///
    /// let mut a: Vector<f64> = Vector::new_row(vec![1.0, 0.0, 3.0, -2.0]);
    ///
    /// assert_eq!(-2.0, *a.get_mut(3))
    /// ```
    pub fn get(self: &Self, i: usize) -> &T
    {
        let (m, n): (usize, usize) = self.data.dim();

        if m == 1
        {
            return self.data.get(0, i)
        }

        if n == 1
        {
            return self.data.get(i, 0)
        }
        panic!("")
    }
}

impl<T> Vector<T>
//where T: One + Zero
{
    /// Returns the mutual component
    ///
    /// # Arguments
    ///
    /// # Panics
    ///
    /// if i out of bounds
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::Vector;
    ///
    /// let mut a: Vector<f64> = Vector::new_row(vec![1.0, 0.0, 3.0, -2.0]);
    ///
    /// *a.get_mut(1) = -4.0;
    /// ```
    pub fn get_mut(self: &mut Self, i: usize) -> &mut T
    {
        let (m, n): (usize, usize) = self.data.dim();
        assert!(m == 1 || n == 1);

        if m == 1
        {
            //row vector
            // assert!(i < n);
            self.data.get_mut(0, i)
        }
        else
        {
            //column vector
            // assert!(i < m);
            self.data.get_mut(i, 0)
        }
    }
}

impl<T> IndexMut<usize> for Vector<T>
{
    fn index_mut(self: &mut Self, index: usize) -> &mut Self::Output
    {
        let (m, n): (usize, usize) = self.data.dim();

        if m == 1
        {
            return self.data.get_mut(0, index)
        }

        if n == 1
        {
            return self.data.get_mut(index, 0)
        }
        panic!("")
    }
}

