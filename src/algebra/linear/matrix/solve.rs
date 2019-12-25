use crate::algebra::linear::{Vector, Matrix};
use crate::algebra::abstr::{Real};

#[cfg(feature = "native")]
use crate::algebra::linear::matrix::{Substitute};

#[cfg(feature = "blaslapack")]
use crate::algebra::abstr::{Zero};


pub trait Solve<T>
{
    /// A * x = b
    ///
    ///
    fn solve(self: &Self, rhs: &T) -> Option<T>;
}

impl<T> Solve<Vector<T>> for  Matrix<T>
    where T: Real
{
    /// Solves Ax = y
    ///  where A \in R^{m * n}, x \in R^n, y \in R^m
    ///
    ///
    fn solve(self: &Self, rhs: &Vector<T>) -> Option<Vector<T>>
    {
        return self.solve_vector_r(rhs);
    }

}

impl<T> Solve<Matrix<T>> for Matrix<T>
    where T: Real
{
    fn solve(self: &Self, rhs: &Matrix<T>) -> Option<Matrix<T>>
    {
        return self.solve_matrix_r(rhs);
    }
}

impl<T> Matrix<T>
    where T: Real
{

    #[cfg(feature = "blaslapack")]
    fn solve_vector_r(self: &Self, y: &Vector<T>) -> Option<Vector<T>>
    {
        let (m, n): (usize, usize) = self.dim();
        let m_i32: i32 = m as i32;
        let n_i32: i32 = n as i32;

        let (y_m, _y_n): (usize, usize) = y.dim();
        let y_m_i32: i32 = y_m as i32;

        let dim_min: i32= m_i32.min(n_i32);
        let mut ipiv: Vec<i32> = vec![Zero::zero(); dim_min as usize];

        let mut info: i32 = 0;

        let mut self_data: Vec<T> = self.clone().data;
        let mut y_data: Vec<T> = y.clone().convert_to_vec();

        T::xgetrf(
            m_i32,
            n_i32,
            self_data.as_mut_slice(),
            m_i32,
            ipiv.as_mut_slice(),
            &mut info,
        );

        if info < 0
        {
            return None
        }

        T::xgetrs(
            m_i32,
            1,
            self_data.as_mut_slice(),
            n_i32,
            ipiv.as_mut_slice(),
            y_data.as_mut_slice(),
            y_m_i32,
            &mut info
        );

        if info != 0
        {
            return None
        }

        return Some(Vector::new_column(y_m, y_data));
    }

    #[cfg(feature = "native")]
    fn solve_vector_r(self: &Self, y: &Vector<T>) -> Option<Vector<T>>
    {
        let (l, u, p): (Matrix<T>, Matrix<T>, Matrix<T>) = self.dec_lu().lup();

        let b_hat: Vector<T> = &p * y;

        let y: Vector<T> = l.substitute_forward(b_hat);

        let x: Vector<T> = u.substitute_backward(y);

        return Some(x);
    }

}



impl<T> Matrix<T>
    where T: Real
{
    #[cfg(feature = "blaslapack")]
    pub fn solve_matrix_r(self: &Self, y: &Matrix<T>) -> Option<Matrix<T>>
    {
        let (m, n): (usize, usize) = self.dim();
        let m_i32: i32 = m as i32;
        let n_i32: i32 = n as i32;

        let (y_m, y_n): (usize, usize) = y.dim();
        let y_n_i32: i32 = y_n as i32;

        let dim_min: i32= m_i32.min(n_i32);
        let mut ipiv: Vec<i32> = vec![Zero::zero(); dim_min as usize];

        let mut info: i32 = 0;

        let mut self_data: Vec<T> = self.clone().data;
        let mut y_data: Vec<T> = y.clone().convert_to_vec();

        T::xgetrf(
            m_i32,
            n_i32,
            self_data.as_mut_slice(),
            m_i32,
            ipiv.as_mut_slice(),
            &mut info,
        );

        if info < 0
        {
            return None;
        }

        T::xgetrs(
            n_i32,
            y_n_i32,
            self_data.as_mut_slice(),
            m_i32,
            ipiv.as_mut_slice(),
            y_data.as_mut_slice(),
            y_n_i32,
            &mut info
        );

        if info != 0
        {
            return None
        }

        return Some(Matrix::new(y_m, y_n, y_data));
    }

    #[cfg(feature = "native")]
    pub fn solve_matrix_r(self: &Self, y: &Matrix<T>) -> Option<Matrix<T>>
    {
        return self.dec_lu().solve(y);
    }
}
