use serde::{Deserialize, Serialize};
use std::clone::Clone;
use crate::algebra::linear::{Vector, Matrix};
use crate::algebra::abstr::{Field, Scalar};
use crate::algebra::linear::matrix::{Inverse, Solve, Substitute};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LUDec<T>
{
    l: Matrix<T>,
    u: Matrix<T>,
    p: Matrix<T>
}

impl<T>  LUDec<T>
{
    pub(super) fn new(l: Matrix<T>, u: Matrix<T>, p: Matrix<T>) -> LUDec<T>
    {
        return
            LUDec{
                l: l,
                u: u,
                p: p
            };
    }

    /// Return l Matrix of LU decomposition
    ///
    pub fn l(self: Self) -> Matrix<T>
    {
        return self.l;
    }

    pub fn u(self: Self) -> Matrix<T>
    {
        return self.u;
    }

    pub fn p(self: Self) -> Matrix<T>
    {
        return self.p;
    }

    /// Return l, u, and p matrix of the LU decomposition
    ///
    pub fn lup(self: Self) -> (Matrix<T>, Matrix<T>, Matrix<T>)
    {
        return (self.l, self.u, self.p);
    }
}

impl<T> Solve<Vector<T>> for LUDec<T>
    where T: Field + Scalar
{
    /// Solves Ax = y
    ///  where A \in R^{m * n}, x \in R^n, y \in R^m
    ///
    ///
    fn solve(self: &Self, rhs: &Vector<T>) -> Result<Vector<T>, ()>
    {

        let b_hat: Vector<T> = &self.p * rhs;

        let y: Vector<T> = self.l.substitute_forward(b_hat);

        let x: Vector<T> = self.u.substitute_backward(y);

        return Ok(x);
    }
}

// TOOD
impl<T> Inverse<T> for LUDec<T>
    where T: Field + Scalar
{
    /// Inverse Matrix
    ///
    /// PAX = LUX = I
    /// X = (PA)^-1
    /// X = A^-1P^-1
    /// XP = A^-1
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::{Matrix};
    /// use mathru::algebra::linear::matrix::Inverse;
    ///
    /// let a: Matrix<f64> = Matrix::new(2, 2, vec![1.0, 0.0, 3.0, -7.0]);
    /// let b_inv: Matrix<f64> = a.inv().unwrap();
    ///
    /// ```
    fn inv(self: &Self) -> Result<Matrix<T>, ()>
    {
        assert_eq!(self.p.nrow(), self.p.ncol());
        let b = Matrix::one(self.p.nrow());
        let x: Matrix<T> = self.solve(&b)?;
        return Ok(x)
    }
}

// TOOD
impl<T> Solve<Matrix<T>> for LUDec<T>
    where T: Field + Scalar
{
    fn solve(self: &Self, rhs: &Matrix<T>) -> Result<Matrix<T>, ()>
    {
        let b_hat: Matrix<T> = &self.p * rhs;

        let y: Matrix<T> = self.l.substitute_forward(b_hat);
        let x: Matrix<T> = self.u.substitute_backward(y);

        return Ok(x);
    }
}