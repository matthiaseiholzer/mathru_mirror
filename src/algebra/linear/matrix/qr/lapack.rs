use crate::algebra::linear::{Matrix};
use crate::algebra::abstr::{Field, Scalar};
use std::clone::Clone;
use serde::{Deserialize, Serialize};
use crate::elementary::Power;
use crate::algebra::abstr::{Zero};

/// QR decomposition
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QRDec<T>
{
    q: Matrix<T>,
    r: Matrix<T>
}

impl<T> QRDec<T>
{
    pub(self) fn new(q: Matrix<T>, r: Matrix<T>) -> QRDec<T>
    {
        QRDec
        {
            q: q,
            r: r
        }
    }

    /// Return the q matrix of the QR decomposition
    ///
    /// # Arguments
    ///
    /// * `self`
    ///
    pub fn q(self: Self) -> Matrix<T>
    {
        return self.q;
    }

    /// Return the r matrix of the qr decomposition
    ///
    /// # Re
    pub fn r(self: Self) -> Matrix<T>
    {
        return self.r;
    }

    pub fn qr(self: Self) -> (Matrix<T>, Matrix<T>)
    {
        return (self.q, self.r);
    }
}

impl<T> Matrix<T>
    where T: Field + Scalar + Power
{
    /// QR Decomposition with Givens rotations
    ///
    /// A = QR \
    /// Q is an orthogonal matrix \
    /// R is an upper triangular matrix \
    ///
    /// # Panics
    ///
    /// if A is not a square matrix
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::{Matrix};
    ///
    /// let a: Matrix<f64> = Matrix::new(2, 2, vec![1.0, -2.0, 3.0, -7.0]);
    ///
    /// let (q, r): (Matrix<f64>, Matrix<f64>) = a.dec_qr().qr();
    ///
    /// ```
    pub fn dec_qr<'a>(self: &'a Self) -> QRDec<T>
    {
        let (m, n) = self.dim();
        assert!(m >= n);

        return self.dec_qr_r()
    }

    fn dec_qr_r<'a>(self: &'a Self) -> QRDec<T>
    {
        let (m, n) : (usize, usize) = self.dim();

        //lapack(fortran) uses column major order
        let mut self_data = self.clone().data;

        let m_i32: i32 = m as i32;
        let n_i32: i32 = n as i32;
        let m_n_min: usize = m.min(n);

        let mut tau: Vec<T> = vec![Zero::zero(); m_n_min];

        let mut info: i32 = 0;

        let lwork: i32 = T::xgeqrf_work_size(m_i32, n_i32, & mut self_data[..], m_i32, &mut tau[..], &mut info);

        assert_eq!(0, info);

       	let mut work: Vec<T> = vec![T::zero(); lwork as usize];

        T::xgeqrf(
            m_i32,
            n_i32,
            &mut self_data[..],
            m_i32,
            tau.as_mut(),
            &mut work,
            lwork,
            &mut info,
        );

        assert_eq!(0, info);
        let a: Matrix<T> = Matrix::new(m, n, self_data.clone());
        let r: Matrix<T> = a.r();

        let lwork = T::xorgqr_work_size(
            m_i32,
            m_n_min as i32,
            tau.len() as i32,
            &mut self_data[..],
            m_i32,
            &mut tau[..],
            &mut info,
        );

        assert_eq!(0, info);

        let mut work = vec![T::zero(); lwork as usize];

        T::xorgqr(
            m_i32,
            m_n_min as i32,
            tau.len() as i32,
            &mut self_data[..],
            m_i32,
            &mut tau[..],
            &mut work,
            lwork,
            &mut info,
        );
        assert_eq!(0, info);

        let q: Matrix<T> = Matrix::new(m, n, self_data);

        return QRDec::new(q, r);
    }

    fn r(mut self: Self) -> Self
    {
        for i in 1..self.m
        {
            for k in 0..(i.min(self.n))
            {
                *self.get_mut(i, k) = T::zero();
            }
        }

        self
    }
}