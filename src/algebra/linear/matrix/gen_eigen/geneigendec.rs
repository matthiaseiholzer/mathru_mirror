use crate::algebra::linear::{Matrix, Vector};
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};
use std::clone::Clone;

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Debug, Clone)]
pub struct GenEigenDec<T>
{
    value: Vector<T>,
    vector: Matrix<T>,
}

impl<T> GenEigenDec<T>
{
    pub(super) fn new(value: Vector<T>, vector: Matrix<T>) -> GenEigenDec<T>
    {
        return GenEigenDec { value, vector };
    }

    pub fn value(self: Self) -> Vector<T>
    {
        return self.value;
    }

    pub fn vector(self: Self) -> Matrix<T>
    {
        return self.vector;
    }

    pub fn pair(self: Self) -> (Vector<T>, Matrix<T>)
    {
        return (self.value, self.vector);
    }
}
