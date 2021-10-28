#[macro_use]
pub mod vector;
pub mod vectorintoiterator;
pub mod vectoriterator;
pub mod vectoriteratormut;

mod mul;
mod add;
mod sub;
mod div;
mod index;

pub use self::{
    vector::Vector, vectorintoiterator::VectorIntoIterator, vectoriterator::VectorIterator,
    vectoriteratormut::VectorIteratorMut,
};
