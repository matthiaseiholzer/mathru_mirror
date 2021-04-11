#[cfg(feature = "blaslapack")]
pub mod lapack;


#[cfg(feature = "blaslapack")]
pub mod geneigendec;
#[cfg(feature = "blaslapack")]
pub use self::geneigendec::GenEigenDec;
