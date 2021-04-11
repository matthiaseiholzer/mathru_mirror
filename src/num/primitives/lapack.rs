use crate::algebra::abstr::{Blas, Lapack, Zero};
use blas;
use lapack;

macro_rules! lapack_impl(
    ($T: ty, $xgehrd: path, $xorghr: path, $xgeev: path, $xgetrf: path, $xgeqrf: path, $xorgqr: path, $xgetri: path, $xpotrf: path,
    $xgetrs: path, $xggev: path)
    => (
        impl Lapack for $T
       	{

       		//Hessenberg
       		fn xgehrd(n: i32, ilo: i32, ihi: i32, a: &mut [Self], lda: i32,
                      tau: &mut [Self], work: &mut [Self], lwork: i32, info: &mut i32)
           	{
                unsafe { $xgehrd(n, ilo, ihi, a, lda, tau, work, lwork, info) }
			}

            fn xgehrd_work_size(n: i32, ilo: i32, ihi: i32, a: &mut [Self], lda: i32,
                                tau: &mut [Self], info: &mut i32) -> i32
            {
                let mut work = [<$T>::zero()];
                let lwork = -1 as i32;

                unsafe { $xgehrd(n, ilo, ihi, a, lda, tau, &mut work, lwork, info) };

                work[0] as i32
            }

         	fn xorghr(n: i32, ilo: i32, ihi: i32, a: &mut [Self], lda: i32, tau: &[Self],
                      work: &mut [Self], lwork: i32, info: &mut i32)
          	{
                unsafe { $xorghr(n, ilo, ihi, a, lda, tau, work, lwork, info) }
            }

            fn xorghr_work_size(n: i32, ilo: i32, ihi: i32, a: &mut [Self], lda: i32,
                                tau: &[Self], info: &mut i32) -> i32 {
                let mut work = [<$T>::zero() ];
                let lwork = -1 as i32;

                unsafe { $xorghr(n, ilo, ihi, a, lda, tau, &mut work, lwork, info) };

                work[0] as i32
            }

            fn xgeev(jobvl: u8, jobvr: u8, n: i32, a: &mut [Self], lda: i32,
                     wr: &mut [Self], wi: &mut [Self],
                     vl: &mut [Self], ldvl: i32, vr: &mut [Self], ldvr: i32,
                     work: &mut [Self], lwork: i32, info: &mut i32)
          	{
                unsafe { $xgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info) }
            }


            fn xgeev_work_size(jobvl: u8, jobvr: u8, n: i32, a: &mut [Self], lda: i32,
                               wr: &mut [Self], wi: &mut [Self], vl: &mut [Self], ldvl: i32,
                               vr: &mut [Self], ldvr: i32, info: &mut i32) -> i32
          	{
                let mut work = [ <$T>::zero() ];
                let lwork = -1 as i32;

                unsafe { $xgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, &mut work, lwork, info) };
                work[0] as i32
			}

			//LU decomposition
			fn xgetrf(m: i32, n: i32, a: &mut [Self], lda: i32, ipiv: &mut [i32], info: &mut i32)
			{
                unsafe { $xgetrf(m, n, a, lda, ipiv, info)}
			}

			fn xgeqrf(m: i32, n: i32, a: &mut [Self], lda: i32, tau: &mut [Self], work: &mut [Self], lwork: i32,
			info: &mut i32)
			{
				unsafe { $xgeqrf(m, n, a, lda, tau, work, lwork, info) };
			}

  			fn xgeqrf_work_size(m: i32, n: i32, a: &mut [Self], lda: i32, tau: &mut [Self], info: &mut i32) -> i32
			{
				let mut work = [<$T>::zero() ];
                let lwork = -1 as i32;

                unsafe { $xgeqrf(m, n, a, lda, tau, &mut work, lwork, info) };
                work[0] as i32
			}

			fn xorgqr(m: i32, n: i32, k: i32, a: &mut [Self], lda: i32, tau: &mut [Self], work: &mut [Self], lwork:
			i32,
			info: &mut i32)
			{
				unsafe { $xorgqr(m, n, k, a, lda, tau, work, lwork, info) };
			}

  			fn xorgqr_work_size(m: i32, n: i32, k: i32, a: &mut [Self], lda: i32, tau: &mut [Self], info: &mut i32) ->
  			 i32
			{
				let mut work = [<$T>::zero() ];
                let lwork = -1 as i32;

                unsafe { $xorgqr(m, n, k, a, lda, tau, &mut work, lwork, info) };
                work[0] as i32
			}

			//
			fn xgetri(n: i32, a: &mut [Self], lda: i32, ipiv: &mut [i32], work: &mut [Self], lwork: i32, info: &mut i32)
			{
                unsafe { $xgetri(n, a, lda, ipiv, work, lwork, info)}
			}

			fn xgetri_work_size(n: i32, a: &mut [Self], lda: i32, ipiv: &mut [i32], info: &mut i32) -> i32
			{
				let mut work = [ <$T>::zero() ];
                let lwork = -1 as i32;
				unsafe { $xgetri(n, a, lda, ipiv, &mut work, lwork, info) };

				work[0] as i32
			}

			//cholsky
			fn xpotrf(uplo: char, n: i32, a: &mut [Self], lda: i32, info: &mut i32)
			{
				unsafe
				{
					$xpotrf(uplo as u8, n, a, lda, info);
				}
			}

			// solve
			fn xgetrs(n: i32, nrhs: i32, a: &mut [Self], lda: i32, ipiv: &mut [i32], b: &mut [Self], ldb: i32, info: &mut i32)
			{
				unsafe
				{
					$xgetrs('N' as u8, n, nrhs, a, lda, ipiv, b, ldb, info);
				}
			}

			// generalized eigenvalue
			fn xggev(n: i32, a: &mut [Self], lda: i32, b: &mut [Self], ldb: i32, alphar: &mut [Self], alphai: &mut [Self], beta: &mut [Self], vl: &mut [Self], ldvl: i32, vr: &mut [Self], ldvr: i32, work: &mut [Self], lwork: i32, info: &mut i32 )
			{
				unsafe
				{
					{ $xggev('V' as u8, 'V' as u8, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work, lwork, info) };
				}
			}

			fn xggev_work_size(n: i32, a: &mut [Self], lda: i32, b: &mut [Self], ldb: i32, alphar: &mut [Self], alphai: &mut [Self], beta: &mut [Self], vl: &mut [Self], ldvl: i32, vr: &mut [Self], ldvr: i32, info: &mut i32 ) -> i32
			{
				let mut work = [ <$T>::zero() ];
                let lwork = -1 as i32;
				unsafe { $xggev('V' as u8, 'V' as u8, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, &mut work, lwork, info) };

				work[0] as i32
			}
      	}
    )
);

macro_rules! blas_impl(
    ($T: ty, $xgemm: path, $xtrsm: path, $xscal: path, $xaxpy: path)
    => (
        impl Blas for $T
       	{
       		fn xgemm(transa: u8, transb: u8, m: i32, n: i32, k: i32, alpha: Self,
    		a: &[Self],
    		lda: i32,
    		b: &[Self],
    		ldb: i32,
    		beta: Self,
    		c: &mut [Self],
    		ldc: i32 )
			{
				unsafe { $xgemm(transa, transb, m, n , k, alpha, a, lda, b, ldb, beta, c, ldc)}
			}

			fn xtrsm(side: char, uplo: char, transa: char, diag: char, m: i32, n: i32, alpha: Self, a: &[Self], lda: i32, b: &mut
			[Self], ldb: i32)
			{
				unsafe
				{
					$xtrsm(side as u8, uplo as u8, transa as u8, diag as u8, m, n, alpha, a, lda, b, ldb);
				}
			}

			fn xscal(n: i32, a: Self, x: &mut [Self], inc: i32)
			{
			    unsafe
			    {
			        $xscal(n, a, x, inc);
			    }
			}

			fn xaxpy(n: i32, a: Self, x: &[Self], incx: i32, y: &mut [Self], incy: i32)
			{
			    unsafe
			    {
			        $xaxpy(n, a, x, incx, y, incy)
			    }
			}
		}
	)
);

lapack_impl!(f32,
             lapack::sgehrd,
             lapack::sorghr,
             lapack::sgeev,
             lapack::sgetrf,
             lapack::sgeqrf,
             lapack::sorgqr,
             lapack::sgetri,
             lapack::spotrf,
             lapack::sgetrs,
             lapack::sggev);

lapack_impl!(f64,
             lapack::dgehrd,
             lapack::dorghr,
             lapack::dgeev,
             lapack::dgetrf,
             lapack::dgeqrf,
             lapack::dorgqr,
             lapack::dgetri,
             lapack::dpotrf,
             lapack::dgetrs,
             lapack::dggev);

blas_impl!(f32, blas::sgemm, blas::strsm, blas::sscal, blas::saxpy);

blas_impl!(f64, blas::dgemm, blas::dtrsm, blas::dscal, blas::daxpy);
