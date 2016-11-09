namespace netlib
{

/// Iterative template routine -- BiCGSTAB
///
/// Original version from [netlib](http://www.netlib.org/templates/cpp) 
/// as of 2016-11-08 (file date: 1998-07-21), modified
/// for numcxx.
///
/// BiCGSTAB solves the unsymmetric linear system Ax = b 
/// using the Preconditioned BiConjugate Gradient Stabilized method
///
/// BiCGSTAB follows the algorithm described on p. 27 of the 
/// SIAM Templates book.
///
/// The return value indicates convergence within max_iter (input)
/// iterations (0), or no convergence within max_iter iterations (1).
///
/// Upon successful return, output arguments have the following values:
///  
/// \param       x  --  approximate solution to Ax = b
/// \param max_iter  --  the number of iterations performed before the
///               tolerance was reached
/// \param     tol  --  the residual after the final iteration


    template < class Matrix, class Vector, class Preconditioner, class Real >
    int 
    BiCGSTAB(const Matrix &A, Vector &x, const Vector &b,
             const Preconditioner &M, int &max_iter, Real &tol)
    {
        Real resid;
        Vector rho_1(1), rho_2(1), alpha(1), beta(1), omega(1);
        Vector p, phat, s, shat, t, v;

        Real normb = norm(b);
        Vector r = b - A * x;
        Vector rtilde = r;

        if (normb == 0.0)
            normb = 1;
  
        if ((resid = norm(r) / normb) <= tol) {
            tol = resid;
            max_iter = 0;
            return 0;
        }

        for (int i = 1; i <= max_iter; i++) {
            rho_1(0) = dot(rtilde, r);
            if (rho_1(0) == 0) {
                tol = norm(r) / normb;
                return 2;
            }
            if (i == 1)
                p = r;
            else {
                beta(0) = (rho_1(0)/rho_2(0)) * (alpha(0)/omega(0));
                p = r + beta(0) * (p - omega(0) * v);
            }
            // numcxx: was
            // phat = M.solve(p);
            M.solve(phat,p);
            v = A * phat;
            alpha(0) = rho_1(0) / dot(rtilde, v);
            s = r - alpha(0) * v;
            if ((resid = norm(s)/normb) < tol) {
                x += alpha(0) * phat;
                tol = resid;
                return 0;
            }
            // numcxx: was
            //shat = M.solve(s);
            M.solve(shat,s);
            t = A * shat;
            omega = dot(t,s) / dot(t,t);
            x += alpha(0) * phat + omega(0) * shat;
            r = s - omega(0) * t;

            rho_2(0) = rho_1(0);
            if ((resid = norm(r) / normb) < tol) {
                tol = resid;
                max_iter = i;
                return 0;
            }
            if (omega(0) == 0) {
                tol = norm(r) / normb;
                return 3;
            }
        }

        tol = resid;
        return 1;
    }
}
