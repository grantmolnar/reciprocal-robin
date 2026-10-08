"""Prime-power formulas for the reciprocal sigma family.

The formulas are written once in the active :class:`._field.RealField`, so the
same expressions evaluate either at points or by directed interval arithmetic.
"""

from ._field import current_field
from ._types import NumericInput, Scalar
from ._validation import (
    require_integer,
    require_kappa,
    require_prime,
    require_real_at_least,
    require_real_greater_than,
)


def sigma_minus_one_pp(p: NumericInput, a: NumericInput) -> Scalar:
    r"""Return the prime-power abundancy factor ``sigma_{-1}(p^a)``.

    .. math::
       \sigma_{-1}(p^a) \;=\; \frac{1 - p^{-(a+1)}}{1 - p^{-1}}.

    Parameters
    ----------
    p : real
        A real value greater than 1. Prime inputs give the intended
        prime-power abundancy factor; real inputs are used by the
        auxiliary threshold analysis.
    a : real
        Real exponent ``> -1``. At ``a = 0`` the result is ``1``; the
        extension to ``-1 < a < 0`` is used only in the analytic auxiliary
        function.

    Returns
    -------
    field element
        Computed in the current :class:`._field.RealField`.

    Examples
    --------
    >>> float(sigma_minus_one_pp(2, 2))
    1.75
    """
    require_real_greater_than(p, 1, "p")
    require_real_greater_than(a, -1, "a")
    F = current_field()
    p = F.number(p)
    one = F.number(1)
    return (one - p ** (-(a + 1))) / (one - one / p)


def rksigma_pp(p: NumericInput, a: NumericInput, k: NumericInput) -> Scalar:
    r"""Reciprocal-sigma functional at a prime power.

    .. math::
       \sigma^{[k]}_{-1}(p^a)
       \;=\;
       \sigma_{-1}(p^a)^k - \sigma_{-1}(p^{a-1})^k.

    This is the prime-power value of the Dirichlet convolution
    ``mu * sigma_{-1}^k``. Multiplication over the prime powers of ``n``
    recovers the full multiplicative function.

    Parameters
    ----------
    p : real
        A real value greater than 1; prime inputs give the intended
        prime-power factor.
    a : real
        Exponent ``>= 0``. By convention ``rksigma_pp(p, 0, k) = 1``.
    k : numeric
        Real exponent ``> 1``.

    Examples
    --------
    >>> from mpmath import mpf
    >>> float(rksigma_pp(2, 2, mpf(2)))
    0.8125
    """
    require_real_greater_than(p, 1, "p")
    require_real_at_least(a, 0, "a")
    require_kappa(k)
    F = current_field()
    if a == 0:
        return F.number(1)
    s_a = sigma_minus_one_pp(p, a)
    s_am1 = sigma_minus_one_pp(p, a - 1)
    return s_a**k - s_am1**k


def F_local(p: NumericInput, a: NumericInput, k: NumericInput) -> Scalar:
    r"""Auxiliary function ``F^[k](p, a)`` from Section 4 of the paper.

    .. math::
       F^{[k]}(p, a)
       \;=\;
       \frac{1}{k \log p} \,\log\!
       \left(\frac{\sigma^{[k]}_{-1}(p^a)}{\sigma^{[k]}_{-1}(p^{a-1})}\right).

    The domain is real ``a >= 1`` and real ``k > 1``. At ``a = 1``, the
    convention ``sigma^{[k]}_{-1}(1)=1`` gives
    ``log(sigma_{-1}(p)^k - 1) / (k log p)``.
    """
    require_real_greater_than(p, 1, "p")
    require_real_at_least(a, 1, "a")
    require_kappa(k)
    F = current_field()
    p_field = F.number(p)
    one = F.number(1)
    # Under MpiField, mpi op mpf widens the point to an interval; mpf op mpi
    # attempts the opposite coercion and fails. Keep the field value on the left.
    if a == 1:
        s1 = sigma_minus_one_pp(p, 1)
        return F.log(s1**k - one) / (F.log(p_field) * k)
    num = rksigma_pp(p, a, k)
    den = rksigma_pp(p, a - 1, k)
    return F.log(num / den) / (F.log(p_field) * k)


def r_kappa_pa(p: int, a: int, k: NumericInput) -> Scalar:
    r"""Gap factor ``r_{p,a}^{(\kappa)} = rho^{[k]}(p^a) / sigma_{-1}(p^a)``.

    The logarithmic form used in the paper is:

    .. math::
       \log r_{p,a}^{(\kappa)}
       \;=\;
       \frac{1}{\kappa - 1}\Bigl[
         \log\!\bigl(p^a\, \sigma_{-1}(p^a)\bigr)
         - \log \kappa
         + \log\!\bigl(1 - u_{p,a}^\kappa\bigr)
       \Bigr],

    where ``u_{p,a} = sigma_{-1}(p^{a-1}) / sigma_{-1}(p^a) in (0, 1)``.
    Valid for integer ``a >= 1`` and real ``kappa > 1``.
    """
    require_prime(p)
    require_integer(a, "a", minimum=1)
    require_kappa(k)
    F = current_field()
    p_a_sigma = F.number(p**a) * sigma_minus_one_pp(p, a)
    u = sigma_minus_one_pp(p, a - 1) / sigma_minus_one_pp(p, a)
    log_r = (F.log(p_a_sigma) - F.log(k) + F.log(F.number(1) - u**k)) / (k - 1)
    return F.exp(log_r)
