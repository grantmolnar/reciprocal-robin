"""Public API shape and input-domain contracts."""

from __future__ import annotations

import pytest

from reciprocal_robin import (
    Factorization,
    MpfField,
    MpiField,
    exceptional_set,
    rho_k,
    robin_ratio,
    using_field,
)
from reciprocal_robin.arithmetic import F_local, rksigma_pp, sigma_minus_one_pp
from reciprocal_robin.exploration import audit_ca_grid, max_robin_ratio_on_ca_grid, reciprocal_ca
from reciprocal_robin.normalized_ca import rho_pp


def test_exploratory_api_has_an_explicit_namespace() -> None:
    import reciprocal_robin as root
    import reciprocal_robin.ca as ca_module
    import reciprocal_robin.exploration as exploration

    exploratory_names = {
        "CAWalkReport",
        "IncrementalRatioWalker",
        "audit_ca_grid",
        "reciprocal_ca",
        "max_robin_ratio_on_ca_grid",
        "walk_ca",
    }
    assert exploratory_names <= set(exploration.__all__)
    assert exploratory_names.isdisjoint(root.__all__)
    assert all(not hasattr(root, name) for name in exploratory_names)

    obsolete_factorization_wrappers = {
        "log_n_from_factorization",
        "rho_k_from_factorization",
        "robin_ratio_from_factorization",
    }
    assert all(not hasattr(root, name) for name in obsolete_factorization_wrappers)
    assert all(not hasattr(ca_module, name) for name in obsolete_factorization_wrappers)
    assert "CAEvent" not in exploration.__all__
    assert not hasattr(exploration, "CAEvent")


def test_private_development_compatibility_module_is_not_shipped() -> None:
    import importlib.util

    assert importlib.util.find_spec("reciprocal_robin.compat") is None


def test_withdrawn_certification_module_is_not_shipped() -> None:
    import importlib.util

    assert importlib.util.find_spec("reciprocal_robin.withdrawn") is None


def test_pointwise_public_domains_are_enforced() -> None:
    with using_field(MpfField()):
        with pytest.raises(ValueError, match="n must be at least 1"):
            rho_k(0, 2)
        with pytest.raises(ValueError, match="kappa must be greater than"):
            rho_k(3, 1)
        with pytest.raises(ValueError, match="n must be at least 3"):
            robin_ratio(2, 2)
        with pytest.raises(ValueError, match="p must be greater than"):
            sigma_minus_one_pp(1, 2)
        with pytest.raises(ValueError, match="a must be greater than"):
            sigma_minus_one_pp(2, -1)
        with pytest.raises(ValueError, match="kappa must be greater than"):
            rksigma_pp(2, 2, 1)
        with pytest.raises(ValueError, match="a must be at least"):
            F_local(2, 0, 2)


def test_factorization_validates_public_mapping() -> None:
    with pytest.raises(ValueError, match="must be prime"):
        Factorization.from_dict({4: 1})
    with pytest.raises(ValueError, match="must be at least 1"):
        Factorization.from_dict({2: 0})
    with pytest.raises(TypeError, match="must be an integer"):
        Factorization.from_dict({2: 1.5})


def test_factorization_private_constructor_cannot_bypass_validation() -> None:
    with pytest.raises(ValueError, match="must be prime"):
        Factorization(_pairs=((4, 1),))
    with pytest.raises(ValueError, match="canonical"):
        Factorization(_pairs=((3, 1), (2, 1)))
    with pytest.raises(ValueError, match="repeated prime"):
        Factorization(_pairs=((2, 1), (2, 2)))


def test_ca_and_normalized_local_domains_are_enforced() -> None:
    with using_field(MpfField()):
        with pytest.raises(ValueError, match="epsilon must be greater"):
            reciprocal_ca(2, "-0.5")
        with pytest.raises(ValueError, match="max_exponent must be at least 1"):
            reciprocal_ca(2, "-0.4", max_exponent=0)
        with pytest.raises(ValueError, match="must be prime"):
            rho_pp(4, 2, 5)


def test_range_and_precision_arguments_are_enforced() -> None:
    with using_field(MpfField()):
        with pytest.raises(ValueError, match="n_max must be at least n_min"):
            exceptional_set(2, n_max=10, n_min=20)
    with pytest.raises(ValueError, match="bumped_dps must be at least base_dps"):
        audit_ca_grid(2, ["-0.4"], base_dps=100, bumped_dps=50)


def test_point_grid_helper_requires_point_field() -> None:
    with using_field(MpiField()):
        with pytest.raises(TypeError, match="point-valued"):
            max_robin_ratio_on_ca_grid(2, ["-0.4"])


def test_factorization_robin_ratio_rejects_n_two() -> None:
    with using_field(MpfField()):
        with pytest.raises(ValueError, match="requires n >= 3"):
            Factorization.from_dict({2: 1}).robin_ratio(2)
