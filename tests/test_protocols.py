"""The named research protocol is usable while explicit settings remain authoritative."""

import json

import pytest

import chembounce
import cli
from chembounce_synthesis.options import apply_expansion_options, has_synthesis_overrides
from chembounce_synthesis.protocols import apply_synthesis_protocol, core_protocol, protocol_provenance
from chembounce_synthesis.workflow import _config


def test_research_preset_enables_robust_sampling_without_activity_cutoffs():
    config = _config(
        {
            "schema_version": 1,
            "catalog": "catalog.tsv",
            "rules": "rules.json",
            "provider": {"kind": "saved", "path": "proposals.json"},
            "protocol": "research-v1",
        }
    )
    assert config["allocation_policy"] == "balanced"
    assert config["minimum_scores"] == {}
    refinement = config["refinement"]
    assert refinement["threshold_policy"] == "all_seeds"
    assert refinement["ranking_statistic"] == "minimum"
    assert refinement["adaptive_max_conformers"] > refinement["num_conformers"]
    assert refinement["max_total_conformers"] == 1680
    assert refinement["max_score_range"] == 0.1
    assert refinement["max_metric_range"] == 0.15
    assert "engineering" in protocol_provenance("research-v1")["tolerance_interpretation"]
    assert core_protocol("research-v1")["search_policy"] == "diverse"


def test_explicit_settings_override_protocol_and_application_is_idempotent():
    initial = {
        "protocol": "research-v1",
        "allocation_policy": "encounter",
        "refinement": {"num_conformers": 16, "max_total_conformers": 100},
        "minimum_scores": {"hit_shape": 0.7},
    }
    result = apply_synthesis_protocol(initial)
    assert result["allocation_policy"] == "encounter"
    assert result["refinement"]["num_conformers"] == 16
    assert result["refinement"]["max_total_conformers"] == 100
    assert result["minimum_scores"] == {"hit_shape": 0.7}
    assert apply_synthesis_protocol(result) == result
    assert "top_n" not in initial["refinement"]


def test_protocol_alone_does_not_require_synthesis_config():
    args = cli.argument_parser().parse_args(["-i", "CCc1ccccc1", "-o", "out", "--protocol", "research-v1"])
    assert not has_synthesis_overrides(args)
    assert apply_expansion_options({}, args) == {"protocol": "research-v1"}


def test_cli_research_protocol_changes_core_policy_and_records_it(tmp_path):
    library = tmp_path / "scaffolds.smi"
    library.write_text("c1ccncc1\nC1CCCCC1\n")
    out = tmp_path / "out"
    result = chembounce.main(
        [
            "-i",
            "CCc1ccccc1",
            "-o",
            str(out),
            "--protocol",
            "research-v1",
            "--scaffold-db",
            str(library),
            "--no-fingerprint-db",
            "--no-shape",
            "--overall_max_n",
            "1",
            "--scaffold_top_n",
            "2",
            "--search-candidate-budget",
            "2",
            "--tanimoto_threshold",
            ".99",
        ]
    )
    assert result == 0
    manifest = json.loads((out / "run_manifest.json").read_text())
    assert manifest["protocol"]["name"] == "research-v1"
    assert manifest["parameters"]["search_policy"] == "diverse"
    assert manifest["parameters"]["search_candidate_budget"] == 2


def test_unknown_protocol_fails_explicitly():
    with pytest.raises(ValueError):
        apply_synthesis_protocol({"protocol": "future"})
