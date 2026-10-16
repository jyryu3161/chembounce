"""Synthetic offline fixtures and fake HTTP bytes only; no live data selection."""

import copy
from dataclasses import asdict
import json
from pathlib import Path
import platform
from urllib.parse import parse_qs, urlparse

import pytest
from rdkit import rdBase

from benchmarks.research import qualified_fresh_panel as fresh
from benchmarks.research.protocol import digest, write_json
from chembounce_synthesis.identity import chemical_identity


def target(ident, alias=None):
    number = int(ident[6:])
    return {
        "target_chembl_id": ident,
        "target_type": "SINGLE PROTEIN",
        "organism": "Homo sapiens",
        "species_group_flag": False,
        "pref_name": alias or ident,
        "target_components": [
            {
                "accession": f"FIXTURE{number}",
                "component_id": number,
                "component_type": "PROTEIN",
                "relationship": "SINGLE PROTEIN",
                "component_description": alias or ident,
                "target_component_synonyms": [],
            }
        ],
    }


def activity(ident=1, target_id="CHEMBL1000", smiles="C1CCCCC1", **changes):
    return {
        "activity_id": ident,
        "target_chembl_id": target_id,
        "canonical_smiles": smiles,
        "standard_units": "nM",
        "standard_type": "Ki",
        "standard_flag": 1,
        "assay_type": "B",
        "standard_value": "100",
        "standard_relation": "=",
        "target_organism": "Homo sapiens",
        "assay_chembl_id": "CHEMBL99999",
        "document_year": 2020,
        **changes,
    }


def assay(target_id="CHEMBL1000", ident="CHEMBL99999", **changes):
    return {
        "assay_chembl_id": ident,
        "target_chembl_id": target_id,
        "confidence_score": 9,
        "relationship_type": "D",
        "assay_organism": "Homo sapiens",
        "assay_type": "B",
        "variant_sequence": None,
        **changes,
    }


def empty_exclusions():
    return {
        key: set()
        for key in (
            "exact",
            "families",
            "scaffolds",
            "generic",
            "unresolved_inputs",
            "target_ids",
            "accessions",
            "components",
        )
    } | {"provenance": {"synthetic": True}}


def entries(*records):
    return [
        {"record": record, "source_sha256": "a" * 64, "source_url": fresh.BASE + "assay.json"}
        for record in records
    ]


@pytest.fixture
def protocol(tmp_path):
    root = Path(fresh.__file__).resolve().parents[2]
    identities = tmp_path / "identities.json"
    write_json(identities, [asdict(chemical_identity("Clc1ccccc1"))])
    failed = tmp_path / "normalization-failures.private.json"
    write_json(failed, [{"input_smiles": "invalid", "type": "ValueError", "message": "synthetic"}])
    manifest = tmp_path / "exclusions.json"
    kernels = {
        name: digest(root / name)
        for name in ("chembounce_synthesis/identity.py", "chembounce_synthesis/common.py")
    }
    write_json(
        manifest,
        {
            "final_union": {
                "excluded_target_ids": sorted(fresh.KNOWN_TARGETS),
                "target_aliases": ["ADCY10"],
                "identities_artifact": identities.name,
            },
            "artifacts": {p.name: {"sha256": digest(p)} for p in (identities, failed)},
            "runtime": {"rdkit": rdBase.rdkitVersion},
            "reviewed_identity_kernel": kernels,
            "failed_unique_inputs_by_reason": {"synthetic": 1},
        },
    )
    mappings = tmp_path / "mapping.json"
    write_json(
        mappings,
        {"targets": [target(key) for key in sorted(fresh.KNOWN_TARGETS)] + [target("CHEMBL999", "ADCY10")]},
    )
    required = {
        "benchmarks/research/qualified_fresh_panel.py",
        "benchmarks/research/import_chembl.py",
        "benchmarks/research/assay_metadata.py",
        *kernels,
    }
    return {
        "schema": "qualified-fresh-panel-protocol-v1",
        "policy": copy.deepcopy(fresh.POLICY),
        "runtime": {"python": platform.python_version(), "rdkit": rdBase.rdkitVersion},
        "source_sha256": {name: digest(root / name) for name in required},
        "acquisition_limits": {
            "max_requests": 1287,
            "max_response_bytes": 100_000_000,
            "max_total_bytes": 1_000_000_000,
            "timeout_seconds": 30,
            "wall_seconds": 3600,
        },
        "exclusions": {
            "preparation_manifest": {"path": str(manifest), "sha256": digest(manifest)},
            "target_sources": [
                {"path": str(mappings), "sha256": digest(mappings), "url": fresh.BASE + "target.json"}
            ],
            "remaining_gaps": ["Synthetic inventory; no training independence"],
        },
    }


def test_hash_bound_exclusions_and_required_alias_fail_closed(protocol):
    exclusions = fresh.validate_protocol(protocol)
    assert "CHEMBL999" in exclusions["target_ids"] and "FIXTURE999" in exclusions["accessions"]
    changed = copy.deepcopy(protocol)
    changed["exclusions"]["target_sources"] = []
    with pytest.raises(ValueError, match="Unresolved required"):
        fresh.validate_protocol(changed)
    Path(protocol["exclusions"]["target_sources"][0]["path"]).write_text("{}")
    with pytest.raises(ValueError, match="hash mismatch"):
        fresh.validate_protocol(protocol)


@pytest.mark.parametrize(
    "changes,reason",
    [
        ({"standard_relation": ">", "standard_value": "100"}, "ambiguous"),
        ({"standard_upper_value": "200"}, "range_measurement"),
        ({"data_validity_comment": "bad"}, "validity_flag"),
        ({"document_year": True}, "publication_year"),
        ({"canonical_smiles": "CC(O)c1ccccc1"}, "unspecified"),
        ({"canonical_smiles": "CC=CCc1ccccc1"}, "unspecified"),
        ({"canonical_smiles": "C1CC1.O"}, "disconnected"),
        ({"canonical_smiles": "CCC"}, "no_normalized_ring"),
    ],
)
def test_strict_activity_and_structure_reasons(changes, reason):
    rows, ledger = fresh.curate_target(
        "CHEMBL1000", [activity(**changes)], entries(assay()), empty_exclusions()
    )
    assert rows == [] and reason in ledger[0]["reason"]


def test_specified_stereo_preserved_and_missing_assay_rejected():
    smiles = "C[C@H](O)c1ccccc1"
    rows, _ = fresh.curate_target(
        "CHEMBL1000", [activity(smiles=smiles)], entries(assay()), empty_exclusions()
    )
    assert rows[0]["smiles"] == chemical_identity(smiles).exact_smiles
    rows, ledger = fresh.curate_target("CHEMBL1000", [activity()], [], empty_exclusions())
    assert not rows and ledger[0]["reason"] == "missing_assay_metadata"


@pytest.mark.parametrize(
    "change",
    [
        {"confidence_score": 8},
        {"relationship_type": "H"},
        {"assay_organism": None},
        {"variant_sequence": {}},
        {"target_chembl_id": "CHEMBL1001"},
    ],
)
def test_assay_independent_assignment_not_activity_claim(change):
    rows, ledger = fresh.curate_target(
        "CHEMBL1000", [activity()], entries(assay(**change)), empty_exclusions()
    )
    assert not rows and ledger[0]["reason"].startswith("assay_")


def test_duplicates_conflicts_and_earliest_exact_representative():
    first = activity(1, document_year=2021)
    second = activity(2, document_year=2018)
    rows, ledger = fresh.curate_target(
        "CHEMBL1000", [first, first, second], entries(assay()), empty_exclusions()
    )
    assert rows[0]["activity_id"] == 2 and ledger[1]["reason"] == "identical_duplicate_activity"
    for field, value in (("standard_value", "20000"), ("standard_type", "IC50")):
        bad = second | {field: value}
        rows, ledger = fresh.curate_target("CHEMBL1000", [first, bad], entries(assay()), empty_exclusions())
        assert not rows and all(row["reason"] == "conflicting_family_labels_or_endpoints" for row in ledger)
    with pytest.raises(ValueError, match="Conflicting duplicate activity"):
        fresh.curate_target(
            "CHEMBL1000", [first, first | {"standard_value": "99"}], entries(assay()), empty_exclusions()
        )
    with pytest.raises(ValueError, match="Conflicting assay"):
        fresh.curate_target(
            "CHEMBL1000", [first], entries(assay(), assay(confidence_score=8)), empty_exclusions()
        )


@pytest.mark.parametrize(
    "field,identity_key,reason",
    [
        ("exact", "exact_smiles", "inspected_exact"),
        ("families", "family_smiles", "inspected_chemical"),
        ("scaffolds", "scaffold_family_smiles", "inspected_normalized"),
    ],
)
def test_each_historical_identity_layer_excludes(field, identity_key, reason):
    exclusions = empty_exclusions()
    exclusions[field].add(getattr(chemical_identity("C1CCCCC1"), identity_key))
    rows, ledger = fresh.curate_target("CHEMBL1000", [activity()], entries(assay()), exclusions)
    assert not rows and reason in ledger[0]["reason"]


def test_target_aliases_and_conflicting_mapping_rejected():
    exclusions = empty_exclusions()
    exclusions["accessions"].add("FIXTURE1000")
    one, two = target("CHEMBL1000"), target("CHEMBL1001")
    eligible, ledger = fresh.ordered_targets([one, two], exclusions)
    assert [r["target_id"] for r in eligible] == ["CHEMBL1001"]
    assert "same_protein_alias" in ledger[0]["reason"]
    two["target_components"][0]["accession"] = "FIXTURE1000"
    eligible, ledger = fresh.ordered_targets([one, two], empty_exclusions())
    assert not eligible and all("conflicting_catalog" in row["reason"] for row in ledger)
    with pytest.raises(ValueError, match="Duplicate"):
        fresh.ordered_targets([one, one], empty_exclusions())


def ring(size):
    return "C1" + "C" * (size - 1) + "1"


def test_full_six_by_two_is_ranked_global_and_without_proxy_scores():
    targets = [target(f"CHEMBL{1000 + i}") for i in range(7)]
    calls = []

    def sample(key):
        calls.append(key)
        index = int(key[6:]) - 1000
        return (
            [activity(i + 1, key, ring(3 + index * 2 + i)) for i in range(2)],
            entries(assay(key)),
            {"synthetic": True},
        )

    selected = fresh.select_panel(targets, sample, empty_exclusions())
    assert selected["status"] == "complete" and len(selected["queries"]) == 12 and len(calls) == 6
    expected = sorted(
        [r["target_chembl_id"] for r in targets], key=lambda key: (fresh.rank("fresh-target-v1", key), key)
    )
    assert calls == expected[:6]
    assert len({chemical_identity(q["smiles"]).scaffold_family_smiles for q in selected["queries"]}) == 12
    assert all(set(q) == {"id", "target_id", "smiles"} for q in selected["queries"])
    assert (
        fresh.select_panel(list(reversed(targets)), sample, empty_exclusions())["queries"]
        == selected["queries"]
    )


def test_global_overlap_rejects_whole_target_without_consuming_tentative_scaffold():
    targets = [target(f"CHEMBL{1000 + i}") for i in range(3)]
    order, _ = fresh.ordered_targets(targets, empty_exclusions())
    supplies = {
        order[0]["target_id"]: [ring(3), ring(4)],
        order[1]["target_id"]: [ring(3), ring(5)],
        order[2]["target_id"]: [ring(5), ring(6)],
    }

    def sample(key):
        return [activity(i + 1, key, s) for i, s in enumerate(supplies[key])], entries(assay(key)), {}

    selected = fresh.select_panel(targets, sample, empty_exclusions())
    assert selected["status"] == "infeasible" and selected["queries"] == []
    assert [r["status"] for r in selected["target_decisions"]] == ["selected", "rejected", "selected"]
    assert len(selected["partial_queries_not_authorized"]) == 4


def test_sixty_target_cap_and_thousand_activity_bound():
    targets = [target(f"CHEMBL{1000 + i}") for i in range(61)]
    calls = []

    def sample(key):
        calls.append(key)
        return [], [], {}

    result = fresh.select_panel(targets, sample, empty_exclusions())
    assert (
        len(calls) == 60
        and result["examined_targets"] == 60
        and len(result["eligible_targets_not_examined"]) == 1
    )
    with pytest.raises(ValueError, match="Activity sample exceeds"):
        fresh.curate_target("CHEMBL1000", [activity()] * 1001, [], empty_exclusions())


def page(key, rows, query, total=None):
    limit, offset = int(query["limit"]), int(query["offset"])
    total = len(rows) if total is None else total
    return {
        key: rows,
        "page_meta": {
            "limit": limit,
            "offset": offset,
            "total_count": total,
            "next": fresh.api_url(
                {"targets": "target", "activities": "activity", "assays": "assay"}[key],
                **(query | {"offset": offset + limit}),
            )
            if offset + len(rows) < total
            else None,
        },
    }


def fake_transport(calls, *, release_change=False):
    def fetch(url, timeout, max_bytes):
        calls.append(url)
        parsed = urlparse(url)
        query = {k: v[0] for k, v in parse_qs(parsed.query).items()}
        for key in ("limit", "offset", "standard_flag"):
            if key in query:
                query[key] = int(query[key])
        resource = parsed.path.rsplit("/", 1)[-1]
        if resource == "status.json":
            payload = {
                "status": "UP",
                "chembl_db_version": "ChEMBL_38" if release_change and len(calls) > 1 else fresh.RELEASE,
            }
        elif resource == "target.json":
            payload = page("targets", [target(f"CHEMBL{1000 + i}") for i in range(6)], query)
        elif resource == "activity.json":
            key = query["target_chembl_id"]
            index = int(key[6:]) - 1000
            rows = [
                activity(i + 1, key, ring(3 + index * 2 + i), assay_chembl_id="CHEMBL" + str(9000 + index))
                for i in range(2)
            ]
            payload = page("activities", rows, query)
        else:
            ids = query["assay_chembl_id__in"].split(",")
            payload = page("assays", [assay("CHEMBL" + str(int(key[6:]) - 8000), key) for key in ids], query)
        raw = json.dumps(payload, indent=1).encode() + b"\n"
        return raw, url, 200, {"content-type": "application/json"}

    return fetch


def test_fake_acquisition_exact_bytes_and_offline_replay(protocol, tmp_path):
    calls = []
    snapshot = fresh.acquire(protocol, tmp_path / "raw", transport=fake_transport(calls))
    assert snapshot["status"] == "complete" and len(calls) == 15
    result = fresh.select_offline(protocol, tmp_path / "raw", tmp_path / "selected")
    assert result["status"] == "complete" and len(result["queries"]) == 12
    public = fresh.read(tmp_path / "selected/public-queries.json")
    assert all(set(q) == {"id", "target_id", "smiles"} for q in public["queries"])
    assert (tmp_path / "raw/raw-0000.json").read_bytes().endswith(b"\n")
    assert fresh.read(tmp_path / "raw/requests.json") == snapshot["sources"]


def test_changed_release_preserves_failed_raw_and_no_offline_selection(protocol, tmp_path):
    calls = []
    with pytest.raises(ValueError, match="release/status changed"):
        fresh.acquire(protocol, tmp_path / "raw", transport=fake_transport(calls, release_change=True))
    snapshot = fresh.read(tmp_path / "raw/snapshot.json")
    assert snapshot["status"] == "incomplete" and len(calls) == 15
    with pytest.raises(ValueError, match="Incomplete"):
        fresh.select_offline(protocol, tmp_path / "raw", tmp_path / "selected")


def test_no_retry_when_request_cap_exhausted(protocol, tmp_path):
    protocol["acquisition_limits"]["max_requests"] = 2
    calls = []
    with pytest.raises(ValueError, match="bound exhausted"):
        fresh.acquire(protocol, tmp_path / "raw", transport=fake_transport(calls))
    assert len(calls) == 2 and fresh.read(tmp_path / "raw/snapshot.json")["status"] == "incomplete"


@pytest.mark.parametrize(
    "mutation", ["raw_bytes", "request_parameters", "missing_batch", "release", "extra_sample"]
)
def test_offline_adversarial_artifacts_fail(protocol, tmp_path, mutation):
    fresh.acquire(protocol, tmp_path / "raw", transport=fake_transport([]))
    path = tmp_path / "raw/snapshot.json"
    snapshot = fresh.read(path)
    if mutation == "raw_bytes":
        (tmp_path / "raw/raw-0000.json").write_text("{}")
    elif mutation == "request_parameters":
        source = snapshot["sources"][2]
        source["parameters"]["standard_units"] = "uM"
        source["requested_url"] = fresh.api_url(source["resource"], **source["parameters"])
    elif mutation == "missing_batch":
        next(iter(snapshot["target_samples"].values()))["assay_sources"] = []
    elif mutation == "release":
        snapshot["release_after"] = 0
    else:
        snapshot["target_samples"]["CHEMBL77777"] = next(iter(snapshot["target_samples"].values()))
    write_json(path, snapshot)
    with pytest.raises(ValueError):
        fresh.select_offline(protocol, tmp_path / "raw", tmp_path / "selected")


@pytest.mark.parametrize(
    "mutation", ["next_missing", "wrong_host", "changed_filter", "wrong_offset", "truncated", "changed_total"]
)
def test_pagination_integrity_and_no_extra_activity_page(mutation):
    query = {"limit": 1000, "offset": 0, "order_by": "activity_id", "target_chembl_id": "CHEMBL1000"}
    url = fresh.api_url("activity", **query)
    data = page("activities", [activity(i) for i in range(1000)], query, total=1001)
    total = None
    if mutation == "next_missing":
        data["page_meta"]["next"] = None
    elif mutation == "wrong_host":
        data["page_meta"]["next"] = data["page_meta"]["next"].replace("www.ebi.ac.uk", "example.org")
    elif mutation == "changed_filter":
        data["page_meta"]["next"] += "&standard_units=uM"
    elif mutation == "wrong_offset":
        data["page_meta"]["offset"] = 1
    elif mutation == "truncated":
        data["activities"].pop()
    else:
        total = 1002
    with pytest.raises(ValueError):
        fresh.verify_page(data, "activities", 0, 1000, url, total=total)


def test_missing_alias_alone_and_conflicting_sourced_alias(protocol):
    spec = protocol["exclusions"]["target_sources"][0]
    path = Path(spec["path"])
    source = fresh.read(path)
    write_json(path, {"targets": source["targets"][:-1]})
    spec["sha256"] = digest(path)
    with pytest.raises(ValueError, match="aliases=.*adcy10"):
        fresh.validate_protocol(protocol)
    write_json(path, {"targets": source["targets"] + [target("CHEMBL998", "ADCY10")]})
    spec["sha256"] = digest(path)
    with pytest.raises(ValueError, match="Conflicting required target alias"):
        fresh.validate_protocol(protocol)


def test_failed_http_preserves_exact_body_without_retry(protocol, tmp_path):
    calls = []

    def failure(url, timeout, maximum):
        calls.append(url)
        return b'{"unavailable": true}\n', url, 503, {"retry-after": "10"}

    with pytest.raises(ValueError, match="response status"):
        fresh.acquire(protocol, tmp_path / "raw", transport=failure)
    assert len(calls) == 1
    assert (tmp_path / "raw/raw-0000.json").read_bytes() == b'{"unavailable": true}\n'
    assert fresh.read(tmp_path / "raw/requests.json")[0]["http_status"] == 503


def test_duplicate_json_fields_fail(tmp_path):
    path = tmp_path / "ambiguous.json"
    path.write_text('{"target":1,"target":2}')
    with pytest.raises(ValueError, match="Duplicate JSON"):
        fresh.read(path)
