# Public reaction records (optional, report-only)

The `chembounce-reaction-records` command imports explicitly selected retained
Open Reaction Database (ORD) records and compares them to an exact product or a
saved final route. It does not change recommendation acceptance, ranking,
scoring, computational support, supplied-observation gates, or any historical
profile. It never assigns synthesis success/failure or experimental validation.

This is a **recorded precedent** interface. Source-byte integrity, chemical
identity, reported reactants, context and measurements remain separate. A
reported compound in a reaction mixture is not a pure isolated product. Missing
or zero measurements are not evidence that a compound is unsynthesizable.

## Runtime and offline example

Import and new matches require Linux, PyArrow and protobuf in an independently
prepared environment (`reaction-records` optional extra). No command installs
packages or retrieves network data. Tested extraction runtime: Python 3.10.12,
PyArrow 23.0.1, protobuf 7.34.1, RDKit 2025.09.6. This is a narrow descriptor
parser, not installation or semantic validation of the complete ORD package.
Default recommendation and report imports do not import these optional packages.

An attributed, unchanged 24,703-byte adverse source fixture is packaged. These
commands use the current environment; choose an interpreter containing the
optional dependencies. `python -m chembounce_synthesis.reaction_record_report`
is equivalent to the installed command.

```bash
chembounce-reaction-records example ./ord-snapshot
chembounce-reaction-records import ./ord-snapshot ./ord-bundle \
  --expected-snapshot-sha256 dbf2dc03a034ba41ece669b2b6b9aa5497fd613198614ddee3e5e3c7f30d574e
chembounce-reaction-records match ./ord-bundle ./product-record-report \
  --product 'C(=C/c1ccccc1)\c1ccccc1' \
  --reactant 'O=S(=O)(Oc1ccccc1)C(F)(F)F' --reactant 'C=Cc1ccccc1' \
  --expected-snapshot-sha256 dbf2dc03a034ba41ece669b2b6b9aa5497fd613198614ddee3e5e3c7f30d574e
```

Open `product-record-report/index.html`; `reaction-records.json` retains all
bindings and typed details. Output directories must be new. The hash above is a
receipt for this packaged snapshot, not a publisher signature or laboratory
certificate. Copying a hash from an untrusted bundle does not create trust.

## Bind evidence to a saved final result and route

```bash
chembounce-reaction-records match ./ord-bundle ./saved-route-record-report \
  --run ./saved-run --route-id ROUTE_ID
```

This reads only bounded `results.json` and `routes.json`. It requires exactly one
selected result/route ID, equal exact final product identities, a connected
acyclic saved graph, bound leaf identities and consistent recorded reactant slots
at every step. Each reaction receives a separate lookup. It saves both file
hashes, route/result hashes, final product and node bindings. It does not replay
reaction rules, load a catalog or a model, generate conformers, inspect other
campaigns, validate scores or requalify the recommendation. Unsupported/missing
saved leaf bindings fail explicitly; no live catalog lookup fills the gap.

Physical context is absent unless supplied with `--context context.json`.
For a product lookup the file contains one context object; for a saved route it
maps reaction node IDs to context objects. `record_context(record, outcome_index)`
returns the exact source protocol view: all reported inputs/preparations,
amounts/order, setup, conditions, notes, workups, outcome time and analyses.
Equality is deliberately strict; a temperature-only object is partial. Even
matching recorded context does not establish that every relevant physical
variable was reported. Unknown roles, quantities, identities or source context
prevent the stronger reported-reaction label.

Python API:

```python
from chembounce_synthesis.reaction_records import import_snapshot, match_product
from chembounce_synthesis.reaction_record_report import match_saved_route, write_report

import_snapshot("ord-snapshot", "ord-bundle")
product_report = match_product("ord-bundle", "C(=C/c1ccccc1)\\c1ccccc1")
route_report = match_saved_route("ord-bundle", "saved-run", "ROUTE_ID")
write_report(route_report, "saved-route-record-report")
```

Every new lookup reloads and rederives source evidence from raw files. There is
no `verify_hashes=False` option and no trusted normalized-object matching API.
`render` can display an already captured JSON report without optional extraction
dependencies, but prominently labels it captured content: rendering provides no
fresh raw-byte, chemistry or saved-artifact verification.

## Snapshot and bundle contract

`ord-retained-snapshot-v1` is illustrated by the packaged `snapshot.json`.
It binds the data provider/revision, the pinned supported schema revision,
dataset basename/hash/size/ID/URL, descriptor basename/hash, explicit
`row_group`, `row`, `reaction_id` selectors and source licensing/attribution.
Paths must be local basenames and regular files; asset symlinks and reserved-name
collisions reject. Acquisition is separate; official-looking URLs are retained
as declared metadata and never authenticated by their spelling.

The bundle preserves original snapshot bytes, dataset, descriptor, each selected
raw protobuf, derived JSON, presence metadata and wire floats. Loading rechecks
dataset/descriptor/record hashes, coordinates, IDs and every normalized record.
An edited identity, role, measurement or provenance projection must equal a fresh
raw derivation or reject. Rewritten coherent source snapshots can pass internal
consistency but cannot acquire an externally trusted receipt automatically.
Importer runtime metadata is explicitly historical reported metadata; new replay
records its own executable, module paths, versions, source hashes and child
resource observations. Full ORD semantic validation is not claimed.

Limits are versioned: 2 MiB dataset bytes, 4,096 rows, 64 selected records,
256 KiB per selected protobuf, 16 MiB declared decompressed column bytes and
actual decoded Arrow buffers. Selected row groups are decoded in single-row
batches; all decoded batches, including rows before the selector, count toward
the decoded-buffer cap. Footer metadata is an early refusal check, **not** an
allocation guarantee. Native decoding runs in a disposable process with at most
256 MiB additional virtual address space after runtime imports, 25 CPU seconds,
35 wall seconds, and a 16 MiB output-file cap. Its observed baseline/cap/peak RSS
are reported. The total address-space cap includes the imported runtime baseline;
it is not a 16 MiB process-memory guarantee. A child crash or resource limit
returns a typed error, never a chemical no-match result. The parent does not
change shared process/environment limits. No downloaded Python is executed.

## Exact identities, missingness and claim boundaries

- Compare full isomeric materials, preserving isotope, charge, disconnected
  components and repeated reactant instances. Only atom maps are cleared in
  comparison copies. No family/scaffold/tautomer matching or fragment selection.
- Unspecified potential stereo, nonabsolute enhanced stereo, unsupported axial
  stereo, conflicting SMILES/MOLBLOCK/CXSMILES or unsupported identity-bearing
  identifiers (including InChI/XYZ/sequence) abstain. Names/accessions are retained
  as unverified descriptions and never resolved to structures.
- Compare explicitly reported REACTANT instances as a multiset; quantities are
  separate. Unknown roles make the complete inventory unresolved. An internal
  standard in the outcomes list is not a product. All outcomes remain separate.
- Retain AREA, YIELD, PURITY, SELECTIVITY and IDENTITY types, value representations,
  units/precision, analysis references, standard and isolation flags. Invalid
  measurements/unresolved analyses produce retained typed record rejections.
  Unknown wire fields/enums, nonfinite values or parsing limits reject the import.
- Preserve explicit false versus absent and zero versus missing when supported by
  protobuf. Implicit scalar defaults and empty repeated/map fields cannot prove
  omission versus an explicitly empty inventory; their metadata states this.
- A complete exact lookup with no retained entry uses
  `no_match_in_retained_snapshot`. If relevant source entries are rejected or
  unassessable, an empty result uses `search_incomplete_no_assessable_match` and
  lists exclusions. Matches can coexist with incomplete source coverage.

The real fixture (`ord-268b5629ab3849f890d911d8df2bb683`) is useful precisely
because it cannot support a synthesis-success claim: the Metal role is unknown,
ligand stereo is unresolved, four products have only names/formulas, and the
analysis is non-isolated GC-FID AREA. Four explicit areas are zero. For
trans-stilbene, the wire float32 expands to `163541584.0` while protobuf JSON uses
the round-trippable short decimal `163541580.0`; both representations and raw
bytes are retained, without inventing measurement precision. There is no reported
yield, DOI or independently inspected analytical attachment.

## Source attribution and scientific scope

Data revision: `704870642e5c93a141e6666436bc14a9b5f880ed` of the
[official ORD repository](https://github.com/open-reaction-database/ord-data).
Schema revision: `51b3a4968835f0090268cab9c0be034b48ac62d4` of
[ord-schema](https://github.com/open-reaction-database/ord-schema).
The small dataset was selected by repository path/size metadata before reading
chemistry, then only its first record was inspected. This is not representative
sampling or measured coverage for generated compounds.

Packaged source data and adapted metadata carry ORD attribution and
CC-BY-SA-4.0 (`LICENSE-DATA`); the schema descriptor carries Apache-2.0
(`LICENSE-SCHEMA`). `NOTICE.txt` distinguishes them and preserves the reported
contributor attribution. These files are included in wheels and copied into
bundles. No rights to unacquired papers/attachments are asserted. Source-backed
records do not prove experimental truth, isolation, route success, activity,
training independence or improved recommendation yield.
