"""Bounded, standard-library-only validation of optional producer display evidence.

These are saved-metadata checks, not chemistry, image semantics or authenticity.
"""

import base64
import hashlib
import json
import math
import struct
import zlib

VERSION = "saved-structure-pose-v1"
ARTIFACT = "review_visuals.json"
LIMITS = {
    "allowance_seconds": 5.0,
    "candidates": 32,
    "seeds_per_candidate": 32,
    "atoms": 512,
    "bonds": 1024,
    "png_bytes": 1048576,
    "image_width": 1024,
    "image_height": 1024,
    "image_pixels": 1048576,
    "artifact_bytes": 16777216,
    "coordinate_magnitude": 1000000.0,
    "zoom_min": 0.25,
    "zoom_max": 8.0,
}
BOND_TYPES = {"SINGLE", "DOUBLE", "TRIPLE", "AROMATIC", "DATIVE", "ZERO", "QUADRUPLE"}
ELEMENTS = (
    "* H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"
).split()


def binding(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def integer(value, lower, upper):
    return type(value) is int and lower <= value <= upper


def number(value, limit):
    return type(value) in (int, float) and math.isfinite(value) and abs(value) <= limit


def digest(value):
    return isinstance(value, str) and len(value) == 64 and all(c in "0123456789abcdef" for c in value)


def checked_png(data):
    """Static 8-bit RGB/RGBA PNG only; CRC/order/decompressed extent bounded."""
    require(isinstance(data, bytes) and 0 < len(data) <= LIMITS["png_bytes"], "PNG byte limit")
    require(data[:8] == b"\x89PNG\r\n\x1a\n", "PNG signature")
    offset, types, compressed, ended = 8, [], bytearray(), False
    width = height = color = None
    while offset < len(data):
        require(offset + 12 <= len(data), "PNG chunk header")
        size = struct.unpack(">I", data[offset : offset + 4])[0]
        kind = data[offset + 4 : offset + 8]
        require(size <= LIMITS["png_bytes"] and offset + 12 + size <= len(data), "PNG chunk length")
        payload = data[offset + 8 : offset + 8 + size]
        crc = struct.unpack(">I", data[offset + 8 + size : offset + 12 + size])[0]
        require(zlib.crc32(kind + payload) & 0xFFFFFFFF == crc, "PNG CRC")
        require(kind not in (b"acTL", b"fcTL", b"fdAT"), "Animated PNG unsupported")
        require(not ended and (types or kind == b"IHDR"), "PNG chunk order")
        if kind == b"IHDR":
            require(not types and size == 13, "PNG IHDR")
            width, height, depth, color, compression, filtering, interlace = struct.unpack(
                ">IIBBBBB", payload
            )
            require(
                0 < width <= LIMITS["image_width"]
                and 0 < height <= LIMITS["image_height"]
                and width * height <= LIMITS["image_pixels"],
                "PNG dimensions",
            )
            require(
                depth == 8 and color in (2, 6) and (compression, filtering, interlace) == (0, 0, 0),
                "Unsupported PNG encoding",
            )
        elif kind == b"IDAT":
            require(b"IDAT" not in types or types[-1] == b"IDAT", "Noncontiguous PNG image data")
            compressed.extend(payload)
        elif kind == b"IEND":
            require(size == 0 and b"IDAT" in types, "PNG IEND")
            ended = True
        else:
            # No compressed metadata/profile chunks: they can expand outside the
            # bounded image raster. The producer disables RDKit zTXt metadata.
            allowed = {b"bKGD": 6, b"pHYs": 9, b"sRGB": 1, b"gAMA": 4, b"cHRM": 32}
            require(
                kind in allowed and size == allowed[kind] and kind not in types and b"IDAT" not in types,
                "Unsupported PNG ancillary/critical chunk",
            )
        types.append(kind)
        offset += size + 12
    require(ended and types[-1] == b"IEND", "PNG end/trailing data")
    expected = height * (1 + width * (3 if color == 2 else 4))
    decoder = zlib.decompressobj()
    try:
        pixels = decoder.decompress(bytes(compressed), expected + 1)
    except zlib.error as exc:
        raise ValueError("Invalid PNG compressed image data") from exc
    require(
        len(pixels) == expected and decoder.eof and not decoder.unused_data and not decoder.unconsumed_tail,
        "PNG decompressed size",
    )
    stride = expected // height
    require(all(pixels[i] <= 4 for i in range(0, expected, stride)), "PNG filter")
    return width, height


def checked_graph(graph, pose=None):
    require(isinstance(graph, dict), "Missing indexed graph")
    require(
        digest(graph.get("sha256"))
        and binding({k: v for k, v in graph.items() if k != "sha256"}) == graph["sha256"],
        "Graph binding",
    )
    atoms, bonds = graph["atoms"], graph["bonds"]
    require(
        isinstance(atoms, list)
        and 0 < len(atoms) <= LIMITS["atoms"]
        and isinstance(bonds, list)
        and len(bonds) <= LIMITS["bonds"],
        "Graph limits",
    )
    for index, atom in enumerate(atoms):
        require(
            type(atom) is dict and atom.get("index") == index and type(atom.get("index")) is int,
            "Atom indices",
        )
        z = atom["atomic_number"]
        require(integer(z, 1, 118) and atom["element"] == ELEMENTS[z], "Element/atomic number")
        require(
            integer(atom["isotope"], 0, 1000) and integer(atom["formal_charge"], -32, 32),
            "Atom isotope/charge",
        )
        require(
            atom["chiral_tag"] in ("CHI_UNSPECIFIED", "CHI_TETRAHEDRAL_CW", "CHI_TETRAHEDRAL_CCW")
            and atom["cip"] in (None, "R", "S", "r", "s"),
            "Atom stereo label",
        )
    pairs = set()
    for bond in bonds:
        a, b = bond["begin"], bond["end"]
        require(
            integer(a, 0, len(atoms) - 1) and integer(b, 0, len(atoms) - 1) and a < b and (a, b) not in pairs,
            "Bond endpoints/duplicates",
        )
        pairs.add((a, b))
        require(
            bond["type"] in BOND_TYPES
            and bond["stereo"]
            in ("STEREONONE", "STEREOANY", "STEREOZ", "STEREOE", "STEREOCIS", "STEREOTRANS"),
            "Bond type/stereo",
        )
        stereo = bond["stereo_atoms"]
        require(
            isinstance(stereo, list)
            and len(stereo) in (0, 2)
            and all(integer(x, 0, len(atoms) - 1) for x in stereo),
            "Bond stereo atoms",
        )
    topology = binding(
        {
            "atoms": [(a["atomic_number"], a["isotope"], a["formal_charge"]) for a in atoms],
            "bonds": sorted((b["begin"], b["end"], b["type"]) for b in bonds),
        }
    )
    require(graph["atom_order_sha256"] == topology, "Indexed topology projection")
    if pose is not None:
        require(
            isinstance(pose, dict)
            and binding({k: v for k, v in pose.items() if k != "sha256"}) == pose["sha256"],
            "Saved pose binding",
        )
        require(
            graph["smiles"] == pose["smiles"] and topology == pose["atom_order_sha256"], "Pose graph identity"
        )
        xyz = pose["coordinates"]
        require(isinstance(xyz, list) and len(xyz) == len(atoms), "Pose coordinate count")
        require(
            all(
                isinstance(p, list)
                and len(p) == 3
                and all(number(x, LIMITS["coordinate_magnitude"]) for x in p)
                for p in xyz
            ),
            "Pose coordinate bounds",
        )
    return graph


def _unavailable(item):
    require(
        item.get("status") == "unavailable"
        and isinstance(item.get("reason"), str)
        and 0 < len(item["reason"]) <= 1024,
        "Unavailable reason",
    )
    require(set(item) == {"status", "reason"}, "Unavailable payload")


def _image(item, smiles):
    if item.get("status") == "unavailable":
        _unavailable(item)
        return
    require(item.get("status") == "available" and item["smiles"] == smiles, "Depiction identity")
    checked_graph(item["graph"])
    require(item["graph"]["smiles"] == smiles, "Depiction graph identity")
    encoded = item["png_base64"]
    require(
        isinstance(encoded, str) and len(encoded) <= (LIMITS["png_bytes"] + 2) // 3 * 4, "PNG base64 limit"
    )
    data = base64.b64decode(encoded, validate=True)
    require(hashlib.sha256(data).hexdigest() == item["png_sha256"], "PNG digest")
    require(checked_png(data) == (item["width"], item["height"]), "PNG declared dimensions")


def validate(artifact, descriptor, inventory, rows, refinements, hit):
    """Return independent candidate display views, or raise a metadata diagnostic."""
    require(
        descriptor["version"] == VERSION and descriptor["status"] in ("available", "unavailable"),
        "Display descriptor",
    )
    if descriptor["status"] == "unavailable":
        require(isinstance(descriptor["reason"], str) and descriptor["reason"], "Display failure reason")
        require("sha256" not in descriptor and "path" not in descriptor, "Unavailable artifact advertised")
        return [{"status": "unavailable", "reason": descriptor["reason"]} for _ in rows]
    require(
        descriptor["publication_completion_within_budget"] is None
        and descriptor["timing_scope"] == "final_persistence_unverified"
        and number(descriptor["prepublication_elapsed_seconds"], LIMITS["allowance_seconds"])
        and descriptor["prepublication_elapsed_seconds"] >= 0,
        "Display publication timing is not verified",
    )
    require(
        descriptor["path"] == ARTIFACT
        and inventory["sha256"] == descriptor["sha256"]
        and inventory["bytes"] <= LIMITS["artifact_bytes"],
        "Display file binding",
    )
    require(
        artifact["version"] == VERSION
        and artifact["schema_version"] == 1
        and type(artifact["schema_version"]) is int
        and artifact["limits"] == LIMITS
        and artifact["hit"] == hit,
        "Display schema/limits/hit",
    )
    require(
        isinstance(artifact["producer"], dict)
        and digest(artifact["producer"]["source_sha256"])
        and isinstance(artifact["producer"]["rdkit_version"], str),
        "Display producer",
    )
    require(
        number(artifact["generation_elapsed_seconds"], LIMITS["allowance_seconds"])
        and artifact["generation_elapsed_seconds"] >= 0,
        "Display elapsed",
    )
    entries = artifact["candidates"]
    require(
        isinstance(entries, list) and len(entries) == len(rows) <= LIMITS["candidates"],
        "Display candidate coverage",
    )
    result = []
    keys = set()
    for entry, row in zip(entries, rows):
        key = (row["route_id"], row["smiles"], row["design"])
        require(key not in keys, "Duplicate display candidate")
        keys.add(key)
        require(
            entry["row_sha256"] == binding(row) and entry["key"] == list(key), "Display selected row binding"
        )
        matches = [
            r for r in refinements if r.get("product") == row["smiles"] and r.get("design") == row["design"]
        ]
        require(len(matches) <= 1, "Ambiguous display refinement")
        trials = matches[0]["trials"] if matches else []
        require(
            isinstance(trials, list)
            and len(trials) <= LIMITS["seeds_per_candidate"]
            and len(entry["seeds"]) == len(trials),
            "Display retained seed coverage",
        )
        require(entry["representative_seed"] == row.get("pose_seed"), "Display representative seed")
        require(set(entry["depictions"]) == {"original_hit", "final_product"}, "Depiction coverage")
        for role, smiles in (("original_hit", hit), ("final_product", row["smiles"])):
            _image(entry["depictions"][role], smiles)
        seen = set()
        seed_views = []
        for saved, trial in zip(entry["seeds"], trials):
            seed = trial["seed"]
            require(integer(seed, 0, 2147483647) and seed not in seen, "Duplicate/invalid display seed")
            seen.add(seed)
            require(
                type(saved["seed"]) is int
                and saved["seed"] == seed
                and saved["trial_sha256"] == binding(trial),
                "Display trial binding",
            )
            pair = saved["pair"]
            view = dict(saved)
            if pair.get("status") == "unavailable":
                _unavailable(pair)
            else:
                require(pair["status"] == "available", "Display pair status")
                witness = trial["similarity_constraints"]["witness"]
                reference, candidate = (
                    witness["reference_pose"],
                    trial["selective_pose_evidence"]["hit_candidate"],
                )
                require(
                    reference["smiles"] == hit
                    and candidate["smiles"] == row["smiles"]
                    and witness["candidate_pose_sha256"] == candidate["sha256"],
                    "Display pair identities/witness",
                )
                require(
                    pair["reference_pose_sha256"] == reference["sha256"]
                    and pair["candidate_pose_sha256"] == candidate["sha256"]
                    and pair["witness_sha256"] == binding(witness),
                    "Display same-seed pose bindings",
                )
                checked_graph(pair["reference_graph"], reference)
                checked_graph(pair["candidate_graph"], candidate)
                view["coordinates"] = {
                    "reference": reference["coordinates"],
                    "candidate": candidate["coordinates"],
                }
                view["metrics"] = {
                    "shape": trial["hit"]["shape"],
                    "pharmacophore": trial["hit"]["pharmacophore"],
                    "polar_recall": witness["polar"]["recall"],
                }
            seed_views.append(view)
        result.append(
            {
                "status": "available",
                "depictions": entry["depictions"],
                "seeds": seed_views,
                "representative_seed": entry["representative_seed"],
            }
        )
    return result
