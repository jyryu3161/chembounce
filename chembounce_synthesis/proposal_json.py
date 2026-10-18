"""Opt-in, auditable JSON dialect for literal SMILES stereo bond backslashes.

This interprets a narrow class of malformed model text. It does not recover the
model's intent, and strict JSON remains the default provider contract.
"""

import hashlib
import json
import math
import time

POLICY = "smiles-bond-backslash-v1"
_SMILES_ATOM_START = frozenset("BCNOPSFIbcnops*[")
_JSON_ESCAPES = frozenset('"\\/bfnrt')
_MAX_REPAIRS = 256
_MAX_OUTPUT = 500_000


def _sha(text):
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _json_text(value):
    _require_scalars(value)
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False)


def _reject_constant(value):
    raise ValueError(f"Nonfinite JSON extension: {value}")


def _require_scalars(value):
    if isinstance(value, str):
        if any(0xD800 <= ord(ch) <= 0xDFFF for ch in value):
            raise ValueError("Non-scalar Unicode in JSON string")
    elif isinstance(value, list):
        for item in value:
            _require_scalars(item)
    elif isinstance(value, dict):
        for key, item in value.items():
            _require_scalars(key)
            _require_scalars(item)


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def _is_smiles_path(path):
    return (
        (len(path) == 3 and path[0] == "reactions" and isinstance(path[1], int) and path[2] == "product")
        or (
            len(path) == 4
            and path[0] == "reactions"
            and isinstance(path[1], int)
            and path[2] == "reactants"
            and isinstance(path[3], int)
        )
        or (len(path) == 2 and path[0] == "building_blocks" and isinstance(path[1], int))
    )


class _Scanner:
    def __init__(self, raw, deadline_monotonic=None):
        if not isinstance(raw, str) or len(raw) > 250_000:
            raise ValueError("Invalid or oversized model response")
        self.raw = raw
        self.start = len(raw) - len(raw.lstrip())
        self.end = len(raw.rstrip())
        if raw.startswith("```json", self.start) and raw[self.end - 3 : self.end] == "```":
            self.start += 7
            self.end -= 3
        elif raw.startswith("```", self.start) and raw[self.end - 3 : self.end] == "```":
            self.start += 3
            self.end -= 3
        self.pos = self.start
        self.repairs = []
        self.parts = []
        self.output_size = 0
        self.ambiguous = []
        self.deadline_monotonic = deadline_monotonic

    def _budget(self):
        if self.deadline_monotonic is not None and time.monotonic() >= self.deadline_monotonic:
            raise ValueError("Proposal decode deadline reached")
        if len(self.repairs) > _MAX_REPAIRS or self.output_size > _MAX_OUTPUT:
            raise ValueError("Proposal decode size limit exceeded")

    def _append(self, value):
        self.parts.append(value)
        self.output_size += len(value)
        if self.output_size > _MAX_OUTPUT:
            raise ValueError("Proposal derived JSON size limit exceeded")

    def _space(self):
        while self.pos < self.end and self.raw[self.pos] in " \t\r\n":
            self.pos += 1

    def _string(self, path, *, key=False):
        if self.pos >= self.end or self.raw[self.pos] != '"':
            raise ValueError("Expected JSON string")
        start = self.pos
        self.pos += 1
        escaped = []
        while self.pos < self.end:
            if self.pos % 1024 == 0:
                self._budget()
            ch = self.raw[self.pos]
            if ch == '"':
                self.pos += 1
                token = '"' + "".join(escaped) + '"'
                value = json.loads(token)
                self._append(token)
                return value
            if ord(ch) < 32:
                raise ValueError("Unescaped control character in JSON string")
            if ch != "\\":
                escaped.append(ch)
                self.pos += 1
                continue
            if self.pos + 1 >= self.end:
                raise ValueError("Unterminated JSON escape")
            follow = self.raw[self.pos + 1]
            if follow in _JSON_ESCAPES:
                if _is_smiles_path(path) and follow != "\\":
                    self.ambiguous.append({"path": list(path), "character_offset": self.pos})
                escaped.append(self.raw[self.pos : self.pos + 2])
                self.pos += 2
                continue
            if follow == "u":
                if self.pos + 6 > self.end or any(
                    digit not in "0123456789abcdefABCDEF" for digit in self.raw[self.pos + 2 : self.pos + 6]
                ):
                    raise ValueError("Invalid Unicode escape in JSON string")
                escaped.append(self.raw[self.pos : self.pos + 6])
                self.pos += 6
                continue
            if key or not _is_smiles_path(path) or follow not in _SMILES_ATOM_START:
                raise ValueError("Invalid JSON escape outside allowed SMILES bond position")
            if self.pos > start + 1 and self.raw[self.pos - 1] == "\\":
                raise ValueError("Ambiguous consecutive JSON backslashes")
            self.repairs.append(
                {
                    "path": list(path),
                    "character_offset": self.pos,
                    "utf8_byte_offset": len(self.raw[: self.pos].encode("utf-8")),
                    "following_character": follow,
                    "original_json_text": "\\" + follow,
                    "replacement_json_text": "\\\\" + follow,
                }
            )
            if len(self.repairs) > _MAX_REPAIRS:
                raise ValueError("Proposal repair count limit exceeded")
            escaped.append("\\\\" + follow)
            self.pos += 2
        raise ValueError("Unterminated JSON string")

    def _value(self, path=(), depth=0):
        self._budget()
        if depth > 128:
            raise ValueError("JSON nesting limit exceeded")
        self._space()
        if self.pos >= self.end:
            raise ValueError("Incomplete JSON value")
        ch = self.raw[self.pos]
        if ch == '"':
            self._string(path)
        elif ch == "{":
            self._append(ch)
            self.pos += 1
            self._space()
            if self.pos < self.end and self.raw[self.pos] == "}":
                self._append("}")
                self.pos += 1
                return
            while True:
                self._space()
                key = self._string(path, key=True)
                self._space()
                if self.pos >= self.end or self.raw[self.pos] != ":":
                    raise ValueError("Expected JSON colon")
                self._append(":")
                self.pos += 1
                self._value((*path, key), depth + 1)
                self._space()
                if self.pos >= self.end:
                    raise ValueError("Incomplete JSON object")
                delimiter = self.raw[self.pos]
                if delimiter not in ",}":
                    raise ValueError("Invalid JSON object delimiter")
                self._append(delimiter)
                self.pos += 1
                if delimiter == "}":
                    return
        elif ch == "[":
            self._append(ch)
            self.pos += 1
            self._space()
            if self.pos < self.end and self.raw[self.pos] == "]":
                self._append("]")
                self.pos += 1
                return
            index = 0
            while True:
                self._value((*path, index), depth + 1)
                index += 1
                self._space()
                if self.pos >= self.end:
                    raise ValueError("Incomplete JSON array")
                delimiter = self.raw[self.pos]
                if delimiter not in ",]":
                    raise ValueError("Invalid JSON array delimiter")
                self._append(delimiter)
                self.pos += 1
                if delimiter == "]":
                    return
        else:
            start = self.pos
            while self.pos < self.end and self.raw[self.pos] not in ",]} \t\r\n":
                self.pos += 1
            if start == self.pos:
                raise ValueError("Invalid JSON value")
            token = self.raw[start : self.pos]
            value = json.loads(token, parse_constant=_reject_constant)
            if isinstance(value, (dict, list, str)):
                raise ValueError("Invalid JSON scalar")
            self._append(token)

    def decode(self):
        self._value()
        self._space()
        if self.pos != self.end:
            raise ValueError("Extra text after JSON value")
        if not self.repairs:
            raise ValueError("No compatible SMILES escape found")
        if self.ambiguous:
            raise ValueError("Ambiguous recognized JSON escape in SMILES value")
        repaired = "".join(self.parts)
        if len(repaired) > _MAX_OUTPUT:
            raise ValueError("Proposal derived JSON size limit exceeded")
        self._budget()
        data = json.loads(repaired, object_pairs_hook=_unique_object, parse_constant=_reject_constant)
        _require_scalars(data)
        if (
            not isinstance(data, dict)
            or not isinstance(data.get("reactions"), list)
            or not isinstance(data.get("building_blocks"), list)
        ):
            raise ValueError("Model response requires reactions and building_blocks lists")
        if not all(isinstance(x, str) for x in data["building_blocks"]) or not all(
            isinstance(row, dict)
            and isinstance(row.get("product"), str)
            and isinstance(row.get("reactants"), list)
            and all(isinstance(x, str) for x in row["reactants"])
            for row in data["reactions"]
        ):
            raise ValueError("Incompatible proposal SMILES field structure")
        canonical = json.dumps(
            data, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False
        )
        self._budget()
        return data, {
            "policy": POLICY,
            "status": "compatibility_decoded",
            "strict_parser_outcome": "malformed_json",
            "raw_utf8_sha256": _sha(self.raw),
            "derived_json_text": repaired,
            "derived_json_sha256": _sha(repaired),
            "parsed_value_sha256": _sha(canonical),
            "repairs": self.repairs,
        }


def _deadline_row(raw):
    return {
        "response": raw,
        "raw_response": raw,
        "compatibility": {
            "policy": POLICY,
            "status": "not_evaluated_deadline",
            "strict_parser_outcome": "not_evaluated_deadline",
            "raw_utf8_sha256": _sha(raw),
            "reason": "Proposal decode deadline reached",
        },
    }


def decode_response(raw, *, deadline_monotonic=None):
    """Return a raw-preserving provider row for a versioned opt-in dialect."""
    if not isinstance(raw, str):
        raise ValueError("Model response must be text")
    if any(0xD800 <= ord(ch) <= 0xDFFF for ch in raw):
        raise ValueError("Non-scalar Unicode in raw model response")
    if deadline_monotonic is not None and (
        isinstance(deadline_monotonic, bool)
        or not isinstance(deadline_monotonic, (int, float))
        or not math.isfinite(deadline_monotonic)
    ):
        raise ValueError("Invalid proposal decode deadline")
    if deadline_monotonic is not None and time.monotonic() >= deadline_monotonic:
        return _deadline_row(raw)
    try:
        from .routes import parse_response

        parse_response(raw)
    except json.JSONDecodeError as strict_error:
        try:
            data, audit = _Scanner(raw, deadline_monotonic).decode()
        except (ValueError, json.JSONDecodeError, RecursionError) as exc:
            if str(exc) == "Proposal decode deadline reached":
                return _deadline_row(raw)
            return {
                "response": raw,
                "raw_response": raw,
                "compatibility": {
                    "policy": POLICY,
                    "status": "rejected",
                    "strict_parser_outcome": "malformed_json",
                    "strict_error": str(strict_error),
                    "raw_utf8_sha256": _sha(raw),
                    "reason": str(exc),
                },
            }
        audit["strict_error"] = str(strict_error)
        return {"response": data, "raw_response": raw, "compatibility": audit}
    except (ValueError, TypeError) as strict_error:
        return {
            "response": raw,
            "raw_response": raw,
            "compatibility": {
                "policy": POLICY,
                "status": "rejected",
                "strict_parser_outcome": "invalid_schema",
                "strict_error": str(strict_error),
                "raw_utf8_sha256": _sha(raw),
                "reason": "Strict JSON has invalid proposal schema",
            },
        }
    if deadline_monotonic is not None and time.monotonic() >= deadline_monotonic:
        return _deadline_row(raw)
    return {
        "response": raw,
        "raw_response": raw,
        "compatibility": {
            "policy": POLICY,
            "status": "strict_valid",
            "strict_parser_outcome": "valid",
            "raw_utf8_sha256": _sha(raw),
        },
    }


def bind_row(raw, sampling, prompt_sha256, decoder_source_sha256, *, design, deadline_monotonic=None):
    """Create a proposal row whose audit binds raw text, sampling and policy source."""
    if not isinstance(design, str) or not design:
        raise ValueError("Compatibility proposal requires a design identity")
    if (
        not isinstance(sampling, dict)
        or type(sampling.get("sample_index")) is not int
        or sampling["sample_index"] < 0
    ):
        raise ValueError("Compatibility proposal requires a sample index")
    row = decode_response(raw, deadline_monotonic=deadline_monotonic)
    row["sampling"] = sampling
    sampling_sha256 = _sha(_json_text(sampling))
    if deadline_monotonic is not None and time.monotonic() >= deadline_monotonic:
        row = _deadline_row(raw)
        row["sampling"] = sampling
    row["compatibility"].update(
        {
            "prompt_sha256": prompt_sha256,
            "decoder_source_sha256": decoder_source_sha256,
            "sampling_sha256": sampling_sha256,
            "design_smiles": design,
            "sample_index": sampling["sample_index"],
        }
    )
    return row


def verify_saved_row(row, provider, *, design):
    """Rederive an opt-in archive row; reject missing or detached audit fields."""
    if not isinstance(row, dict) or not all(
        key in row for key in ("response", "raw_response", "compatibility", "sampling")
    ):
        raise ValueError("Saved compatibility proposal is missing raw response or audit")
    if not isinstance(row["raw_response"], str) or not isinstance(row["sampling"], dict):
        raise ValueError("Saved compatibility proposal has invalid raw or sampling type")
    if not isinstance(design, str) or not design:
        raise ValueError("Saved compatibility proposal requires a design identity")
    if type(row["sampling"].get("sample_index")) is not int or row["sampling"]["sample_index"] < 0:
        raise ValueError("Saved compatibility proposal requires a sample index")
    if (
        isinstance(row["compatibility"], dict)
        and row["compatibility"].get("status") == "not_evaluated_deadline"
    ):
        expected = _deadline_row(row["raw_response"])
        expected["sampling"] = row["sampling"]
        expected["compatibility"].update(
            {
                "prompt_sha256": provider["prompt_sha256"],
                "decoder_source_sha256": provider["response_decoder_source_sha256"],
                "sampling_sha256": _sha(_json_text(row["sampling"])),
                "design_smiles": design,
                "sample_index": row["sampling"].get("sample_index"),
            }
        )
    else:
        expected = bind_row(
            row["raw_response"],
            row["sampling"],
            provider["prompt_sha256"],
            provider["response_decoder_source_sha256"],
            design=design,
        )
    try:
        actual_response = _json_text(row["response"])
        expected_response = _json_text(expected["response"])
        actual_audit = _json_text(row["compatibility"])
        expected_audit = _json_text(expected["compatibility"])
    except (TypeError, ValueError, OverflowError, RecursionError) as exc:
        raise ValueError("Saved compatibility proposal has invalid JSON values") from exc
    if actual_response != expected_response or actual_audit != expected_audit:
        raise ValueError("Saved compatibility proposal differs from raw rederivation")
    return expected
