"""Observe official decode inputs without altering its return, errors or RNG."""

from contextlib import contextmanager


@contextmanager
def observe_decode(vocabulary):
    original = vocabulary.decode
    had_instance_attribute = "decode" in vars(vocabulary)
    previous = vars(vocabulary).get("decode")
    captured = []

    def decode(sequence, *args, **kwargs):
        # Snapshot only plain integer IDs. The original call and its exception pass through.
        result = original(sequence, *args, **kwargs)
        captured.append([int(x) for x in sequence])
        return result

    vocabulary.decode = decode
    try:
        yield captured
    finally:
        if had_instance_attribute:
            vocabulary.decode = previous
        else:
            del vocabulary.decode


def bind_outputs(
    captured, source_rows, inputs, outputs, nlls, *, eos, bos, max_length, tokenizer, vocabulary
):
    count = len(source_rows)
    if not (len(captured) == 2 * count and len(inputs) == len(outputs) == len(nlls) == count):
        raise ValueError("Official sample/decode row count changed")
    if captured[:count] != source_rows:
        raise ValueError("Official input decode order changed")
    records = []
    for index, (tokens, text, nll) in enumerate(zip(captured[count:], outputs, nlls)):
        if not tokens or tokens[0] != bos or len(tokens) > max_length:
            raise ValueError("Unexpected official generated token shape")
        terminated = eos in tokens
        end = tokens.index(eos) if terminated else None
        if terminated and any(token != 0 for token in tokens[end + 1 :]):
            raise ValueError("Unexpected non-padding tokens after EOS")
        if not terminated and len(tokens) != max_length:
            raise ValueError("Unterminated output did not reach official length limit")
        if tokenizer.untokenize(vocabulary.decode(tokens)) != text:
            raise ValueError("Generated row/token binding mismatch")
        records.append(
            {
                "batch_row": index,
                "smiles": text,
                "native_nll": float(nll),
                "token_ids": tokens,
                "terminated": terminated,
                "truncated": not terminated,
                "eos_position": end,
                "padded_length": len(tokens),
            }
        )
    return records


def validate_tokens(smiles, tokenizer, vocabulary, max_length):
    tokens = tokenizer.tokenize(smiles)
    unknown = sorted({token for token in tokens if token not in vocabulary})
    if unknown:
        raise ValueError(f"Unsupported input tokens (official encode would drop/map them): {unknown}")
    if len(tokens) > max_length:
        raise ValueError("Input exceeds official encoder positional capacity including BOS/EOS")
    return tokens
