"""Synthesis-aware design, independent of the optional SynLlama runtime."""


def recommend(hit, designs, config, output_dir, *, deadline_monotonic=None):
    """Lazily request similar, scaffold-changed, computationally supported products."""
    from .workflow import recommend as implementation

    return implementation(hit, designs, config, output_dir, deadline_monotonic=deadline_monotonic)
