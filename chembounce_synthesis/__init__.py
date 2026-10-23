"""Synthesis-aware design, independent of the optional SynLlama runtime."""


def recommend(hit, designs, config, output_dir, *, deadline_monotonic=None):
    """Lazily request similar, scaffold-changed, computationally supported products."""
    from .workflow import recommend as implementation

    return implementation(hit, designs, config, output_dir, deadline_monotonic=deadline_monotonic)


def recommend_routes(hit, config, output_dir, *, deadline_monotonic=None):
    """Opt into bounded catalog/rule route-first generation without a design/model provider."""
    from .workflow import recommend_routes as implementation

    return implementation(hit, config, output_dir, deadline_monotonic=deadline_monotonic)
