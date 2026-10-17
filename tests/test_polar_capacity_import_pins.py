"""Imported-runtime pinning with synthetic module objects; no Torch/native imports."""

from importlib.machinery import ModuleSpec
from importlib.util import spec_from_file_location
from pathlib import Path
from types import ModuleType, SimpleNamespace

import pytest

from benchmarks.research import polar_capacity_prepare as p


def file_module(tmp_path, name, filename=None):
    path = tmp_path / (filename or (name.replace(".", "_") + ".py"))
    path.write_text("# synthetic implementation, never imported\n")
    module = ModuleType(name)
    module.__file__ = str(path)
    module.__spec__ = spec_from_file_location(name, path)
    module.__loader__ = module.__spec__.loader
    return module


def torch_namespaces(tmp_path):
    """Mirror installed _Ops/_Classes and child class metadata without their code."""
    ops = file_module(tmp_path, "torch._ops")
    classes = file_module(tmp_path, "torch._classes")

    def dynamic_lookup(self, name):
        raise AssertionError("Pinning must never invoke dynamic __getattr__: " + name)

    definitions = [
        (ops, "_Ops", "_ops.py", "torch.ops", "ops"),
        (classes, "_Classes", "_classes.py", "torch.classes", "classes"),
        (ops, "_OpNamespace", "torch.ops", "torch.ops.synthetic", "synthetic"),
        (ops, "_HigherOrderNamespace", "torch.ops", "torch.ops.higher_order", "higher_order"),
        (classes, "_ClassNamespace", None, "torch.classessynthetic", "synthetic"),
    ]
    modules = {"torch._ops": ops, "torch._classes": classes}
    for owner, cls_name, label, declared, attr in definitions:
        members = {"__module__": owner.__name__, "__getattr__": dynamic_lookup}
        if label is not None:
            members["__file__"] = label
        cls = type(cls_name, (ModuleType,), members)
        setattr(owner, cls_name, cls)
        namespace = cls(declared)
        if cls_name in {"_Ops", "_Classes"}:
            setattr(owner, attr, namespace)
            modules[declared] = namespace
        else:
            parent = "torch.classes" if cls_name == "_ClassNamespace" else "torch.ops"
            namespace.name = attr
            setattr(modules[parent], attr, namespace)
            modules[parent + "." + attr] = namespace
    return modules


def pin(monkeypatch, modules):
    # Replace only this administrative module's sys reference, not global sys.modules.
    monkeypatch.setattr(p, "sys", SimpleNamespace(modules=modules))
    return p.imported_file_pins()


def test_virtual_namespaces_pin_implementation_not_cwd_labels(tmp_path, monkeypatch):
    modules = torch_namespaces(tmp_path)
    monkeypatch.chdir(tmp_path)
    # Misleading local files must not become the virtual namespace's authority.
    for name in ("_ops.py", "_classes.py", "torch.ops"):
        (tmp_path / name).write_text("wrong local implementation\n")
    expected = {
        m.__file__: p.d.digest(m.__file__)
        for n, m in modules.items()
        if n in {"torch._ops", "torch._classes"}
    }
    assert pin(monkeypatch, modules) == expected
    before = {n: dict(vars(m)) for n, m in modules.items()}
    assert pin(monkeypatch, modules) == expected
    assert {n: dict(vars(m)) for n, m in modules.items()} == before


@pytest.mark.parametrize(
    "attack",
    [
        "class",
        "name",
        "singleton",
        "spec",
        "loader",
        "own_file",
        "label",
        "child_parent",
        "child_name",
        "missing_implementation",
    ],
)
def test_virtual_namespace_spoof_or_mismatched_provenance_fails(tmp_path, monkeypatch, attack):
    modules = torch_namespaces(tmp_path)
    namespace = modules["torch.ops"]
    if attack == "class":
        modules["torch.ops"] = ModuleType("torch.ops")
    elif attack == "name":
        namespace.__name__ = "torch.other"
    elif attack == "singleton":
        modules["torch._ops"].ops = ModuleType("torch.ops")
    elif attack == "spec":
        namespace.__spec__ = ModuleSpec("torch.ops", loader=None, origin="_ops.py")
    elif attack == "loader":
        namespace.__loader__ = object()
    elif attack == "own_file":
        namespace.__file__ = "_ops.py"
    elif attack == "label":
        type(namespace).__file__ = "other.py"
    elif attack == "child_parent":
        namespace.synthetic = ModuleType("torch.ops.synthetic")
    elif attack == "child_name":
        modules["torch.ops.synthetic"].name = "wrong"
    else:
        modules.pop("torch._ops")
    with pytest.raises(ValueError, match="Torch namespace"):
        pin(monkeypatch, modules)


@pytest.mark.parametrize(
    "attack",
    [
        "missing_file",
        "missing_origin",
        "different_origin",
        "absent_origin",
        "unrecognized_spec",
        "spec_name",
        "module_name",
        "loader_identity",
    ],
)
def test_real_implementation_missing_file_or_spec_is_not_exempted(tmp_path, monkeypatch, attack):
    modules = torch_namespaces(tmp_path)
    implementation = modules["torch._ops"]
    if attack == "missing_file":
        Path(implementation.__file__).unlink()
    elif attack == "missing_origin":
        implementation.__file__ = str(tmp_path / "missing.py")
        implementation.__spec__ = spec_from_file_location("torch._ops", implementation.__file__)
    elif attack == "different_origin":
        implementation.__spec__ = spec_from_file_location("torch._ops", tmp_path / "different.py")
    elif attack == "absent_origin":
        implementation.__spec__ = None
    elif attack == "spec_name":
        implementation.__spec__.name = "different"
    elif attack == "module_name":
        implementation.__name__ = "different"
    elif attack == "loader_identity":
        implementation.__loader__ = object()
    else:
        implementation.__spec__ = SimpleNamespace(origin=implementation.__file__)
    with pytest.raises(ValueError, match="Imported|Torch namespace implementation"):
        pin(monkeypatch, modules)


def test_real_files_relative_symlinks_spec_only_and_fileless_modules(tmp_path, monkeypatch):
    module = file_module(tmp_path, "real", "real.py")
    monkeypatch.chdir(tmp_path)
    module.__file__ = "real.py"
    alias = tmp_path / "alias.py"
    alias.symlink_to(tmp_path / "real.py")
    other = ModuleType("alias")
    other.__file__ = str(alias)
    other.__spec__ = spec_from_file_location("alias", alias)
    origin_only = ModuleType("origin_only")
    origin_only.__spec__ = spec_from_file_location("origin_only", tmp_path / "real.py")
    builtin = ModuleType("builtin")
    builtin.__spec__ = ModuleSpec("builtin", loader=None, origin="built-in")
    namespace = ModuleType("namespace")
    namespace.__spec__ = ModuleSpec("namespace", loader=None, is_package=True)
    modules = {m.__name__: m for m in [module, other, origin_only, builtin, namespace]}
    assert pin(monkeypatch, modules) == {str(tmp_path / "real.py"): p.d.digest(tmp_path / "real.py")}


@pytest.mark.parametrize("declared", ["file", "spec", "nonlocation_spec"])
def test_missing_real_path_still_fails(tmp_path, monkeypatch, declared):
    module = ModuleType("missing")
    missing = tmp_path / "does-not-exist.py"
    if declared == "file":
        module.__file__ = str(missing)
    elif declared == "spec":
        module.__spec__ = spec_from_file_location("missing", missing)
    else:
        module.__spec__ = ModuleSpec("missing", loader=None, origin=str(missing))
    with pytest.raises(ValueError, match="Imported file cannot be pinned"):
        pin(monkeypatch, {"missing": module})


def test_unreadable_real_file_is_not_ignored(tmp_path, monkeypatch):
    module = file_module(tmp_path, "real")

    def unreadable(path):
        raise PermissionError("synthetic unreadable file")

    monkeypatch.setattr(p.d, "digest", unreadable)
    with pytest.raises(PermissionError, match="unreadable"):
        pin(monkeypatch, {"real": module})


def test_unknown_inherited_file_claim_fails_without_dynamic_getattr(tmp_path, monkeypatch):
    class Unknown(ModuleType):
        __file__ = "virtual.py"

        def __getattr__(self, name):
            raise AssertionError("dynamic lookup")

    with pytest.raises(ValueError, match="Unknown virtual imported file"):
        pin(monkeypatch, {"unknown": Unknown("unknown")})
