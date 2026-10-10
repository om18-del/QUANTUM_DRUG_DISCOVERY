from __future__ import annotations

from app.config import Settings, project_root


def test_project_root_is_repo_root():
    settings = Settings()
    repo_root = project_root()
    settings_paths = [
        settings.cache_dir,
        settings.job_store_path,
        settings.h2_hamiltonian_dir,
        settings.complexes_manifest_path,
        settings.models_manifest_path,
    ]
    for path in settings_paths:
        resolved = path.resolve()
        assert repo_root in resolved.parents or resolved == repo_root, (
            f"Path {path} is outside repo root {repo_root}"
        )


def test_env_file_points_at_project_root():
    settings = Settings()
    expected_env = project_root() / ".env"
    assert settings.model_config["env_file"] == str(expected_env)


def test_project_root_is_independent_of_cwd(tmp_path, monkeypatch):
    """project_root() must resolve from this file's location, not from os.getcwd().

    Otherwise every data path moves when the server is started from a different
    directory (uvicorn --reload, pytest from a subdirectory, a service manager).
    """
    expected = project_root()
    assert expected.is_absolute()
    assert (expected / "app").is_dir()
    assert (expected / "SPEC.md").is_file()

    monkeypatch.chdir(tmp_path)

    assert project_root() == expected
    assert Settings().h2_hamiltonian_dir == expected / "app" / "data" / "hamiltonians"
    assert Settings().complexes_manifest_path == (
        expected / "app" / "data" / "curated_complexes" / "manifest.json"
    )
    assert Settings().model_config["env_file"] == str(expected / ".env")
