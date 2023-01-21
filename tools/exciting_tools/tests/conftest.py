"""
File containing pytest fixtures. They can be seen from all other files and subdirectories.
That's because of the special name of the file.
Needed for testing environment variables.
"""
import pytest


@pytest.fixture
def mock_env_jobflow(monkeypatch):
    monkeypatch.setenv("USE_JOBFLOW", "true")


@pytest.fixture
def mock_env_jobflow_missing(monkeypatch):
    monkeypatch.delenv("USE_JOBFLOW", raising=False)
