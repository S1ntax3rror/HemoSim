import pytest

from HemoAnalysisLibrary import fe


def test_fe_returns():
    assert fe.get_fe_distances() is None


if __name__ == "__main__":
    pytest.main()
