import pytest

from delayed_entry.adapters import genie_bpc
from helpers import FIXTURE, crc_settings


@pytest.fixture(scope="session")
def fixture_result():
    """The made-up BPC fixture read with the CRC settings."""
    return genie_bpc.read(FIXTURE, None, crc_settings())
