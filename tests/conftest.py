# Test configuration. Imported by pytest automatically.

def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "slow: long-running tests",
    )
