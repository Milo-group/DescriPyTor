# Path helpers for CS3 recreation tests. Imported by pytest automatically.

def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "slow: full CS3 GOAT table recreation (thousands of conformers)",
    )
