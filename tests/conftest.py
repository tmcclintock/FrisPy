def pytest_addoption(parser):
    """Options for pytest.

    run-slow: used to mark tests that are slow to run, and should
        only be run on CI
    """
    parser.addoption(
        "--run-slow", action="store_true", default=False, help="Run slow tests."
    )
