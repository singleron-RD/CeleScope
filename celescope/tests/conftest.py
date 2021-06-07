def pytest_addoption(parser):
    parser.addoption("--assays", action="store", default="")


def pytest_generate_tests(metafunc):
    # This is called for every test. Only get/set command line arguments
    # if the argument is specified in the list of test "fixturenames".
    option_value = metafunc.config.option.assays
    if 'assays' in metafunc.fixturenames and option_value is not None:
        metafunc.parametrize("assays", [option_value])