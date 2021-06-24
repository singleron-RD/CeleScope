def pytest_addoption(parser):
    parser.addoption("--assays", action="store", default="")
    parser.addoption("--test_dir", help='test_dir', action="store", default="/SGRNJ03/randd/user/zhouyiqi/multi_tests/")


def pytest_generate_tests(metafunc):
    # This is called for every test. Only get/set command line arguments
    # if the argument is specified in the list of test "fixturenames".
    assays_value = metafunc.config.option.assays
    test_dir_value = metafunc.config.option.test_dir
    if 'assays' in metafunc.fixturenames and assays_value is not None:
        metafunc.parametrize("assays", [assays_value])
    if 'test_dir' in metafunc.fixturenames and test_dir_value is not None:
        metafunc.parametrize("test_dir", [test_dir_value])
