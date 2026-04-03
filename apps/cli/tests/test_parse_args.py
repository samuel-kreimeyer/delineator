from delineator_cli.main import parse_args


def test_parse_args_accepts_basic_inputs() -> None:
    parsed = parse_args(["input.xml", "-o", "output", "--use-cogo-points"])
    assert parsed.input_file.name == "input.xml"
    assert parsed.output_dir.name == "output"
    assert parsed.use_cogo_points is True
