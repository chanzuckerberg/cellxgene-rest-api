import json


def parse_schema(filename):
	with open(filename) as fh:
		schema = json.load(fh)
	return schema


if __name__ == "__main__":
	parse_schema("data/test_data_schema.json")
