import json
import tomli
import os

# import yaml


def toml_to_json(toml_file):
    with open(toml_file, mode="rb") as fileObj:
        content = tomli.load(fileObj)
        with open(f"{os.path.splitext(toml_file)[0]}.json", "w", encoding="utf-8") as f:
            json.dump(content, f, ensure_ascii=False, indent=4)

        # to yaml if  I change my mind
        # with open(f"{os.path.splitext(toml_file)[0]}.yaml", "w", encoding="utf-8") as f:
        #     config_file = yaml.dump(content, f)
        #     print(config_file)

    return


toml_to_json("some_toml_file")
