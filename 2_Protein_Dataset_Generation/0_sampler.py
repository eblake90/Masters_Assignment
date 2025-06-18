#!/usr/bin/env python3

import requests
import json


def main():
    accession = "P22626"
    data = requests.get(f"https://rest.uniprot.org/uniprotkb/{accession}.json").json()

    with open(f"data/0_{accession}_sample.json", 'w') as f:
        json.dump(data, f, indent=2)

    print("Done")


if __name__ == "__main__":
    main()